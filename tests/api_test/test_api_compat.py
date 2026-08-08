"""
Backwards-compatibility guard for the public be_scan API.

be_scan is published and other labs have code calling it. Existing calls must
keep working across versions. This module compares the live signatures against
the v2.0.0 snapshot in api_baseline_v2.json and enforces five rules:

  1. Every baseline function still exists in the same module namespace.
  2. Every baseline parameter still exists, spelled identically.
  3. Every baseline default is unchanged.
  4. Baseline parameters keep their relative order, and anything that was
     positionally reachable is still reachable at the same index.
  5. New parameters must have a default, and must be appended after the
     existing ones (or be keyword-only).

Together these mean a new feature may add `guide_len=20` to count_reads, but
may not rename df_filepath, reorder comparisons, or flip save_df to False.

If a test here fails, the default assumption is that the source change is
wrong, not that the baseline is stale. See generate_baseline.py.
"""

import importlib
import inspect
import json
from pathlib import Path

import pytest

from tests.api_test.generate_baseline import encode_default

BASELINE_PATH = Path(__file__).resolve().parent / "api_baseline_v2.json"
BASELINE = json.loads(BASELINE_PATH.read_text())

POSITIONAL_KINDS = {"POSITIONAL_ONLY", "POSITIONAL_OR_KEYWORD"}


def _live_params(modname, fnname):
    """Return the current parameter list for a function, in baseline format."""
    module = importlib.import_module(modname)
    fn = getattr(module, fnname)
    return [
        {
            "name": p.name,
            "kind": str(p.kind),
            "default": encode_default(p.default),
        }
        for p in inspect.signature(fn).parameters.values()
    ]


def _all_functions():
    """Yield (modname, fnname) for every function in the baseline."""
    for modname, fns in sorted(BASELINE.items()):
        for fnname in sorted(fns):
            yield modname, fnname


ALL_FUNCTIONS = list(_all_functions())
FUNCTION_IDS = [f"{m.split('.')[-1]}.{f}" for m, f in ALL_FUNCTIONS]


def test_baseline_is_populated():
    """Guard against an empty or truncated baseline silently passing everything."""
    n_fns = sum(len(v) for v in BASELINE.values())
    n_params = sum(len(f["params"]) for v in BASELINE.values() for f in v.values())
    assert n_fns == 40, f"baseline has {n_fns} functions, expected 40"
    assert n_params == 468, f"baseline has {n_params} parameters, expected 468"


# --- Rule 1: functions still exist ---------------------------------------- #


@pytest.mark.parametrize("modname,fnname", ALL_FUNCTIONS, ids=FUNCTION_IDS)
def test_function_still_exported(modname, fnname):
    module = importlib.import_module(modname)
    assert hasattr(module, fnname), (
        f"{modname}.{fnname} existed in v2.0.0 and is now missing. "
        f"Removing or moving a public function breaks existing user code."
    )
    assert inspect.isfunction(getattr(module, fnname)), (
        f"{modname}.{fnname} is no longer a function."
    )


# --- Rules 2 and 3: parameters and their defaults -------------------------- #


@pytest.mark.parametrize("modname,fnname", ALL_FUNCTIONS, ids=FUNCTION_IDS)
def test_parameters_preserved(modname, fnname):
    expected = BASELINE[modname][fnname]["params"]
    actual = _live_params(modname, fnname)
    actual_by_name = {p["name"]: p for p in actual}

    for param in expected:
        name = param["name"]
        assert name in actual_by_name, (
            f"{modname}.{fnname}: parameter '{name}' existed in v2.0.0 and is "
            f"now missing. Renaming or removing a parameter breaks callers "
            f"that pass it by keyword."
        )
        assert actual_by_name[name]["default"] == param["default"], (
            f"{modname}.{fnname}: default for '{name}' changed from "
            f"{param['default']!r} to {actual_by_name[name]['default']!r}. "
            f"Callers relying on the old default would silently get different "
            f"results."
        )


# --- Rule 4: ordering and positional reachability -------------------------- #


@pytest.mark.parametrize("modname,fnname", ALL_FUNCTIONS, ids=FUNCTION_IDS)
def test_parameter_order_preserved(modname, fnname):
    expected = BASELINE[modname][fnname]["params"]
    actual = _live_params(modname, fnname)

    expected_names = [p["name"] for p in expected]
    actual_names = [p["name"] for p in actual]

    # The baseline names must appear in the live signature in the same order.
    surviving = [n for n in actual_names if n in set(expected_names)]
    assert surviving == expected_names, (
        f"{modname}.{fnname}: parameter order changed.\n"
        f"  v2.0.0 order: {expected_names}\n"
        f"  current order: {surviving}\n"
        f"Reordering breaks callers that pass arguments positionally."
    )


@pytest.mark.parametrize("modname,fnname", ALL_FUNCTIONS, ids=FUNCTION_IDS)
def test_positional_arguments_still_reachable(modname, fnname):
    """A param passed positionally in v2.0.0 must still land in the same slot."""
    expected = BASELINE[modname][fnname]["params"]
    actual = _live_params(modname, fnname)

    expected_positional = [p["name"] for p in expected if p["kind"] in POSITIONAL_KINDS]
    actual_positional = [p["name"] for p in actual if p["kind"] in POSITIONAL_KINDS]

    for index, name in enumerate(expected_positional):
        assert index < len(actual_positional), (
            f"{modname}.{fnname}: '{name}' was positional argument {index} in "
            f"v2.0.0 but the function now accepts only "
            f"{len(actual_positional)} positional arguments."
        )
        assert actual_positional[index] == name, (
            f"{modname}.{fnname}: positional slot {index} was '{name}' in "
            f"v2.0.0 but is now '{actual_positional[index]}'. "
            f"A caller passing arguments positionally would pass the wrong "
            f"value to the wrong parameter."
        )


# --- Rule 5: new parameters must be safe additions ------------------------- #


@pytest.mark.parametrize("modname,fnname", ALL_FUNCTIONS, ids=FUNCTION_IDS)
def test_new_parameters_are_additive(modname, fnname):
    expected = BASELINE[modname][fnname]["params"]
    actual = _live_params(modname, fnname)

    expected_names = {p["name"] for p in expected}
    new_params = [p for p in actual if p["name"] not in expected_names]

    for param in new_params:
        if param["kind"] in ("VAR_POSITIONAL", "VAR_KEYWORD"):
            continue  # *args / **kwargs never break an existing call
        assert param["default"] != "<required>", (
            f"{modname}.{fnname}: new parameter '{param['name']}' has no "
            f"default. Every existing call would immediately fail with a "
            f"TypeError. Give it a default that preserves current behavior."
        )

    # New parameters must not be interleaved among the v2.0.0 ones, or they
    # would shift positional slots. Keyword-only params are exempt: they are
    # unreachable positionally by construction.
    actual_names = [p["name"] for p in actual]
    kw_only = {p["name"] for p in actual if p["kind"] == "KEYWORD_ONLY"}
    last_baseline_index = max(
        (i for i, n in enumerate(actual_names) if n in expected_names),
        default=-1,
    )
    for index, name in enumerate(actual_names):
        if name in expected_names or name in kw_only:
            continue
        assert index > last_baseline_index, (
            f"{modname}.{fnname}: new parameter '{name}' was inserted at "
            f"position {index}, among the v2.0.0 parameters. Append it after "
            f"the existing parameters, or make it keyword-only."
        )
