"""
Regenerate the public-API baseline snapshot.

    python tests/api_test/generate_baseline.py

DO NOT run this to "fix" a failing test_api_compat.py.

be_scan is published and in use by other labs. The baseline exists so that a
change which would break somebody's existing call fails loudly in CI. If
test_api_compat.py fails, the correct response is almost always to change the
code back, not to regenerate the baseline.

Regenerate only when a breaking change has been deliberately decided on and
the major version is being bumped. When you do, commit the regenerated JSON in
its own commit with an explanation of what broke and why.
"""

import argparse
import importlib
import inspect
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_OUT = Path(__file__).resolve().parent / "api_baseline_v2.json"

MODULES = [
    "be_scan.sgrna",
    "be_scan.analysis",
    "be_scan.plot",
    "be_scan.figure_plot",
]


def encode_default(value):
    """
    Render a parameter default as something JSON can round-trip.

    Plain JSON types are stored as-is so that a change from e.g. [4, 8] to
    [4, 9] is visible in a diff. Anything else (dataclass instances, compiled
    regexes, functions) falls back to its repr, prefixed so the comparison
    logic knows not to treat it as a literal.

    The value is round-tripped through JSON so that live signatures and the
    stored baseline are normalized the same way. Without this, a tuple default
    such as (600, 300) would read back from the baseline as a list and compare
    unequal to itself.
    """
    if value is inspect.Parameter.empty:
        return "<required>"
    try:
        return json.loads(json.dumps(value))
    except (TypeError, ValueError):
        return f"<repr>{value!r}"


def snapshot():
    """Walk the public modules and record every exported be_scan function."""
    api = {}
    for modname in MODULES:
        module = importlib.import_module(modname)
        functions = {}
        for name in sorted(dir(module)):
            if name.startswith("_"):
                continue
            obj = getattr(module, name)
            if not inspect.isfunction(obj):
                continue
            # figure_plot does `from figure_classes import *`, which pulls in
            # third-party functions. Only record what be_scan actually defines.
            if not getattr(obj, "__module__", "").startswith("be_scan"):
                continue
            functions[name] = {
                "module": obj.__module__,
                "params": [
                    {
                        "name": p.name,
                        "kind": str(p.kind),
                        "default": encode_default(p.default),
                    }
                    for p in inspect.signature(obj).parameters.values()
                ],
            }
        api[modname] = functions
    return api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    api = snapshot()
    args.out.write_text(json.dumps(api, indent=2, sort_keys=True) + "\n")

    n_fns = sum(len(v) for v in api.values())
    n_params = sum(len(f["params"]) for v in api.values() for f in v.values())
    print(f"Wrote {args.out}")
    print(f"  {n_fns} public functions, {n_params} parameters")


if __name__ == "__main__":
    import sys

    sys.path.insert(0, str(REPO_ROOT))
    main()
