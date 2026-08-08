"""
Shared pytest configuration for the be_scan test suite.

Two things happen here before any test module is imported:

1. The repository root is pushed to the FRONT of sys.path. A copy of be_scan
   is often present in site-packages, and without this the suite may silently
   validate the installed copy instead of the working tree.
2. Matplotlib is switched to the non-interactive "Agg" backend. Several
   be_scan modules import pyplot at import time, so this has to happen before
   the first `import be_scan`.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# --- 1. make sure we test the working tree, not an installed copy ---------- #
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

# --- 2. headless plotting -------------------------------------------------- #
import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

TEST_DATA = REPO_ROOT / "tests" / "test_data"


def pytest_report_header(config):
    """Print which be_scan is under test, so a shadowed import is obvious."""
    import be_scan

    return f"be_scan under test: {be_scan.__file__}"


@pytest.fixture(scope="session")
def repo_root():
    """Absolute path to the repository root."""
    return REPO_ROOT


@pytest.fixture(scope="session")
def test_data():
    """Absolute path to tests/test_data."""
    return TEST_DATA


@pytest.fixture
def tmp_out(tmp_path):
    """
    A scratch output directory for functions that write files.

    Every be_scan function that saves output takes some form of out_dir /
    output_dir. Point it here so no test ever writes into the repository.
    """
    out = tmp_path / "out"
    out.mkdir(parents=True, exist_ok=True)
    return out


@pytest.fixture(autouse=True)
def _close_figures():
    """Close any figures a test leaves open, so plots never leak between tests."""
    yield
    plt.close("all")
