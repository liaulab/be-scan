import matplotlib.pyplot as plt
import pytest

@pytest.fixture(autouse=True)
def matplotlib_hide_plots(monkeypatch):
    # prevent plots from opening during test
    monkeypatch.setattr(plt, "show", lambda: None)
