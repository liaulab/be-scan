import subprocess
import re

# make sure every subcommand has help text
def test_subcommand_help():
    out = subprocess.run("python -m be_scan -h", shell=True, capture_output=True)
    match = re.search(r"positional arguments:\n\s+\{([A-Za-z_,]+)\}", out.stdout.decode(), re.MULTILINE)
    assert match and match[1]
    subcommands = match[1].split(",")
    for subcommand in subcommands:
        assert re.search(rf"{subcommand}\s+.+", out.stdout.decode()), f"Missing help text for subcommand `{subcommand}`"
