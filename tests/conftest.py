import logging
import sys
from pathlib import Path

import pytest

# ``pytest`` from a plain checkout must be able to import the package and the
# fixture module that live next to this file.
ROOT = Path(__file__).resolve().parent.parent
for path in (str(ROOT), str(ROOT / "tests")):
    if path not in sys.path:
        sys.path.insert(0, path)


@pytest.fixture
def logs(caplog, monkeypatch):
    """Capture ``xtraee`` records.

    The package tree sets ``propagate = False`` and writes to its own handler,
    so pytest's ``caplog`` cannot see it by default; let the records through
    to the root logger for the duration of the test.
    """
    from xtraee import config

    logger = logging.getLogger("xtraee")
    previous = logger.propagate

    # the CLI reconfigures the tree, which turns propagation back off
    real_configure = config.configure

    def configure_and_propagate(verbose: bool = False) -> None:
        real_configure(verbose)
        logger.propagate = True

    monkeypatch.setattr("xtraee.config.configure", configure_and_propagate)
    monkeypatch.setattr("xtraee.cli.configure", configure_and_propagate)

    logger.propagate = True
    caplog.set_level(logging.DEBUG, logger="xtraee")
    try:
        yield caplog
    finally:
        logger.propagate = previous
        logger.setLevel(logging.WARNING)
