"""Logging for xtraee.

Importing the package never changes the level of the root logger: everything
xtraee logs goes through the ``xtraee`` logger tree, which defaults to
warnings and errors.  The CLI calls :func:`configure` once it knows whether
``--debug`` was requested, so ``--debug`` actually controls verbosity instead
of every parser shouting to stdout.
"""

import logging
import sys

LOGGER_NAME = "xtraee"

#: True when debug logging (and the internal consistency assertions) is enabled.
debug = False
_LEVEL = logging.WARNING

_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
_DATEFMT = "%Y-%m-%d %H:%M:%S"


class ColorFormatter(logging.Formatter):
    COLORS = {
        "DEBUG": "\033[94m",  # Blue
        "INFO": "\033[92m",  # Green
        "WARNING": "\033[93m",  # Yellow
        "ERROR": "\033[91m",  # Red
        "CRITICAL": "\033[95m",  # Magenta
    }
    RESET = "\033[0m"

    def format(self, record):
        # wrap the *message* in color codes, not the entire record
        msg = super().format(record)
        color = self.COLORS.get(record.levelname, "")
        return f"{color}{msg}{self.RESET}"


def configure(verbose: bool = False) -> None:
    """Set the verbosity of the whole ``xtraee`` logger tree.

    ``verbose=False`` keeps only warnings and errors, so a normal run stays
    quiet.  ``verbose=True`` (the CLI's ``--debug``) turns on debug messages
    and the internal consistency checks that read :data:`debug`.
    """
    global _LEVEL, debug
    debug = bool(verbose)
    _LEVEL = logging.DEBUG if debug else logging.WARNING

    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(_LEVEL)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(ColorFormatter(_FORMAT, _DATEFMT))
        logger.addHandler(handler)
    logger.propagate = False


def get_logger(name: str | None = None) -> logging.Logger:
    """Return the logger for one module, under the shared ``xtraee`` tree."""
    if not name:
        return logging.getLogger(LOGGER_NAME)
    return logging.getLogger(f"{LOGGER_NAME}.{name}")


configure()
