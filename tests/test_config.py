"""Logging: quiet by default, ``--debug`` turns the detail on."""

import logging

from xtraee import config


def teardown_function():
    config.configure(False)


def test_package_logger_is_quiet_by_default():
    assert config.debug is False
    logger = logging.getLogger("xtraee")
    assert logger.level == logging.WARNING
    # children inherit, so a parser debug message is dropped
    assert config.get_logger("QCCSDParser").isEnabledFor(logging.DEBUG) is False


def test_configure_debug_enables_debug_messages():
    config.configure(True)
    assert config.debug is True
    assert config.get_logger("QCCSDParser").isEnabledFor(logging.DEBUG) is True


def test_configure_is_reversible():
    config.configure(True)
    config.configure(False)
    assert config.debug is False
    assert logging.getLogger("xtraee").level == logging.WARNING


def test_get_logger_returns_children_of_one_tree(logs):
    logger = config.get_logger("unit")
    assert logger.name == "xtraee.unit"
    logger.warning("hello")
    assert "hello" in logs.text
    assert "WARNING" in logs.text


def test_importing_the_package_does_not_touch_the_root_logger():
    assert logging.getLogger().level == logging.WARNING
    assert logging.getLogger("xtraee").propagate is False
