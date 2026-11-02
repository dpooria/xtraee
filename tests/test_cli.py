"""End-to-end CLI behaviour, driven through ``main(argv)``."""

from pathlib import Path

import fixtures
import pytest

from xtraee import __version__
from xtraee.cli import inject_default_command, main

REAL_LOG = Path(__file__).parent / "test_qccsd" / "input2.log"


def test_extract_writes_the_overview_file(tmp_path):
    exit_code = main([str(REAL_LOG), "--outdir", str(tmp_path)])
    assert exit_code == 0
    produced = tmp_path / "data_EOM-CCSD.txt"
    assert produced.exists()
    assert "---- Happy family ----" in produced.read_text()


def test_extract_full_writes_both_files(tmp_path):
    assert main([str(REAL_LOG), "--full", "--outdir", str(tmp_path)]) == 0
    assert (tmp_path / "fulldata_EOM-CCSD.txt").exists()
    assert (tmp_path / "data_EOM-CCSD.txt").exists()


def test_compare_prints_the_accuracy(tmp_path, capsys):
    """The whole point of ``compare``: its result must reach the terminal.

    It used to be logged at INFO while the logger sat at WARNING, so the
    command printed nothing at all and exited 0.
    """
    exit_code = main(
        [
            "compare",
            str(REAL_LOG),
            str(REAL_LOG),
            "CCSD:singlet-1/A",
            "CCSD:singlet-1/A",
            "--acc-method",
            "cosine",
            "--outdir",
            str(tmp_path),
        ]
    )
    out = capsys.readouterr().out
    assert exit_code == 0
    assert "Comparing singlet-1/A and singlet-1/A" in out
    assert "accuracy" in out
    assert "1.0000000000" in out


def test_compare_reports_a_missing_state(tmp_path, logs):
    exit_code = main(
        [
            "compare",
            str(REAL_LOG),
            str(REAL_LOG),
            "CCSD:singlet-99/A",
            "CCSD:singlet-1/A",
            "--outdir",
            str(tmp_path),
        ]
    )
    assert exit_code == 1
    assert "not found" in logs.text


def test_compare_reports_a_malformed_state_name(tmp_path, logs):
    exit_code = main(
        ["compare", str(REAL_LOG), str(REAL_LOG), "noseparator", "CIS:singlet-1/A",
         "--outdir", str(tmp_path)]
    )
    assert exit_code == 1
    assert "format is" in logs.text


def test_compareall_writes_the_csv_matrix(tmp_path):
    assert main(["compareall", str(REAL_LOG), "--outdir", str(tmp_path)]) == 0
    csv = tmp_path / "A_A_EOM-CCSD_compare.csv"
    assert csv.exists()
    header = csv.read_text().splitlines()[0]
    assert header.startswith(",triplet-1/A")
    assert "singlet-1/A" in csv.read_text()


def test_descriptors_writes_the_csv_into_outdir(tmp_path, capsys):
    """``--descriptors-file`` used to be resolved against the CWD."""
    exit_code = main(["descriptors", str(REAL_LOG), "--outdir", str(tmp_path)])
    out = capsys.readouterr().out
    produced = tmp_path / "descriptors.csv"
    assert exit_code == 0
    assert produced.exists()
    assert f"Wrote 15 states with descriptors to {produced}" in out
    header = produced.read_text().splitlines()[0].split(",")
    assert "identifier" in header
    assert "frontier_no_1" in header


def test_descriptors_accepts_a_path_and_a_custom_name(tmp_path):
    target = tmp_path / "nested" / "meta.csv"
    assert main(
        ["descriptors", str(REAL_LOG), "--outdir", str(tmp_path / "out"),
         "--descriptors-file", str(target)]
    ) == 0
    assert target.exists()

    renamed = main(
        ["descriptors", str(REAL_LOG), "--outdir", str(tmp_path),
         "--descriptors-file", "renamed.csv"]
    )
    assert renamed == 0
    assert (tmp_path / "renamed.csv").exists()


def test_descriptors_refuses_a_non_ccsd_input(tmp_path, logs):
    cis = fixtures.write(tmp_path, "cis.log", fixtures.cis_log())
    exit_code = main(["descriptors", str(cis), "--outdir", str(tmp_path / "out")])
    assert exit_code == 1
    assert "only implemented for EOM-CCSD" in logs.text


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_default_command_injection():
    assert inject_default_command(["a.log"]) == ["extract", "a.log"]
    assert inject_default_command([]) == ["extract"]
    # an option value must not be mistaken for the command
    assert inject_default_command(["--threshold", "0.2", "a.log"]) == [
        "extract", "--threshold", "0.2", "a.log"
    ]
    assert inject_default_command(["--out-data", "compare", "a.log"]) == [
        "extract", "--out-data", "compare", "a.log"
    ]
    # a command further right is hoisted in front of the shared options
    assert inject_default_command(["--outdir", "out", "compare", "a.log"]) == [
        "compare", "--outdir", "out", "a.log"
    ]
    assert inject_default_command(["a.log", "--full"]) == ["extract", "a.log", "--full"]
    # help and version keep their own meaning
    assert inject_default_command(["--help"]) == ["--help"]


def test_options_before_the_command_still_work(tmp_path):
    exit_code = main(["--threshold", "0.2", "--outdir", str(tmp_path), str(REAL_LOG)])
    assert exit_code == 0
    assert (tmp_path / "data_EOM-CCSD.txt").exists()


def test_unknown_accuracy_method_fails_loudly(tmp_path):
    with pytest.raises(ValueError, match="Method not recognized"):
        main(
            ["compare", str(REAL_LOG), str(REAL_LOG),
             "CCSD:singlet-1/A", "CCSD:singlet-1/A",
             "--acc-method", "nope", "--outdir", str(tmp_path)]
        )
