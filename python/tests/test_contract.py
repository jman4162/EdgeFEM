"""Tests for the producer-owned artifact contract fixture.

These pin the format EdgeFEM guarantees to consumers. If a change here breaks
one of these, the consumer loaders (opensatcom's edgefem_json_loader) break
with it — bump edgefem.contract.FIXTURE_REVISION and refresh the vendored
fixtures downstream instead of editing the assertions.
"""

import json
import math

import pytest

pyedgefem = pytest.importorskip("pyedgefem", reason="pyedgefem binding not built")

from edgefem import contract  # noqa: E402


@pytest.fixture(scope="module")
def written(tmp_path_factory):
    out = tmp_path_factory.mktemp("contract")
    json_path, csv_path = contract.write_golden_array_package(out)
    return json_path, csv_path


class TestArrayPackageJson:
    def test_handshake_fields(self, written):
        doc = json.loads(written[0].read_text())
        assert doc["format"] == "EdgeFEM_ArrayPackage"
        assert doc["version"] == "1.0"
        assert doc["model_name"] == f"golden-2x2-rev{contract.FIXTURE_REVISION}"

    def test_shapes_and_units(self, written):
        doc = json.loads(written[0].read_text())
        assert doc["num_elements"] == 4
        assert all(len(p) == 3 for p in doc["element_positions"])  # metres, 3-vectors
        assert doc["frequencies"] == sorted(doc["frequencies"])  # Hz, ascending
        assert doc["num_frequencies"] == len(doc["frequencies"]) == 2
        assert len(doc["S_matrices"]) == 2
        for entry in doc["S_matrices"]:
            assert entry["frequency"] in doc["frequencies"]
            data = entry["data"]
            assert len(data) == 4 and all(len(row) == 4 for row in data)
            # complex entries serialize as [re, im] pairs
            assert all(len(cell) == 2 for row in data for cell in row)

    def test_pinned_s_values(self, written):
        """The diagonal is -15 dB at +30 deg; adjacent coupling -20 dB at -45 deg."""
        doc = json.loads(written[0].read_text())
        data = doc["S_matrices"][0]["data"]
        s00 = complex(data[0][0][0], data[0][0][1])
        assert abs(s00) == pytest.approx(10 ** (-15 / 20), rel=1e-12)
        assert math.degrees(math.atan2(s00.imag, s00.real)) == pytest.approx(30.0, abs=1e-9)
        s01 = complex(data[0][1][0], data[0][1][1])
        assert abs(s01) == pytest.approx(10 ** (-20 / 20), rel=1e-12)

    def test_embedded_patterns_are_summary_only(self, written):
        """The JSON deliberately carries only a pattern summary; the grid data
        travels in the companion CSV. Consumers must not expect more."""
        doc = json.loads(written[0].read_text())
        assert doc["has_embedded_patterns"] is False
        assert "embedded_patterns" not in doc


class TestPatternsCsv:
    def test_columns_and_units(self, written):
        lines = written[1].read_text().splitlines()
        assert lines[0] == (
            "port_index,theta_deg,phi0_mag,phi0_phase_deg,phi90_mag,phi90_phase_deg"
        )
        # 4 ports x 37 theta samples
        assert len(lines) - 1 == 4 * 37
        first = lines[1].split(",")
        assert first[0] == "0"
        assert float(first[1]) == -90.0  # degrees in the CSV, radians internally

    def test_boresight_magnitude(self, written):
        """cos(0) = 1 for the E-plane cut, 0.9 for the H-plane cut."""
        rows = [line.split(",") for line in written[1].read_text().splitlines()[1:]]
        boresight = [r for r in rows if r[0] == "0" and float(r[1]) == 0.0]
        assert len(boresight) == 1
        assert float(boresight[0][2]) == pytest.approx(1.0, abs=1e-12)
        assert float(boresight[0][4]) == pytest.approx(0.9, abs=1e-12)


def test_deterministic(tmp_path):
    a = contract.write_golden_array_package(tmp_path / "a")
    b = contract.write_golden_array_package(tmp_path / "b")
    assert a[0].read_bytes() == b[0].read_bytes()
    assert a[1].read_bytes() == b[1].read_bytes()
