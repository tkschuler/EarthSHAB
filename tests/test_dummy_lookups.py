"""Lookup-helper tests against synthetic forecasts.

Per Q9(c) "lookup tier": getNearestLat / getNearestLon / getNearestAlt are
exercised via the reader_adapter against the scenarios that actually stress
the lookup logic for each axis.

Scenarios used here intentionally cover only the lookups (not the
interpolation composite, which has its own file). For each, the assertion is:

    "the reader returns the index of the grid point closest to the query,
    and a query at exactly a grid value returns that grid value's index."

The cross-reader equality assertion verifies that GFS and ERA5 agree on
nearest-index in the GFS-ascending-lat convention (the ERA5 adapter handles
the descending convention internally — at the end of the day both readers
should point at the same physical lat/lon).
"""

from __future__ import annotations

import numpy as np
import pytest

from EarthSHAB.GFS import GFS

from tests.conftest import ON_GRID_EPS
from tests.tools.make_dummy_forecasts import (
    LAT_VALUES_ASC,
    LON_VALUES,
    ALT_M_GFS_ORDER,
)
from tests.tools.reader_adapter import (
    is_gfs,
    nearest_alt,
    nearest_lat,
    nearest_lon,
)


# ---------------------------------------------------------------------------
# Helpers to compare an index returned by either reader to a physical lat/lon.
# ---------------------------------------------------------------------------

def _stored_lat(reader, idx):
    """Return the physical lat (degrees) at the reader's stored lat[idx]."""
    return float(reader.lat[idx])


def _stored_lon(reader, idx):
    return float(reader.lon[idx])


# ---------------------------------------------------------------------------
# getNearestLat
# ---------------------------------------------------------------------------

# Use only the scenarios that exercise lat-axis logic.
@pytest.fixture(params=["lat_ramp_only", "linear_ramp_all_dims", "all_static"])
def scenario(request):
    return request.param


class TestNearestLat:
    def test_query_at_exact_grid_lat_returns_that_index(self, reader):
        for target in LAT_VALUES_ASC:
            idx = nearest_lat(reader, float(target))
            assert _stored_lat(reader, idx) == pytest.approx(
                float(target), abs=ON_GRID_EPS
            )

    def test_query_below_min_lat_returns_min_endpoint(self, reader):
        idx = nearest_lat(reader, float(LAT_VALUES_ASC.min()) - 10.0)
        # GFS lat ascending: nearest to a far-below query is index 0 (= min).
        # ERA5 lat descending: nearest is the LAST stored index (= min).
        # In either case the PHYSICAL lat is the minimum stored value.
        assert _stored_lat(reader, idx) == pytest.approx(
            float(LAT_VALUES_ASC.min()), abs=ON_GRID_EPS
        )

    def test_query_above_max_lat_returns_max_endpoint(self, reader):
        idx = nearest_lat(reader, float(LAT_VALUES_ASC.max()) + 10.0)
        assert _stored_lat(reader, idx) == pytest.approx(
            float(LAT_VALUES_ASC.max()), abs=ON_GRID_EPS
        )

    def test_query_off_grid_picks_nearest(self, reader):
        # Halfway between two grid points → either is acceptable.
        # Pick a query that's 60%/40% between grid[2] and grid[3]:
        target = 0.6 * LAT_VALUES_ASC[2] + 0.4 * LAT_VALUES_ASC[3]
        idx = nearest_lat(reader, float(target))
        # Should resolve to grid[2] (the closer one).
        assert _stored_lat(reader, idx) == pytest.approx(
            float(LAT_VALUES_ASC[2]), abs=ON_GRID_EPS
        )


# ---------------------------------------------------------------------------
# getNearestLon
# ---------------------------------------------------------------------------

class TestNearestLon:
    def test_query_at_exact_grid_lon_returns_that_index(self, reader):
        for target in LON_VALUES:
            idx = nearest_lon(reader, float(target))
            assert _stored_lon(reader, idx) == pytest.approx(
                float(target), abs=ON_GRID_EPS
            )

    def test_query_off_grid_picks_nearest(self, reader):
        target = 0.7 * LON_VALUES[4] + 0.3 * LON_VALUES[5]
        idx = nearest_lon(reader, float(target))
        assert _stored_lon(reader, idx) == pytest.approx(
            float(LON_VALUES[4]), abs=ON_GRID_EPS
        )


class TestNearestLonGfsWrap:
    """GFS-specific: lon is in [0, 360) and queries should be normalized
    via % 360. This isn't an ERA5 concern, so it doesn't parametrize."""

    @pytest.fixture(params=["lon_ramp_only", "all_static"])
    def scenario(self, request):
        return request.param

    @pytest.fixture(params=["gfs"])
    def reader_type(self, request):
        return request.param

    def test_negative_360_equivalent_query(self, reader):
        # Synthetic grid lon=1.0 exists. Query at lon = -359.0 (% 360 = 1.0)
        # should resolve to the same index as query at lon = 1.0.
        assert is_gfs(reader)
        idx_pos = nearest_lon(reader, 1.0)
        idx_neg = nearest_lon(reader, -359.0)
        assert idx_pos == idx_neg


# ---------------------------------------------------------------------------
# getNearestAlt
# ---------------------------------------------------------------------------

class TestNearestAlt:
    @pytest.fixture(params=["altitude_ramp_only", "geopotential_vs_height"])
    def scenario(self, request):
        return request.param

    def test_query_at_exact_grid_altitude(self, reader, synthetic_start_dt):
        """For each stored altitude, querying that altitude resolves to its
        own index. Note hour_idx 0 — altitude column is constant in time
        for the synthetic scenarios, so this is fine.
        """
        # ALT_M_GFS_ORDER is the physical altitude column. For ERA5 the
        # stored z is altitude*g, but getNearestAlt resolves against
        # whatever altitude the reader uses internally — which is z/g for
        # ERA5 after the constructor's /g conversion. So the test
        # interrogates the reader's *internal* altitude column.
        # The lat/lon arguments here pick a grid point; geopotential_vs_height
        # has constant altitude in lat/lon so the choice doesn't matter.
        lat = float(reader.lat[len(reader.lat) // 2])
        lon = float(reader.lon[len(reader.lon) // 2])
        for target_alt in ALT_M_GFS_ORDER:
            idx = nearest_alt(reader, 0, lat, lon, float(target_alt))
            # Verify the reader's stored hgtprs at that index is within
            # ON_GRID_EPS (small float32 round-trip noise) of the queried
            # altitude. The reader's hgtprs is shape
            # (time, level, lat, lon). ERA5 has already divided by g in
            # __init__, so hgtprs is in meters either way.
            lat_i = nearest_lat_helper_index(reader, lat)
            lon_i = nearest_lon_helper_index(reader, lon)
            stored_alt = float(reader.hgtprs[0, idx, lat_i, lon_i])
            assert stored_alt == pytest.approx(
                float(target_alt), abs=ON_GRID_EPS * max(abs(target_alt), 1.0)
            )

    def test_below_lowest_altitude_returns_bottom_index(
        self, reader, synthetic_start_dt
    ):
        lat = float(reader.lat[0])
        lon = float(reader.lon[0])
        # Synthetic min altitude = 100 m. Query at -1000 should clamp.
        idx = nearest_alt(reader, 0, lat, lon, -1000.0)
        lat_i = nearest_lat_helper_index(reader, lat)
        lon_i = nearest_lon_helper_index(reader, lon)
        stored_alt = float(reader.hgtprs[0, idx, lat_i, lon_i])
        # Should resolve to the bottom-most stored altitude (closest to
        # below-min query).
        assert stored_alt == pytest.approx(
            float(ALT_M_GFS_ORDER.min()), abs=ON_GRID_EPS * 100.0
        )


def nearest_lat_helper_index(reader, lat):
    return nearest_lat(reader, lat)


def nearest_lon_helper_index(reader, lon):
    return nearest_lon(reader, lon)


# ---------------------------------------------------------------------------
# Cross-reader index agreement on physical lat/lon
# ---------------------------------------------------------------------------

class TestCrossReaderLookupAgreement:
    """Both readers, fed the same physical query, must point at the same
    physical grid cell — even though their underlying indices may differ
    because of the GFS-asc vs ERA5-desc lat convention.

    The assertion is on the PHYSICAL value of the resolved cell, not on
    the integer index (which is schema-dependent).
    """

    @pytest.fixture(params=["all_static"])  # one scenario is enough
    def scenario(self, request):
        return request.param

    def test_same_lat_resolves_to_same_physical_lat(
        self,
        scenario,
        all_dummies,
        monkeypatch,
        synthetic_start_dt,
    ):
        # We need both a GFS and ERA5 reader for the same scenario. The
        # reader fixture only gives us one — so build both manually here.
        from EarthSHAB.ERA5 import ERA5
        from EarthSHAB.GFS import GFS
        import EarthSHAB.config_earth as config_earth
        from tests.conftest import _patch_config_for_gfs, _patch_config_for_era5

        path_gfs = all_dummies[(scenario, "gfs")]
        path_era5 = all_dummies[(scenario, "era5")]

        _patch_config_for_gfs(monkeypatch, path_gfs, synthetic_start_dt)
        gfs_reader = GFS(config_earth.simulation["start_coord"])

        # Re-patch for ERA5 (monkeypatch will undo all of it after the test).
        _patch_config_for_era5(monkeypatch, path_era5, synthetic_start_dt)
        era5_reader = ERA5(config_earth.simulation["start_coord"])

        for target_lat in LAT_VALUES_ASC:
            gfs_idx = nearest_lat(gfs_reader, float(target_lat))
            era5_idx = nearest_lat(era5_reader, float(target_lat))
            assert _stored_lat(gfs_reader, gfs_idx) == pytest.approx(
                _stored_lat(era5_reader, era5_idx), abs=ON_GRID_EPS
            )
