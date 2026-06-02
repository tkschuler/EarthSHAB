"""Lookup-helper tests against synthetic forecasts.

Per Q9(c) "lookup tier": getNearestLatIdx / getNearestLonIdx /
getNearestAlt are exercised against the scenarios that actually stress the
lookup logic for each axis.

Scenarios used here intentionally cover only the lookups (not the
interpolation composite, which has its own file). For each, the assertion is:

    "the reader returns the index of the grid point closest to the query,
    and a query at exactly a grid value returns that grid value's index."

v2 collapsed the two readers into a single Forecast class, so the
cross-reader equality assertion that used to live here is now redundant
and has been removed. The GFS-specific lon-wrap test is also gone — the
v2 canonical schema stores longitude in [-180, 180), so a query at lon
= -359 is no longer expected to wrap to lon = 1.
"""

from __future__ import annotations

import pytest

from tests.conftest import ON_GRID_EPS
from tests.tools.make_dummy_forecasts import (
    LAT_VALUES_ASC,
    LON_VALUES,
    ALT_M_GFS_ORDER,
)


# ---------------------------------------------------------------------------
# Helpers to compare an index returned by the reader to a physical lat/lon.
# ---------------------------------------------------------------------------

def _stored_lat(reader, idx):
    """Return the physical lat (degrees) at the reader's stored lat[idx]."""
    return float(reader.lat[idx])


def _stored_lon(reader, idx):
    return float(reader.lon[idx])


def _nearest_lat(reader, lat):
    return reader.getNearestLatIdx(float(lat))


def _nearest_lon(reader, lon):
    return reader.getNearestLonIdx(float(lon))


def _nearest_alt(reader, hour_idx, lat, lon, alt_m):
    lat_i = _nearest_lat(reader, lat)
    lon_i = _nearest_lon(reader, lon)
    return reader.getNearestAltbyIndex(int(hour_idx), lat_i, lon_i, alt_m)


# ---------------------------------------------------------------------------
# getNearestLatIdx
# ---------------------------------------------------------------------------

# Use only the scenarios that exercise lat-axis logic.
@pytest.fixture(params=["lat_ramp_only", "linear_ramp_all_dims", "all_static"])
def scenario(request):
    return request.param


class TestNearestLat:
    def test_query_at_exact_grid_lat_returns_that_index(self, reader):
        for target in LAT_VALUES_ASC:
            idx = _nearest_lat(reader, float(target))
            assert _stored_lat(reader, idx) == pytest.approx(
                float(target), abs=ON_GRID_EPS
            )

    def test_query_below_min_lat_returns_min_endpoint(self, reader):
        idx = _nearest_lat(reader, float(LAT_VALUES_ASC.min()) - 10.0)
        # Canonical v2 lat is descending: the LAST stored index is the
        # minimum physical lat. closest() of a far-below query lands there.
        assert _stored_lat(reader, idx) == pytest.approx(
            float(LAT_VALUES_ASC.min()), abs=ON_GRID_EPS
        )

    def test_query_above_max_lat_returns_max_endpoint(self, reader):
        idx = _nearest_lat(reader, float(LAT_VALUES_ASC.max()) + 10.0)
        assert _stored_lat(reader, idx) == pytest.approx(
            float(LAT_VALUES_ASC.max()), abs=ON_GRID_EPS
        )

    def test_query_off_grid_picks_nearest(self, reader):
        # Off-grid: pick a query that's 60%/40% between grid[2] and grid[3]:
        target = 0.6 * LAT_VALUES_ASC[2] + 0.4 * LAT_VALUES_ASC[3]
        idx = _nearest_lat(reader, float(target))
        # Should resolve to grid[2] (the closer one).
        assert _stored_lat(reader, idx) == pytest.approx(
            float(LAT_VALUES_ASC[2]), abs=ON_GRID_EPS
        )


# ---------------------------------------------------------------------------
# getNearestLonIdx
# ---------------------------------------------------------------------------

class TestNearestLon:
    def test_query_at_exact_grid_lon_returns_that_index(self, reader):
        for target in LON_VALUES:
            idx = _nearest_lon(reader, float(target))
            assert _stored_lon(reader, idx) == pytest.approx(
                float(target), abs=ON_GRID_EPS
            )

    def test_query_off_grid_picks_nearest(self, reader):
        target = 0.7 * LON_VALUES[4] + 0.3 * LON_VALUES[5]
        idx = _nearest_lon(reader, float(target))
        assert _stored_lon(reader, idx) == pytest.approx(
            float(LON_VALUES[4]), abs=ON_GRID_EPS
        )


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
        # ALT_M_GFS_ORDER is the physical altitude column. The reader stores
        # z (geopotential, m²/s²) and divides by g in __init__, so hgtprs
        # is altitude in meters regardless of source.
        lat = float(reader.lat[len(reader.lat) // 2])
        lon = float(reader.lon[len(reader.lon) // 2])
        for target_alt in ALT_M_GFS_ORDER:
            idx = _nearest_alt(reader, 0, lat, lon, float(target_alt))
            lat_i = _nearest_lat(reader, lat)
            lon_i = _nearest_lon(reader, lon)
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
        idx = _nearest_alt(reader, 0, lat, lon, -1000.0)
        lat_i = _nearest_lat(reader, lat)
        lon_i = _nearest_lon(reader, lon)
        stored_alt = float(reader.hgtprs[0, idx, lat_i, lon_i])
        # Should resolve to the bottom-most stored altitude (closest to
        # below-min query).
        assert stored_alt == pytest.approx(
            float(ALT_M_GFS_ORDER.min()), abs=ON_GRID_EPS * 100.0
        )
