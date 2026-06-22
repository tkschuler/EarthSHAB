"""Shared forecast-core utilities (vendored from HAB-COM).

wind_interp / forecast_backends / vectorization / transform are kept byte-
identical to HAB-COM except for package-prefixed cross-imports, so stencil/
interpolation edits stay portable between the two repos until the shared core
is extracted into a standalone package.
"""
