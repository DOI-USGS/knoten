# Changelog

All changes that impact users of this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!---
This document is intended for users of the applications and API. Changes to things
like tests should not be noted in this document.

When updating this file for a PR, add an entry for your change under Unreleased
and one of the following headings:
 - Added - for new features.
 - Changed - for changes in existing functionality.
 - Deprecated - for soon-to-be removed features.
 - Removed - for now removed features.
 - Fixed - for any bug fixes.
 - Security - in case of vulnerabilities.

If the heading does not yet exist under Unreleased, then add it as a 3rd heading,
with three #.


When preparing for a public release candidate add a new 2nd heading, with two #, under
Unreleased with the version number and the release date, in year-month-day
format. Then, add a link for the new version at the bottom of this document and
update the Unreleased link so that it compares against the latest release tag.


When preparing for a bug fix release create a new 2nd heading above the Fixed
heading to indicate that only the bug fixes and security fixes are in the bug fix
release.
-->

## Unreleased

### Added
- `create_csm` now dispatches to `_from_isd` and `_from_state` to test whether the sensor model can be instantiated from either and ISD or a state file.
- `generate_image_coordinate` to `csm.py`. This provides a similar interface to `generate_ground_coordinate` and abstracts away the `csmapi` from the user.
- A surface class (moved from AutoCNet; credit @jessemapel) with support for Ellipsoid DEMs and basic support for raster DEMs readable by the plio.io.io_gdal.GeoDataset. Support is basic because it uses a single pixel intersection and not an interpolated elevation like ISIS does.
- A check to `generate_ground_point` when a GeoDataset is used to raise a `ValueError` if the algorithm intersects a no data value in the passed DEM. This ensures that valid heights are used in the intersection computation. Fixes [#120](https://github.com/DOI-USGS/knoten/issues/120)

### Changed
- Removed all `pyproj` calls from csm.py, abstracting them into the reprojection and pyproj.Transformer code inside utils.py. Updated the transformations to use the new pipeline style syntax to avoid deprecation warnings about old syntax.

### Fixed
- The init method that searches for the libusgscsm to support searching in the `csmplugins` subdirectory. This approach depends on being able to find `csmapi` in a standard location and then assumes that the `libusgscsm` shared library is in a subdirectoy of that `lib` directory. Fixes [#118](https://github.com/DOI-USGS/knoten/issues/118)

