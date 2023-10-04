import pyproj

def reproject(record, semi_major, semi_minor, source_proj, dest_proj, **kwargs):
    """
    Thin wrapper around PyProj's Transform() function to transform 1 or more three-dimensional
    point from one coordinate system to another. If converting between Cartesian
    body-centered body-fixed (BCBF) coordinates and Longitude/Latitude/Altitude coordinates,
    the values input for semi-major and semi-minor axes determine whether latitudes are
    planetographic or planetocentric and determine the shape of the datum for altitudes.
    If semi_major == semi_minor, then latitudes are interpreted/created as planetocentric
    and altitudes are interpreted/created as referenced to a spherical datum.
    If semi_major != semi_minor, then latitudes are interpreted/created as planetographic
    and altitudes are interpreted/created as referenced to an ellipsoidal datum.

    Parameters
    ----------
    record : object
          Pandas series object

    semi_major : float
              Radius from the center of the body to the equater

    semi_minor : float
              Radius from the pole to the center of mass

    source_proj : str
                      Pyproj string that defines a projection space ie. 'geocent'

    dest_proj : str
                   Pyproj string that defines a project space ie. 'latlon'

    Returns
    -------
    : list
    Transformed coordinates as y, x, z

    """
    transformer = pyproj.Transformer.from_crs(f'+proj={source_proj} +a={semi_major} +b={semi_minor}',
                                              f'+proj={dest_proj} +a={semi_major} +b={semi_minor}',
                                              always_xy=True)
    source_proj = f'+proj={source_proj} +a={semi_major} +b={semi_minor}'
    dest_proj = f'+proj={dest_proj} +a={semi_major} +b={semi_minor}'
    transformer = create_transformer(source_proj, dest_proj)
    x, y, z = transformer.transform(record[0], record[1], record[2], errcheck=True)

    return y, x, z

def create_transformer(source_proj, dest_proj):
    return pyproj.Transformer.from_crs(source_proj,
                                       dest_proj,
                                       always_xy=True)