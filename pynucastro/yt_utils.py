import re

from pynucastro.networks import Composition


def get_point(ds, pos):
    """Get a specific point from a yt dataset object.

    Field values can be extracted by subscripting: ``point['Temp']``.

    Works around various annoyances with the argument handling in
    ``ds.point()``, so things like ``get_point(ds, ds.all_data().argmax('enuc'))``
    work as expected for all coordinate systems.

    Parameters
    ----------

    ds : yt.data_objects.static_output.Dataset
        a yt dataset
    pos : Iterable(float, ~unyt.array.unyt_quantity)
        The position to get; must have at least ``ds.dimensionality`` elements.

    Returns
    -------
    yt.data_objects.selection_objects.point.YTPoint

    """
    # lazy import so we don't add a hard dependency on unyt
    # pylint: disable=import-outside-toplevel
    import unyt

    # for a cylindrical dataset, ds.all_data().argmax(f) returns a tuple of
    # unyt_quantity objects with units (code_length, code_length, dimensionless),
    # which ds.point() doesn't like, so convert everything to floats
    def to_value(p) -> float:
        if isinstance(p, unyt.array.unyt_quantity):
            if p.units.is_dimensionless:
                return p.to_value()
            return p.to_value(ds.units.code_length)
        return float(p)

    pos = list(map(to_value, pos))
    # ds.point() requires a position with 3 coordinate values, regardless of
    # the dimensionality
    if len(pos) < ds.dimensionality:
        msg = f"need at least {ds.dimensionality} coordinate values in pos"
        raise ValueError(msg)
    if len(pos) < 3:
        pos.extend(ds.domain_center[len(pos):].to_value(ds.units.code_length))
    return ds.point(pos)


def to_conditions(point):
    """Extract the thermodynamic conditions needed to evaluate a network from a
    :py:class:`~yt.data_objects.selection_objects.point.YTPoint` object.

    Returns
    -------
    :
        tuple of (rho, T, composition)

    """
    rho = point["density"][0].to_value("g/cm**3")
    T = point["temperature"][0].to_value("K")
    comp = Composition([])
    pat = re.compile(r"X\(([a-zA-Z]+\d+)\)")
    for ns, field in point.ds.field_list:
        m = pat.match(field)
        if m is not None:
            comp[m[1]] = point[ns, field][0].to_value()

    return rho, T, comp
