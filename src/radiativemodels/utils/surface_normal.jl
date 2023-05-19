function set_unit_surface_normal!(vector, position, metric, metric_inverse, model, coords_top)

    set_surface_differential!(vector, position, model, coords_top)
    vector .= raise_index(vector,metric_inverse)
    normalize_spacelike!(vector, metric)

end