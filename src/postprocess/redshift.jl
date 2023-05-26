function get_energies_quotient(ki, kf, cache)
    q = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric) /
        scalar_product(kf, cache.emitter_four_velocity, cache.emitter_metric)
    return q
end