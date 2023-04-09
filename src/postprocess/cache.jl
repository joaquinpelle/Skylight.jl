abstract type PostProcessCache end

@with_kw mutable struct OTEPostProcessCache <: PostProcessCache
    
    observer_metric::Matrix{Float64} = zeros(4,4)
    emitter_metric::Matrix{Float64} = zeros(4,4)
    observer_four_velocity::Vector{Float64} = zeros(4)
    emitter_four_velocity::Vector{Float64} = zeros(4)

end

function dump_metrics_and_emitter_four_velocity_in!(cache::OTEPostProcessCache, initial_position, final_position, spacetime, model, coord_system)

    set_metric!(cache.observer_metric, initial_position, spacetime) 
    set_metric!(cache.emitter_metric, final_position, spacetime)

    set_emitter_four_velocity!(cache.emitter_four_velocity, final_position, cache.emitter_metric, spacetime, model, coord_system)

end

function dump_observer_four_velocity_in!(cache::OTEPostProcessCache)
    cache.observer_four_velocity = ∂t()
end

function unpack_views(cache::OTEPostProcessCache)

    @views begin

        observer_metric = cache.observer_metric
        emitter_metric = cache.emitter_metric
        observer_four_velocity = cache.observer_four_velocity
        emitter_four_velocity = cache.emitter_four_velocity

    end

    return observer_metric, emitter_metric, observer_four_velocity, emitter_four_velocity

end