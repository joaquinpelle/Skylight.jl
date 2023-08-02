@with_kw struct LamppostCorona <: AbstractRadiativeModel
    height::Float64
    theta_offset::Float64 = 1e-8
    spectral_index::Float64 = 2.0
    @assert height > 0.0 "height must be positive"
end

function space_positions(npoints, model::LamppostCorona, ::SphericalTopology)
    space_pos = zeros(3, npoints)
    space_pos[1,:] .= model.height
    space_pos[2,:] .= model.theta_offset
    return space_pos
end

function emitter_four_velocity!(vector, position, metric, spacetime, ::LamppostCorona, coords_top)
    static_four_velocity!(vector, metric)
    return nothing
end