@with_kw struct LamppostCorona <: AbstractCorona
    height::Float64
    theta_offset::Float64 = 1e-8
    spectral_index::Float64 = 2.0
    @assert height>0.0 "height must be positive"
end

stationarity(::LamppostCorona) = IsStationary()

function space_positions(npoints, model::LamppostCorona, ::SphericalTopology)
    space_pos = zeros(3, npoints)
    space_pos[1, :] .= model.height
    space_pos[2, :] .= model.theta_offset
    return space_pos
end

function rest_frame_four_velocity!(vector,
    position,
    metric,
    spacetime,
    ::LamppostCorona,
    coords_top)
    static_four_velocity!(vector, metric)
    return nothing
end


"""Assumes unit energy in the rest frame of the emitter"""
function emissivity_profile(output_data::AbstractMatrix,
    spacetime::AbstractSpacetime,
    disk::AbstractAccretionDisk,
    corona::LamppostCorona)
    at_source = map(ray -> is_final_position_at_source(ray[1:4], spacetime, disk) && ray[3] ≈ π/2 && abs(Skylight.norm_squared(ray[5:8], metric(ray[1:4], spacetime))) < 1e-2, eachcol(output_data))
    radii = output_data[2,at_source]
    q = energies_quotients(output_data[:,at_source], spacetime, disk)
    # bins = radial_bins(disk, nbins=100)
    bins = range(cbrt(disk.inner_radius), stop=cbrt(disk.outer_radius), length=nbins).^3
    centers = midpoints(bins)
    A = equatorial_ring_areas(bins, spacetime)
    positions = (hcat∘map)(r -> equatorial_position(r, coordinates_topology(spacetime)), radii)
    γ = lorentz_factors(positions, spacetime, disk)
    h = StatsBase.fit(Histogram, radii, bins)
    h = LinearAlgebra.normalize(h, mode=:probability)
    𝓝 = h.weights
    qavg = average_inside_bins(q, radii, bins)
    Γ = corona.spectral_index
    n = 𝓝./(A.*γ)
    ϵ = qavg.^Γ.*n
    return ϵ, centers
end