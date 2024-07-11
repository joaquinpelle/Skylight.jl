"""
    LamppostCorona <: AbstractCorona

Lamppost corona model as described in [Dauser et al. (2013)](https://academic.oup.com/mnras/article/430/3/1694/978006).

# Fields
- `height::Float64`: The height of the lamppost corona. 
- `spectral_index::Float64`: The spectral index of the power-law photon emission. Default is `2.0`.
- `theta_offset::Float64`: The offset of the lamppost corona in the polar angle. Default is `1e-8`. Must be small but nonzero to avoid the coordinate singularity at the polar axis.

# Examples
```julia
corona = LamppostCorona(height=2.5, spectral_index = 2.0)
```
"""
@with_kw struct LamppostCorona <: AbstractCorona
    height::Float64
    theta_offset::Float64 = 1e-8
    spectral_index::Float64 = 2.0
    @assert height>0.0 "height must be positive"
end

stationarity(::LamppostCorona) = IsStationary()
axisymmetry(::LamppostCorona) = IsAxisymmetric()

isvacuum(::LamppostCorona) = Vacuum()

function space_positions(npoints, spacetime, model::LamppostCorona, ::SphericalTopology, cache)
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
    corona::LamppostCorona;
    nbins)
    at_source = map(ray -> is_final_position_at_source(ray[1:4], spacetime, disk) && ray[3] ≈ π/2 && abs(Skylight.norm_squared(ray[5:8], metric(ray[1:4], spacetime))) < 1e-2, eachcol(output_data))
    radii = output_data[2,at_source]
    q = energies_quotients(output_data[:,at_source], spacetime, disk)
    # bins = radial_bins(disk, nbins=100)
    edges = range(cbrt(disk.inner_radius), stop=cbrt(disk.outer_radius), length=nbins+1).^3
    centers = midpoints(edges)
    A = equatorial_ring_areas(edges, spacetime)
    positions = (hcat∘map)(r -> equatorial_position(r, coordinates_topology(spacetime)), centers)
    γ = lorentz_factors(positions, spacetime, disk)
    h = StatsBase.fit(Histogram, radii, edges)
    h = LinearAlgebra.normalize(h, mode=:probability)
    𝓝 = h.weights
    qavg = average_inside_bins(q, radii, edges)
    Γ = corona.spectral_index
    n = 𝓝./(A.*γ)
    ϵ = qavg.^Γ.*n
    interp = LinearInterpolation(ϵ, centers; extrapolate=true)
    return interp(edges), edges
end