using Skylight
using CairoMakie

function radiative_transfer(;ν, θ, Δθ, Np, Npp)
    spacetime = SchwarzschildSpacetimeSphericalCoordinates(M = 1.0)
    model = CircularHotSpot(
        star_radius_in_km = 12,
        spin_frequency_in_Hz = ν,
        center_colatitude_in_degrees = θ,
        angular_radius_in_radians = Δθ,
        M1 = 1.4,
        temperature_in_keV = 0.35)
    configurations = VacuumETOConfigurations(spacetime = spacetime,
        radiative_model = model,
        number_of_points = Np,
        number_of_packets_per_point = Npp, 
        max_radius = 750.0,
        unit_mass_in_solar_masses = 1.4)
    initial_data = initialize(configurations)
    cb, cbp = callback_setup(configurations) #... or, define your own cb and cbp
    run = integrate(initial_data,
        configurations,
        cb,
        cbp;
        method = VCABM(),
        reltol = 1e-13,
        abstol = 1e-21)
    output_data = run.output_data
    return initial_data, output_data, configurations
end

function create_skymap(initial_data, output_data, configurations; Nξ, Nϕ)
    observation_energies = [keV_to_erg(1.0)]
    phase, skymap = specific_flux_skymap(initial_data, 
        output_data, 
        configurations, 
        observation_energies, 
        Nθ = Nξ, 
        Nϕ = Nϕ)
    skymap ./= maximum(skymap)
    return phase, skymap[1, :, :]
end

function plot_light_curve(skymap; ξ, Δϕ=0.0)
    Nξ, Nϕ = size(skymap)
    ξbins = create_bins(num_bins = Nξ, start = 0.0, stop = 180)
    ϕbins = create_bins(num_bins = Nϕ, start = 0.0, stop = 1.0)
    phase = midpoints(ϕbins)
    idξ = searchsortedlast(ξbins, ξ)

    light_curve = skymap[idξ, :]

    Δn = round(Int, Δϕ*Nϕ)
    light_curve = circshift(light_curve, Δn)

    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1])
    lines!(ax, phase, light_curve, label = "ξ = $(ξ)°")
    ax.ylabel = "Flux at 1 keV"
    ax.xlabel = "Phase"
    ax.xlabelsize = 22
    ax.ylabelsize = 22
    ax.xticklabelsize = 18
    ax.yticklabelsize = 18
    ax.xtickalign = 1
    ax.ytickalign = 1
    ax.xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    ax.titlesize = 22
    leg = axislegend(ax, position = :ct, labelsize = 22)
    CairoMakie.save("light_curve.png", fig)
end

function plot_skymap(skymap)
    Nξ, Nϕ = size(skymap)
    ξbins = create_bins(num_bins = Nξ, start = 0.0, stop = 180)
    ϕbins = create_bins(num_bins = Nϕ, start = 0.0, stop = 1.0)
    phase = midpoints(ϕbins)
    inclinations = midpoints(ξbins)

    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1])
    heatmap!(ax, phase, inclinations, transpose(skymap), label = "Skymap")
    ax.ylabel = "Specific flux"
    ax.xlabel = "Phase"
    CairoMakie.save("skymap.png", fig)
end

initial_data, output_data, configurations = radiative_transfer(ν=200, θ=90, Δθ=0.01, Np=100, Npp=20000)
phase, skymap = create_skymap(initial_data, output_data, configurations; Nξ=40, Nϕ=60)
plot_light_curve(skymap; ξ=90)