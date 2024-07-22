using Skylight
using CairoMakie
using DelimitedFiles

#Plot fonts
set_theme!(; fonts = (; regular = "Times New Roman"))

struct Resolutions
    npixels::Int
    num_bins::Int
    reltol::Float64
    abstol::Float64
end

struct ResolutionsSet
    npixels::Vector{Int}
    num_bins::Vector{Int}
    reltol::Vector{Float64}
    abstol::Vector{Float64}
end

function highest_resolution(resolutions_set::ResolutionsSet)
    return Resolutions(resolutions_set.npixels[end],
        resolutions_set.num_bins[end],
        resolutions_set.reltol[end],
        resolutions_set.abstol[end])
end

function index_npixels(resolutions_set::ResolutionsSet, idx::Int) 
    return Resolutions(resolutions_set.npixels[idx],
        resolutions_set.num_bins[end],
        resolutions_set.reltol[end],
        resolutions_set.abstol[end])
end

function index_num_bins(resolutions_set::ResolutionsSet, idx::Int) 
    return Resolutions(resolutions_set.npixels[end],
        resolutions_set.num_bins[idx],
        resolutions_set.reltol[end],
        resolutions_set.abstol[end])
end

function index_reltol(resolutions_set::ResolutionsSet, idx::Int) 
    return Resolutions(resolutions_set.npixels[end],
        resolutions_set.num_bins[end],
        resolutions_set.reltol[idx],
        resolutions_set.abstol[end])
end

function index_abstol(resolutions_set::ResolutionsSet, idx::Int) 
    return Resolutions(resolutions_set.npixels[end],
        resolutions_set.num_bins[end],
        resolutions_set.reltol[end],
        resolutions_set.abstol[idx])
end

function line_broadening_precision_test(resolutions_set::ResolutionsSet, rotation_sense::AbstractRotationSense)

    spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0, a=0.5)
    rISCO = isco_radius(spacetime, rotation_sense)
    model = DummyDisk(inner_radius = rISCO, outer_radius = 15.0, rotation_sense = rotation_sense)

    #TODO implement timing

    ref_resolutions = highest_resolution(resolutions_set)
    ref_binned_fluxes, ref_bins = line_broadening_precision_test(spacetime, model, ref_resolutions) 

    for i in 1:length(resolutions_set.npixels)
        resolutions = index_npixels(resolutions_set, i)
        binned_fluxes, bins = line_broadening_precision_test(spacetime, model, resolutions) 

        if i == 1
            max_diff = maximum(abs.(binned_fluxes - ref_binned_fluxes))
        else
            max_diff = max(max_diff, maximum(abs.(binned_fluxes - ref_binned_fluxes)))
        end
    end
    # We calculate midpoints of x to use as x coordinates for y
    bins_midpoints = midpoints(bins)
    return bins_midpoints, binned_fluxes, dexter 
end

function line_broadening_precision_test(spacetime, model, resolutions::Resolutions)
    
    npixels = resolutions.npixels
    num_bins = resolutions.num_bins
    reltol = resolutions.reltol
    abstol = resolutions.abstol

    camera = ImagePlane(distance = 500.0,
        observer_inclination_in_degrees = 30,
        horizontal_side = 31.5,
        vertical_side = 31.5,
        horizontal_number_of_pixels = npixels,
        vertical_number_of_pixels = npixels)

    configurations = VacuumOTEConfigurations(spacetime = spacetime,
        camera = camera,
        radiative_model = model,
        unit_mass_in_solar_masses = 1.0)

    initial_data = initialize(configurations)
    cb, cbp = callback_setup(configurations; rhorizon_bound = 2e-1)
    run = integrate(initial_data,
        configurations,
        cb,
        cbp;
        method = VCABM(),
        reltol = reltol,
        abstol = abstol)

    binned_fluxes, bins = line_emission_spectrum(initial_data,
        run.output_data,
        configurations;
        num_bins = num_bins,
        stop = 1.05)

    return binned_fluxes, bins 
end

function test_line_broadening(figname, rotation_sense)
    spacetime = KerrSpacetimeBoyerLindquistCoordinates(M = 1.0, a = 0.5)

    camera = ImagePlane(distance = 500.0,
        observer_inclination_in_degrees = 30,
        horizontal_side = 31.5,
        vertical_side = 31.5,
        horizontal_number_of_pixels = 120,
        vertical_number_of_pixels = 120)

    rISCO = isco_radius(spacetime, rotation_sense)
    model = DummyDisk(inner_radius = rISCO,
        outer_radius = 15.0,
        rotation_sense = rotation_sense)

    configurations = VacuumOTEConfigurations(spacetime = spacetime,
        camera = camera,
        radiative_model = model,
        unit_mass_in_solar_masses = 1.0)

    initial_data = initialize(configurations)
    cb, cbp = callback_setup(configurations; rhorizon_bound = 2e-1)
    run = integrate(initial_data,
        configurations,
        cb,
        cbp;
        method = VCABM(),
        reltol = 1e-5,
        abstol = 1e-5)

    binned_fluxes, bins = line_emission_spectrum(initial_data,
        run.output_data,
        configurations;
        num_bins = 30,
        stop = 1.05)

    dexter = readdlm(dexter_filename(rotation_sense), ',', Float64)

    max_dex = maximum(dexter[:, 2])
    max_flux = maximum(binned_fluxes)

    # We calculate midpoints of x to use as x coordinates for y
    bins_midpoints = 0.5 * (bins[1:(end - 1)] + bins[2:end])

    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1],
        xlabel = L"E/E_0",
        ylabel = "Flux (arbitrary)",
        title = "Relativistic line broadening",
        titlefont = :regular)
    skl = lines!(ax,
        bins_midpoints,
        binned_fluxes / max_flux,
        linewidth = 3,
        color = :black)
    dex = scatter!(ax,
        dexter[:, 1],
        dexter[:, 2] / max_dex,
        linewidth = 2,
        marker = :utriangle,
        color = :red,
        markersize = 16)

    ax.titlesize = 22
    ax.xlabelsize = 22
    ax.ylabelsize = 22
    ax.xticklabelsize = 15
    ax.yticklabelsize = 15

    skl.label = "Skylight"
    dex.label = "Dexter & Agol"

    axislegend(; labelsize = 18, position = :lt)

    # Save the figure
    save(figname, fig; dpi = 300)
end

dexter_filename(::ProgradeRotation) = "dexter_a05_pro.txt"
dexter_filename(::RetrogradeRotation) = "dexter_a05_ret.txt"

test_line_broadening("linebroadening_pro.png", ProgradeRotation())
test_line_broadening("linebroadening_ret.png", RetrogradeRotation())
