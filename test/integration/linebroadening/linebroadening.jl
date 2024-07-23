using Pkg
Pkg.activate("../../../run/Dev")
using Skylight
using CairoMakie
using DelimitedFiles

#Plot fonts
set_theme!(; fonts = (; regular = "Times New Roman"))

dexter_filename(::ProgradeRotation) = "dexter_a05_pro.txt"
dexter_filename(::RetrogradeRotation) = "dexter_a05_ret.txt"

function data_log(rotation_sense, npixels, reltol, abstol)
    println("Running line broadening data with:")
    println("  - Rotation sense: $rotation_sense")
    println("  - Number of pixels: $npixels")
    println("  - Relative tolerance: $reltol")
    println("  - Absolute tolerance: $abstol")
end

function line_broadening_data(rotation_sense; npixels, reltol, abstol)

    data_log(rotation_sense, npixels, reltol, abstol)

    spacetime = KerrSpacetimeBoyerLindquistCoordinates(M = 1.0, a = 0.5)
    rotation_sense = ProgradeRotation()
    camera = ImagePlane(distance = 500.0,
        observer_inclination_in_degrees = 30,
        horizontal_side = 31.5,
        vertical_side = 31.5,
        horizontal_number_of_pixels = npixels,
        vertical_number_of_pixels = npixels)

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
        reltol = reltol,
        abstol = abstol)

    return initial_data, run, configurations
end

function line_broadening_plot(initial_data, run, configurations; num_bins)

    rotation_sense = configurations.radiative_model.rotation_sense
    npixels = configurations.camera.horizontal_number_of_pixels
    reltol = run.reltol
    abstol = run.abstol
    output_data = run.output_data

    binned_fluxes, bins = line_emission_spectrum(initial_data,
        output_data,
        configurations;
        num_bins = num_bins,
        stop = 1.05)

    dexter = readdlm(dexter_filename(rotation_sense), ',', Float64)

    max_dex = maximum(dexter[:, 2])
    max_flux = maximum(binned_fluxes)

    # We calculate midpoints of x to use as x coordinates for y
    bins_midpoints = 0.5 * (bins[1:(end - 1)] + bins[2:end])

    fig = Figure(size = (600, 400))
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
    rstring = rotation_sense == ProgradeRotation() ? "prograde" : "retrograde"
    figname = "plots/$(rstring)_npx$(npixels)_reltol$(reltol)_abstol$(abstol)_nbins$(num_bins).png"
    save(figname, fig; dpi = 300)
end

function line_broadening_benchmark(rotation_sense)
    for npixels in [200, 400, 600, 800, 1000, 1200]
        for tol in [1e-6, 1e-10, 1e-14]
            initial_data, run, configurations = line_broadening_data(rotation_sense; npixels = npixels, reltol = tol, abstol = tol)
            for num_bins in [30, 40, 50, 60, 70, 80, 90, 100]
                line_broadeening_plot(initial_data, run, configurations; num_bins = num_bins)
            end
        end
    end
end

# test_line_broadening(ProgradeRotation(), npixels = 200, reltol = 1e-6, abstol = 1e-6, num_bins = 30)
# test_line_broadening(RetrogradeRotation(), npixels = 200, reltol = 1e-6, abstol = 1e-6, num_bins = 30)

line_broadening_benchmark(ProgradeRotation())