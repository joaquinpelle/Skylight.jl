function line_broadening_data(rotation_sense; npixels, reltol, abstol)

    data_log(rotation_sense, npixels, reltol, abstol)

    spacetime = KerrSpacetimeBoyerLindquistCoordinates(M = 1.0, a = 0.5)
    
    camera = ImagePlane(distance = 500.0,
        observer_inclination_in_degrees = 30,
        horizontal_side = 31.5,
        vertical_side = 31.5,
        horizontal_number_of_pixels = npixels,
        vertical_number_of_pixels = npixels)

    rISCO = isco_radius(spacetime, rotation_sense)
    model = DexterAgolDisk(inner_radius = rISCO,
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

    println("num_bins: $num_bins")

    rotation_sense = configurations.radiative_model.rotation_sense
    npixels = configurations.camera.horizontal_number_of_pixels
    reltol = run.kwargs[:reltol]
    abstol = run.kwargs[:abstol]
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
    figname = figure_name(rotation_sense, npixels, reltol, abstol, num_bins)
    save(figname, fig; dpi = 300)
end

function line_broadening_test(rotation_sense; npixels, tols, num_bins)
    println("Running in $(Threads.nthreads()) threads")
    println()
    for npx in npixels
        for tol in tols
            initial_data, run, configurations = line_broadening_data(rotation_sense; npixels = npx, reltol = tol, abstol = tol)
            for nbins in num_bins
                line_broadening_plot(initial_data, run, configurations; num_bins = nbins)
            end
            test_end_log()
        end
    end
end

dexter_filename(::ProgradeRotation) = "dexter/a05_pro.txt"
dexter_filename(::RetrogradeRotation) = "dexter/a05_ret.txt"

figure_name(rotation_sense, npixels, reltol, abstol, num_bins) = "plots/$(rotation_string(rotation_sense))/npx$(npixels)_reltol$(reltol)_abstol$(abstol)_nbins$(num_bins).png"

rotation_string(::ProgradeRotation) = "prograde"
rotation_string(::RetrogradeRotation) = "retrograde"

function data_log(rotation_sense, npixels, reltol, abstol)
    println("Running line broadening data with:")
    println("  - Rotation sense: $rotation_sense")
    println("  - Number of pixels: $npixels")
    println("  - Relative tolerance: $reltol")
    println("  - Absolute tolerance: $abstol")
    println("")
end

function test_end_log()
    println("")
    println("-------------------------------------------")
    println("")
end