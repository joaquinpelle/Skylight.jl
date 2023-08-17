using Skylight
using CairoMakie
using Printf

M1 = 1.02e7 #Unit mass in solar masses
energies = [200] #energies in keV of the fermion
inclinations = [5] #inclinations in degrees of the observer
Nres = 100 #number of pixels per side of the image
srsat = 75 #side of the image plane in units of Rsat

inner_radii = Dict("r0"=>1.654835e+06, "rsat"=>1.980475e+12)

for (rname, rd_in) in inner_radii
    inner_radii[rname] = CGS_to_geometrized(rd_in, Dimensions.length, M1=M1)
end
    
rsat = inner_radii["rsat"]
rd_out = 1e3*rsat

spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
rd_in = isco_radius(spacetime)

for E in energies

    rar_spacetime = RARSpacetime(data_dir = "io/$(E)keV")
    r_inf, r_sup = radial_bounds(rar_spacetime)

    for ξ in inclinations
            
        ξstr = string(@sprintf("%02d", ξ)) 
        Nstr = string(@sprintf("%04d", Nres))
        filename = "$(E)keV_ideg$(ξstr)_s$(srsat)rsat_N$(Nstr)"

        camera = ImagePlane(distance = 1e-5*r_sup,
                                observer_inclination_in_degrees = ξ,
                                horizontal_side = srsat*rsat,
                                vertical_side = srsat*rsat,
                                horizontal_number_of_pixels = Nres,
                                vertical_number_of_pixels = Nres)

        model = Skylight.ShakuraSunyaevDisk(inner_radius=rd_in, outer_radius=rd_out, alpha=0.5, M1=M1)
                
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                                camera = camera,
                                                radiative_model = model,
                                                unit_mass_in_solar_masses=model.M1)

        initial_data = initialize(configurations)

        cb, cbp = callback_setup(configurations, rhorizon_bound=2e-1) #... or, define your own cb and cbp

        run = integrate(initial_data, configurations, cb, cbp; method=VCABM(), reltol=1e-13, abstol=1e-21)

        save_to_hdf5("io/schw/$(filename).h5", configurations, initial_data, run)        
        
        output_data = run.output_data
        Iobs = observed_bolometric_intensities(initial_data, output_data, configurations)

        xs,ys = axes_ranges(camera)

        zs = grid_view(Iobs, configurations)

        fig = Figure(font = "Times New Roman")
        ax = Axis(fig[1,1], xlabel=L"\alpha/(GM/c^2)", ylabel=L"\beta/(GM/c^2)", ylabelsize = 26, xlabelsize = 26) 
        hmap = heatmap!(xs, ys, zs/maximum(zs); colormap=:gist_heat, interpolate=true)
        Colorbar(fig[:, end+1], hmap, label=L"I \text{(arbitrary)}", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
        colsize!(fig.layout, 1, Aspect(1, 1.0))
        colgap!(fig.layout, 7)
        CairoMakie.save("plots/schw/bolometric/$(filename).png", fig)
        
        λ_EHT_Apr17 = 0.13
        ε = PhysicalConstants.h*PhysicalConstants.c/λ_EHT_Apr17

        Iobs = observed_specific_intensities(initial_data, output_data, configurations, [ε])

        zs = grid_view(Iobs, configurations; energy_index=1)

        fig = Figure(font = "Times New Roman")
        ax = Axis(fig[1,1], xlabel=L"\alpha/(GM/c^2)", ylabel=L"\beta/(GM/c^2)", ylabelsize = 26, xlabelsize = 26) 
        hmap = heatmap!(xs, ys, zs/maximum(zs); colormap=:gist_heat, interpolate=true)
        Colorbar(fig[:, end+1], hmap, label=L"I \text{(arbitrary)}", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
        colsize!(fig.layout, 1, Aspect(1, 1.0))
        colgap!(fig.layout, 7)
        CairoMakie.save("plots/schw/specific/$(filename).png", fig)

    end
end