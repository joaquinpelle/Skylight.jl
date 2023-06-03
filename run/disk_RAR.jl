using Skylight
using CairoMakie
using Printf

M1 = 1.02e7 #Unit mass in solar masses
energies = [200] #energies in keV of the fermion
inclinations = [5, 25, 45, 65, 85] #inclinations in degrees of the observer
Nres = 10 #number of pixels per side of the image
srsat = 150 #side of the image plane in units of Rsat

rsat = 1.980475e+12
rsat = CGS_to_geometrized(rsat, Dimensions.length, M1=M1) 
rd_out = 1e3*rsat

for E in energies

    spacetime = RARSpacetime(data_dir = "io/$(E)keV")
    r_inf, r_sup = radial_bounds(spacetime)
    inner_radii = Dict("r0"=>2*r_inf, "rsat"=>rsat)

    for ξ in inclinations
        for (rname, rd_in) in inner_radii
            
            ξstr = string(@sprintf("%02d", ξ)) 
            Nstr = string(@sprintf("%04d", Nres))
            filename = "$(E)keV_ideg$(ξstr)_$(rname)_s$(srsat)rsat_N$(Nstr)"

            image_plane = ImagePlane(distance = 1e-5*r_sup,
                                    observer_inclination_in_degrees = ξ,
                                    horizontal_side = srsat*rsat,
                                    vertical_side = srsat*rsat,
                                    horizontal_number_of_pixels = Nres,
                                    vertical_number_of_pixels = Nres)

            model = RARDisk(inner_radius=rd_in, outer_radius=rd_out, alpha=0.5, M1=M1)
                    
            configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                                    image_plane = image_plane,
                                                    observed_times = [0.0],
                                                    radiative_model = model,
                                                    unit_mass_in_solar_masses=model.M1)

            initial_data = get_initial_data(configurations)

            cb, cb_params = get_callback_and_params(configurations) #... or, define your own cb and cb_params

            run = integrate(initial_data, configurations, cb, cb_params; method=VCABM(), reltol=1e-13, abstol=1e-21)

            save_to_hdf5("io/RAR/$(filename).h5", configurations, initial_data, run; mode="w")        
            
            output_data = run.output_data

            Iobs, q = observed_bolometric_intensities(initial_data, output_data, configurations)

            xs,ys = get_pixel_coordinates_vectors(configurations)

            zs = view_as_grid(Iobs, configurations)

            fig = Figure(font = "Times New Roman")
            ax = Axis(fig[1,1], xlabel=L"\alpha/(GM/c^2)", ylabel=L"\beta/(GM/c^2)", ylabelsize = 26, xlabelsize = 26) 
            hmap = heatmap!(xs, ys, zs/maximum(zs); colormap=:gist_heat, interpolate=true)
            Colorbar(fig[:, end+1], hmap, label=L"I \text{(arbitrary)}", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
            colsize!(fig.layout, 1, Aspect(1, 1.0))
            colgap!(fig.layout, 7)
            CairoMakie.save("plots/RAR/bolometric/$(filename).png", fig)
            
            λ_EHT_Apr17 = 0.13
            ε = PhysicalConstants.h*PhysicalConstants.c/λ_EHT_Apr17

            Iobs, q = observed_specific_intensities(initial_data, output_data, configurations, [ε])

            zs = view_as_grid(Iobs, configurations; E_idx=1)

            fig = Figure(font = "Times New Roman")
            ax = Axis(fig[1,1], xlabel=L"\alpha/(GM/c^2)", ylabel=L"\beta/(GM/c^2)", ylabelsize = 26, xlabelsize = 26) 
            hmap = heatmap!(xs, ys, zs/maximum(zs); colormap=:gist_heat, interpolate=true)
            Colorbar(fig[:, end+1], hmap, label=L"I \text{(arbitrary)}", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
            colsize!(fig.layout, 1, Aspect(1, 1.0))
            colgap!(fig.layout, 7)
            CairoMakie.save("plots/RAR/specific/$(filename).png", fig)

        end
    end
end