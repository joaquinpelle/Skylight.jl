using Skylight

a_LBS1 = [-0.12835651132329287, 0.013682816437869968, -0.0000481174714805307, 
      0.00005645401289538938, -4.950385183447073e-6, 1.2257668397193955e-6, 
      2.234495931121217, -0.286705001384064, 0.0535085065184157, 
      -0.002796436948776261, 0.00043953642633610485, 1.5216185242475638e-5, 
      -2.2791584430232886e-6, 6.116664949485882e-7]

b_LBS1 = [0.002169200852262319, -0.19546743355115126, 0.06287130575728342, 
      -0.010651955590611452, 0.0009203500154294379, -3.751481658539803e-5, 
      13.99470385461245, -3.910737098471579, 0.3926444282051113, 
      0.1007509282796355, -0.03817380739469051, 0.00609572731408874, 
      -0.0004946554674560391, 1.8733095859225056e-5]


spacetime = BosonStarSpacetime(a=a_LBS1,b=b_LBS1)

image_plane = ImagePlane(distance = 500.0,
                         observer_inclination_in_degrees = 45,
                         horizontal_side_image_plane = 90.0,
                         vertical_side_image_plane = 90.0,
                         horizontal_number_of_nodes = 50,
                         vertical_number_of_nodes = 50)

model = BosonStarAccretionDisk(inner_radius=9.89994, outer_radius=79.89994, temperature_file="TempLBS1.dat")
        
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                         image_plane = image_plane,
                                         observed_times = [0.0],
                                         radiative_model = model,
                                         unit_mass_in_solar_masses=1.0)

initial_data = get_initial_data(configurations)

cb, cb_params = get_callback_and_params(configurations) #... or, define your own cb and cb_params

output_data = integrate(initial_data, configurations, cb, cb_params; method=VCABM(), reltol=1e-13, abstol=1e-21)

bolometric_intensities = get_observed_bolometric_intensities(initial_data, output_data, configurations)


xs, ys = get_coordinate_arrays(configurations) 
zs = view_intensities_grid(output_data, configurations, E_idx=1)

fig = Figure(font = "CMU Serif") #resolution=(600,400)
ax = Axis(fig[1,1], xlabel=L"\alpha", ylabel=L"\beta", ylabelsize = 26, xlabelsize = 26) 
hmap = heatmap!(xs, ys, zs; colormap=:gist_heat, interpolate=true)
Colorbar(fig[:, end+1], hmap, label=L"I", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
colsize!(fig.layout, 1, Aspect(1, 1.0))
colgap!(fig.layout, 7)
CairoMakie.save("plot.png", fig)
