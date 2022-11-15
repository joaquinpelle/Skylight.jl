using Skylight
using DifferentialEquations
using CairoMakie

spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0,a=0.9)

image_plane = ImagePlane(observer_distance = 500.0,
                         observer_inclination_in_degrees = 45,
                         horizontal_side_image_plane = 20.0,
                         vertical_side_image_plane = 20.0,
                         horizontal_number_of_nodes = 200,
                         vertical_number_of_nodes = 200)

model = Skylight.DummyExtendedRegion(rbound=0.3)

configurations = Skylight.NonVacuumOTEConfigurations(spacetime=spacetime,
                                   image_plane = image_plane,
                                   radiative_model = model,
                                   initial_times = [0.0],
                                   observed_energies = [1.0],
                                   τmax = 2.0)

initial_data = get_initial_data(configurations)

cb, cb_params = get_callback_and_params(configurations) #... or, define your own cb and cb_params

solver_options = SolverOptions(method=VCABM(), reltol=1e-13, abstol=1e-21, output_type = SaveEndstate())

output_data = integrate_transfer(initial_data, configurations, cb, cb_params, solver_options)



function intensities_on(image_plane, output_data, E_idx)

    NE = Int((size(output_data, 1)-8)/2)
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    return [output_data[8+NE+E_idx,j+(i-1)*Nβ] for i in 1:Nα, j in 1:Nβ]

end

function get_horizontal_coordinates(image_plane)

    sα = image_plane.horizontal_side_image_plane
    Nα = image_plane.horizontal_number_of_nodes

    return range(-0.5*sα, stop=0.5*sα; length=Nα)

end

function get_vertical_coordinates(image_plane)

    sβ = image_plane.horizontal_side_image_plane
    Nβ = image_plane.horizontal_number_of_nodes

    return range(-0.5*sβ, stop=0.5*sβ; length=Nβ)

end

xs = get_horizontal_coordinates(image_plane)
ys = get_vertical_coordinates(image_plane)
zs = intensities_on(image_plane, output_data, 1)

fig, ax, hm = heatmap(xs, ys, zs)
Colorbar(fig[:, end+1], hm)

CairoMakie.save("plot.png", fig)
