using Parameters
using Test

include("geometry.jl")
include("spacetimes.jl")
include("initialdataOTE.jl")

spacetime = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0,
                                                                                    a=0.5))

image_plane = ImagePlane(observer_distance = 500.0,
                         observer_inclination_in_degrees = 45,
                         horizontal_side_image_plane = 10.0,
                         vertical_side_image_plane = 10.0,
                         horizontal_number_of_nodes = 50,
                         vertical_number_of_nodes = 50)

configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                              image_plane = image_plane,
                                              initial_times = [0.0])

initial_data = initialize_OTE(configurations)