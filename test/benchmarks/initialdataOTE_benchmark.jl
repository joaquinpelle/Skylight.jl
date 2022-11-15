using Skylight, BenchmarkTools, BenchmarkPlots, StatsPlots

suite = BenchmarkGroup()

spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

numbers = [100,200,400,800]

for n in numbers

    image_plane = ImagePlane(observer_distance = 500.0,
                        observer_inclination_in_degrees = 90,
                        horizontal_side_image_plane = 10.0,
                        vertical_side_image_plane = 10.0,
                        horizontal_number_of_nodes = n,
                        vertical_number_of_nodes = n)

    configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                            image_plane = image_plane,
                                            initial_times = [0.0])

    suite[string(n)] = @benchmarkable get_initial_data($configurations)
    
end

tune!(suite)
t = run(suite)

p = plot()

for n in numbers
    plot!(p, t[string(n)],yaxis=:log10)
end

plot(p)