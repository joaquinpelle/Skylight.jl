using Skylight
using BenchmarkTools

function time_intensities(; N, tasks_per_thread)
    spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
    camera = ImagePlane(distance = 200.0,
                            observer_inclination_in_degrees = 5,
                            horizontal_side = 20.0,
                            vertical_side = 20.0,
                            horizontal_number_of_pixels = N,
                            vertical_number_of_pixels = N)
    model = NovikovThorneDisk(inner_radius=6.0, outer_radius=18.0)
    configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                    camera = camera,
                                    radiative_model = model,
                                    unit_mass_in_solar_masses=1.0)
    initial_data = initialize(configurations)
    cb, cbp = callback_setup(configurations; rhorizon_bound=2e-1) #... or, define your own cb and cbp
    run = integrate(initial_data, configurations, cb, cbp; method=VCABM(), reltol=1e-5, abstol=1e-5)
    output_data = run.output_data
    println("N = $N")
    println("Serial")
    @time _ = Skylight.observed_bolometric_intensities_serial(initial_data, output_data, configurations, camera)
    println("Multithread")
    @time _ = observed_bolometric_intensities(initial_data, output_data, configurations, camera; tasks_per_thread=tasks_per_thread)
    println("Multithread tmap")
    @time _ = Skylight.observed_bolometric_intensities_tmap(initial_data, output_data, configurations, camera; tasks_per_thread=tasks_per_thread)
    return nothing
end

function test_intensities(; N, tasks_per_thread)
    spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
    camera = ImagePlane(distance = 200.0,
                            observer_inclination_in_degrees = 5,
                            horizontal_side = 20.0,
                            vertical_side = 20.0,
                            horizontal_number_of_pixels = N,
                            vertical_number_of_pixels = N)
    model = NovikovThorneDisk(inner_radius=6.0, outer_radius=18.0)
    configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                    camera = camera,
                                    radiative_model = model,
                                    unit_mass_in_solar_masses=1.0)
    initial_data = initialize(configurations)
    cb, cbp = callback_setup(configurations; rhorizon_bound=2e-1) #... or, define your own cb and cbp
    run = integrate(initial_data, configurations, cb, cbp; method=VCABM(), reltol=1e-5, abstol=1e-5)
    output_data = run.output_data
    println("N = $N")
    Iobs = observed_bolometric_intensities(initial_data, output_data, configurations, camera; tasks_per_thread = tasks_per_thread)
    Itmap = Skylight.observed_bolometric_intensities_tmap(initial_data, output_data, configurations, camera; tasks_per_thread = tasks_per_thread)
    all(Iobs .== Itmap) || error("Iobs != Itmap") 
end

function btime_intensities(; N, tasks_per_thread)
    spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
    camera = ImagePlane(distance = 200.0,
                            observer_inclination_in_degrees = 5,
                            horizontal_side = 20.0,
                            vertical_side = 20.0,
                            horizontal_number_of_pixels = N,
                            vertical_number_of_pixels = N)
    model = NovikovThorneDisk(inner_radius=6.0, outer_radius=18.0)
    configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                    camera = camera,
                                    radiative_model = model,
                                    unit_mass_in_solar_masses=1.0)
    initial_data = initialize(configurations)
    cb, cbp = callback_setup(configurations; rhorizon_bound=2e-1) #... or, define your own cb and cbp
    run = integrate(initial_data, configurations, cb, cbp; method=VCABM(), reltol=1e-5, abstol=1e-5)
    output_data = run.output_data
    println("N = $N")
    println("Serial")
    @btime _ = Skylight.observed_bolometric_intensities_serial($initial_data, $output_data, $configurations, $camera)
    println("Multithread")
    @btime _ = observed_bolometric_intensities($initial_data, $output_data, $configurations, $camera; tasks_per_thread = $tasks_per_thread)
    println("Multithread tmap")
    @btime _ = Skylight.observed_bolometric_intensities_tmap($initial_data, $output_data, $configurations, $camera; tasks_per_thread = $tasks_per_thread)
    return nothing
end

time_intensities(N=3, tasks_per_thread=2)
time_intensities(N=3, tasks_per_thread=2)
btime_intensities(N=3, tasks_per_thread=2)
btime_intensities(N=10, tasks_per_thread=2)
btime_intensities(N=50, tasks_per_thread=2)
btime_intensities(N=50, tasks_per_thread=1) 
btime_intensities(N=100, tasks_per_thread=2) 
time_intensities(N=600, tasks_per_thread=2) 
time_intensities(N=1000, tasks_per_thread=2)
test_intensities(N=200, tasks_per_thread=2)