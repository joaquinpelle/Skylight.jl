using Skylight, Test

@testset "Solver options" begin
    
    solver_options = SolverOptions(output_type = SaveEndstate())

    func1 = Skylight.output_func(SaveEndstate())
    func2 = Skylight.output_func(SaveSolution()) 
    @test func1([0,1],1) == (1, false)
    @test func2([0,1],1) == ([0,1], false)

end

@testset "Collect output" begin
    
    spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

    image_plane = ImagePlane(observer_distance = 500.0,
                            observer_inclination_in_degrees = 45,
                            horizontal_side_image_plane = 10.0,
                            vertical_side_image_plane = 10.0,
                            horizontal_number_of_nodes = 10,
                            vertical_number_of_nodes = 10)

    model = SyntheticPolarCap(
                            star_radius=5.0,
                            angular_speed = 0.05, 
                            misalignment_angle_in_degrees=90,
                            angular_radius_in_degrees=60, 
                            temperature=1.0)
            
    configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                    image_plane = image_plane,
                                    initial_times = [0.0,1.0],
                                    radiative_model = model)

    initial_data = get_initial_data(configurations)

    cb, cb_params = get_callback_and_params(configurations) #... or, define your own cb and cb_params

    solver_options = SolverOptions(method=VCABM(), reltol=1e-13, abstol=1e-21, output_type = SaveEndstate())

    ensembleprob = Skylight.set_ensemble_problem(initial_data, configurations, cb_params, solver_options)
    
    method = solver_options.method
    reltol = solver_options.reltol
    abstol = solver_options.abstol
        
    #Also consider EnsembleSplitThreads() for multinodes and EnsembleGPUArray() for GPU

    stats = @timed sim = solve(ensembleprob, method, reltol=reltol, abstol=abstol, callback = cb, EnsembleThreads(); trajectories = 100)

    Skylight.collect_output(sim, solver_options.output_type)
    
end
