using Skylight, Test

@testset "Callbacks" begin

    spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)
    
    rhor = 1.0+sqrt(0.19)

    image_plane = ImagePlane(observer_distance = 50.0,
                            observer_inclination_in_degrees = 45,
                            horizontal_side_image_plane = 30.0,
                            vertical_side_image_plane = 40.0,
                            horizontal_number_of_nodes = 50,
                            vertical_number_of_nodes = 50)
    
    rmax = 1.1*sqrt(5000)

    @testset "Neutron star hot spots" begin

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

        
        @test Skylight.get_rmax(configurations) ≈ rmax
        
        cb_params = Skylight.get_cb_params(model, configurations)

        @test cb_params.rmin == 5.0
        @test cb_params.rmax == rmax

        p = Skylight.allocate_geodesics_cache(spacetime, cb_params)
        prob = ODEProblem(Skylight.geodesic_equations!, zeros(8), (0.0,1.0), p)
        integrator = init(prob, VCABM(),save_everystep=false,dt=1.0)
        u = [rand(), 10.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0]
        

        @test Skylight.star_cartesian_coordinates_condition(u,0.0,integrator) ≈ (rmax^2-100)*75
        @test Skylight.star_spherical_coordinates_condition(u,0.0,integrator) ≈ (rmax-10)*5.0

        Skylight.star_affect!(integrator)

        @test integrator.sol.retcode == :Terminated

        get_callback_and_params(configurations) 

    end

    @testset "Black hole accretion disk" begin

        model = NovikovThorneDisk(inner_radius=6.0, outer_radius=15.0)

        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                   image_plane = image_plane,
                                   initial_times = [0.0,1.0],
                                   radiative_model = model)

        cb_params = Skylight.get_cb_params(model, configurations)

        rmin = rhor + 1e-6
        @test cb_params.rmin ≈ rmin
        @test cb_params.rmax ≈ rmax
        @test cb_params.inner_radius == 6.0
        @test cb_params.outer_radius == 15.0
        
        p = Skylight.allocate_geodesics_cache(spacetime, cb_params)
        prob = ODEProblem(Skylight.geodesic_equations!, zeros(8), (0.0,1.0), p)
        integrator = init(prob, VCABM(),save_everystep=false,dt=1.0)

        u = [rand(), 10.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0]
        set_u!(integrator, u) 
        out = [0.0,0.0]

        r_kerr = Skylight.get_kerr_radius(u[1:4], spacetime)

        Skylight.black_hole_accretion_disk_cartesian_coordinates_condition(out, u, 0.0, integrator)

        @test out ≈ [0.0, (rmax^2-r_kerr^2)*(r_kerr^2-rmin*rmin)]
        
        Skylight.black_hole_accretion_disk_spherical_coordinates_condition(out, u, 0.0, integrator)

        @test out ≈ [-π/2, (rmax-10.0)*(10.0-rmin)]
        
        Skylight.black_hole_accretion_disk_cartesian_coordinates_affect!(integrator, 1)

        @test integrator.sol.retcode == :Terminated

        integrator = init(prob, VCABM(),save_everystep=false,dt=1.0)
        set_u!(integrator, u) 

        Skylight.black_hole_accretion_disk_cartesian_coordinates_affect!(integrator, 2)

        @test integrator.sol.retcode == :Terminated

        integrator = init(prob, VCABM(),save_everystep=false,dt=1.0)
        set_u!(integrator, u) 

        Skylight.black_hole_accretion_disk_spherical_coordinates_affect!(integrator, 1)

        @test integrator.sol.retcode == :Terminated

        integrator = init(prob, VCABM(),save_everystep=false,dt=1.0)
        set_u!(integrator, u) 

        Skylight.black_hole_accretion_disk_spherical_coordinates_affect!(integrator, 2)

        @test integrator.sol.retcode == :Terminated
        
        get_callback_and_params(configurations) 

    end

    @testset "Star across wormhole" begin
        
        spacetime = ChargedWormholeSpacetimeRegularCoordinates(b0=1.0, Q=0.5)

        model = StarAcrossWormhole(l_center=10.0, star_radius=5.0)

        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                   image_plane = image_plane,
                                   initial_times = [0.0,1.0],
                                   radiative_model = model)

        cb_params = Skylight.get_cb_params(model, configurations)

        @test cb_params.rmax ≈ rmax
        @test cb_params.star_radius == 5.0
        @test cb_params.l_center == 10.0
        
        p = Skylight.allocate_geodesics_cache(spacetime, cb_params)
        prob = ODEProblem(Skylight.geodesic_equations!, zeros(8), (0.0,1.0), p)
        integrator = init(prob, VCABM(),save_everystep=false,dt=1.0)

        u = [rand(), 10.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0]
        set_u!(integrator, u) 
        out = [0.0,0.0]

        Skylight.star_across_wormhole_condition(out, u, 0.0, integrator)

        Skylight.star_across_wormhole_affect!(integrator, 1)
        @test integrator.sol.retcode != :Terminated

        get_callback_and_params(configurations) 

    end

end