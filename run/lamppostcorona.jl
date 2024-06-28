using Skylight
using CairoMakie
using Printf

function main(; height, npp)

    spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0, a=0.0)
    disk = DummyDisk(inner_radius = isco_radius(spacetime, ProgradeRotation()), outer_radius = 100.0)
    corona = LamppostCorona(height=height, theta_offset=1e-5, spectral_index = 2.0)
    configurations = VacuumETOConfigurations(spacetime=spacetime,
                                    radiative_model = corona,
                                    number_of_points=1,
                                    number_of_packets_per_point = npp, 
                                    max_radius = 110.0,
                                    unit_mass_in_solar_masses=1.0)
    initial_data = initialize(configurations)
    cbp = callback_parameters(spacetime, disk, configurations; rhorizon_bound=2e-3)
    cb = callback(spacetime, disk)
    sim = integrate(initial_data, configurations, cb, cbp; method=VCABM(), reltol=1e-5, abstol=1e-5)
    output_data = sim.output_data

    I, bins_midpoints = emissivity_profile(output_data, spacetime, disk, corona)
end

main(height=2.5, npp=5000)