module Skylight

using DataInterpolations
using DelimitedFiles
using ForwardDiff
using HDF5
using KeywordDispatch
using LinearAlgebra
using Parameters
using Random
using Reexport

@reexport using DifferentialEquations
@reexport using PreallocationTools

include("utils/types.jl")
include("spacetimes/types.jl")
include("radiativemodels/types.jl")
include("configurations/types.jl")
include("initialdata/types.jl")
include("transfer/types.jl")
include("postprocess/types.jl")

include("utils/utils.jl")
include("spacetimes/spacetimes.jl")
include("radiativemodels/radiativemodels.jl")
include("configurations/configurations.jl")
include("initialdata/initialdata.jl")
include("transfer/transfer.jl")
include("postprocess/postprocess.jl")

export AbstractSpacetime,
    AbstractBlackHoleSpacetime,
    AbstractRegularCompactObjectSpacetime,
    MinkowskiSpacetimeCartesianCoordinates,
    MinkowskiSpacetimeSphericalCoordinates,
    AbstractSchwarzschildSpacetime,
    SchwarzschildSpacetimeKerrSchildCoordinates, 
    SchwarzschildSpacetimeSphericalCoordinates,
    AbstractKerrSpacetime,
    KerrSpacetimeKerrSchildCoordinates,
    KerrSpacetimeBoyerLindquistCoordinates,
    JohannsenSpacetime,
    ChargedWormholeSpacetimeSphericalCoordinates,
    ChargedWormholeSpacetimeRegularCoordinates,
    RARSpacetime,
    BosonStarSpacetime,
    NumericalSpacetime

export AbstractChristoffelCache

export AbstractCoordinatesTopology, 
    CartesianTopology,
    SphericalTopology

export AbstractRotationSense,
    ProgradeRotation,
    RetrogradeRotation

export IsStationary,
    IsSphericallySymmetric,
    IsAxiallySymmetric

export is_stationary,
    is_spherically_symmetric,
    is_axially_symmetric 

export energy,
    angular_momentum,
    z_angular_momentum

export set_metric!,
    set_metric_inverse!,
    volume_element, 
    allocate_christoffel_cache, 
    set_christoffel!,
    coordinates_topology

export radius, 
    event_horizon_radius,
    isco_radius,
    circular_geodesic_angular_speed

export AbstractRadiativeModel,
    AbstractSurfaceEmissionModel,
    DummyExtendedRegion,
    DummyModel,
    AbstractAccretionDisk,
    NovikovThorneDisk,
    AccretionDiskWithTabulatedTemperature,
    SyntheticPolarCap,
    StarAcrossWormhole

export set_emitter_four_velocity!,
    set_surface_differential!,
    emitted_bolometric_intensity,
    emitted_specific_intensity,
    is_final_position_at_source

export AbstractCamera,
    PinholeCamera,
    ImagePlane,
    AbstractConfigurations,
    NonVacuumOTEConfigurations,
    VacuumOTEConfigurations,
    VacuumETOConfigurations

export get_initial_data, 
    integrate, 
    output_data, 
    get_callback_and_params

export is_final_position_at_observer,
    energies_quotients, 
    observed_bolometric_intensities, 
    observed_specific_intensities,
    line_emission_spectrum, 
    rescale_intensities_normalization_at_real_observer!, 
    get_pixel_coordinates_vectors, 
    view_as_grid

export save_to_hdf5,
    append_runs_to_hdf5,
    load_initial_data_from_hdf5,
    load_configurations_from_hdf5,
    load_runs_from_hdf5,
    load_output_data_from_hdf5,
    load_callback_params_from_hdf5,
    load_kwargs_from_hdf5,
    load_callback_from_hdf5
    
export tangent_vector_zaxis_rotation!

export geometrized_to_CGS,
    CGS_to_geometrized

export TT, XX, YY, ZZ,
    RR, TH, PH,
    KTT, KXX, KYY, KZZ,
    KRR, KTH, KPH,
    LL, KLL

end
