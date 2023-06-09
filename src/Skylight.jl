module Skylight

using DataInterpolations
using DelimitedFiles
using ForwardDiff
using HDF5
using KeywordDispatch
using LinearAlgebra
using Parameters
using Polynomials
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
    FRKerrSpacetime,
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

export stationarity,
    spherical_symmetry,
    axial_symmetry 

export energy,
    angular_momentum,
    z_angular_momentum

export metric!,
    metric_inverse!,
    volume_element, 
    allocate_christoffel_cache, 
    christoffel!,
    coordinates_topology

export radius,
    event_horizon_radius,
    isco_radius,
    circular_geodesic_angular_speed,
    radial_bounds,
    mass,
    mass_enclosed,
    mass_enclosed_derivative,
    horizons,
    horizon_parameter,
    cosmologic_horizon_radius

export AbstractRadiativeModel,
    AbstractSurfaceEmissionModel,
    DummyExtendedRegion,
    DummyModel,
    AbstractAccretionDisk,
    NovikovThorneDisk,
    ShakuraSunyaevDisk,
    RARDisk,
    AccretionDiskWithTabulatedTemperature,
    SyntheticPolarCap,
    StarAcrossWormhole

export emitter_four_velocity!,
    surface_differential!,
    emitted_bolometric_intensity,
    emitted_specific_intensity,
    fluxes,
    is_final_position_at_source,
    opaque_interior_surface_trait

export AbstractCamera,
    PinholeCamera,
    ImagePlane,
    AbstractConfigurations,
    NonVacuumOTEConfigurations,
    VacuumOTEConfigurations,
    VacuumETOConfigurations

export max_radius

export initialize, 
    integrate, 
    output_data, 
    callback_setup

export is_final_position_at_observer,
    energies_quotients, 
    observed_bolometric_intensities, 
    observed_specific_intensities,
    line_emission_spectrum, 
    rescale_intensities_normalization_at_real_observer!, 
    axes_ranges, 
    grid_view

export contract
    lower_index,
    raise_index,
    scalar_product,
    norm_squared,
    normalize_timelike!,
    normalize_spacelike!,
    orthogonal_projection,
    cos_angle_between_vectors,
    cos_angle_between_null_vectors
    
export inverse_4x4_symmetric!,
    determinant_4x4_symmetric

export planck_function,
    planck_integral,
    thermal_emission_window_intensity,
    thermal_emission_window_specific_intensity,
    thermal_emission_window_bolometric_intensity

export time_translation_generator,
    time_translation_generator!,
    rotation_generators,
    rotation_generators!,
    zaxis_rotation_generator,
    zaxis_rotation_generator!, 
    tangent_vector_zaxis_rotation!

export PhysicalConstants,
    Dimensions

export geometrized_to_CGS,
    CGS_to_geometrized

export TT, XX, YY, ZZ,
    RR, TH, PH,
    KTT, KXX, KYY, KZZ,
    KRR, KTH, KPH,
    LL, KLL

export save_to_hdf5,
    append_runs_to_hdf5,
    load_initial_data_from_hdf5,
    load_configurations_from_hdf5,
    load_runs_from_hdf5,
    load_output_data_from_hdf5,
    load_callback_params_from_hdf5,
    load_kwargs_from_hdf5,
    load_callback_from_hdf5
    
end
