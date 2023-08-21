module Skylight

using Base.Threads: @threads, @spawn, nthreads, threadid
using DataInterpolations
using DelimitedFiles
using ForwardDiff
using HDF5
using KeywordDispatch
using LinearAlgebra
using Optim
using Parameters
using Polynomials
using Random
using Reexport
using Roots
using SpecialFunctions
using StaticArrays

@reexport using OrdinaryDiffEq
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
    BosonStarSpacetime

export AbstractSpacetimeCache,
    AbstractChristoffelCache

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
    axial_angular_momentum,
    time_translation_killing_vector,
    spherical_rotation_killing_vectors,
    axial_rotation_killing_vector

export metric!,
    metric,
    metric_inverse!,
    metric_inverse,
    volume_element,
    allocate_christoffel_cache,
    christoffel!,
    christoffel,
    coordinates_topology

export radius,
    event_horizon_radius,
    isco_radius,
    mbco_radius,
    circular_geodesic_angular_speed,
    circular_geodesic_specific_angular_momentum,
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
    StarAcrossWormhole,
    VerticalScreen,
    LamppostCorona,
    IonTorus

export torus_specific_angular_momentum!,
    cusp_and_center_radius!,
    torus_potentials_at_center_and_surface!

export Bremsstrahlung,
    Synchrotron,
    InverseCompton,
    ThermalEmission,
    SynchrotronAndBremsstrahlung

export allocate_cache,
    rest_frame_four_velocity!,
    surface_differential!,
    rest_frame_bolometric_intensity,
    rest_frame_specific_intensity,
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
    callback_setup,
    integrate

export callback,
    callback_parameters

export is_final_position_at_observer,
    is_final_position_at_source,
    energies_quotients,
    observed_bolometric_intensities,
    observed_specific_intensities,
    fluxes!,
    spectrum,
    line_emission_spectrum,
    rescale_intensities_normalization_at_real_observer!,
    axes_ranges,
    grid_view

export contract,
    lower_index,
    raise_index,
    scalar_product,
    norm_squared,
    normalize!,
    normalize_timelike!,
    normalize_spacelike!,
    orthogonal_projection,
    cos_angle_between_vectors,
    cos_angle_between_null_vectors

export inv4x4,
    inv4x4sym!,
    det4x4,
    det4x4sym

export planck_function,
    planck_integral,
    thermal_emission_window_intensity,
    thermal_emission_specific_intensity,
    thermal_emission_bolometric_intensity

export cartesian_from_spherical,
    spherical_from_cartesian

export time_translation_generator,
    time_translation_generator!,
    rotation_generators,
    rotation_generators!,
    zaxis_rotation_generator,
    zaxis_rotation_generator!,
    static_four_velocity,
    static_four_velocity!,
    circular_motion_four_velocity,
    circular_motion_four_velocity!

export PhysicalConstants,
    Dimensions

export geometrized_to_CGS,
    CGS_to_geometrized,
    pc_to_cm,
    cm_to_pc,
    Hz_to_erg,
    erg_to_Hz,
    eV_to_erg,
    erg_to_eV,
    keV_to_erg,
    erg_to_keV,
    per_erg_to_per_eV,
    per_eV_to_per_erg,
    per_erg_to_per_keV,
    per_keV_to_per_erg

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
