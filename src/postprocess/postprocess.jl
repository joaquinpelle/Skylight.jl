include("cache.jl")
include("views.jl")
include("finalpositions.jl")
include("normalization.jl")
include("redshift.jl")
include("intensities.jl")
include("fluxes.jl")
include("spectrum.jl")
include("lineemission.jl")

energies_quotients(initial_data, output_data, configurations::VacuumOTEConfigurations; kwargs...) = energies_quotients(initial_data, output_data, configurations, configurations.camera; kwargs...)
observed_bolometric_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations; kwargs...) = observed_bolometric_intensities(initial_data, output_data, configurations, configurations.camera; kwargs...)
observed_specific_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations, energy::Number; kwargs...) = observed_specific_intensities(initial_data, output_data, configurations, configurations.camera, [energy]; kwargs...)
observed_specific_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations, energies::AbstractVector; kwargs...) = observed_specific_intensities(initial_data, output_data, configurations, configurations.camera, energies; kwargs...)
spectrum(initial_data, output_data, configurations::VacuumOTEConfigurations; Emin::Number, Emax::Number, NE::Integer, kwargs...) = spectrum(initial_data, output_data, configurations, configurations.camera, range(Emin, stop=Emax, length=NE); kwargs...)

"""
    line_emission_spectrum(initial_data, output_data, configurations::VacuumOTEConfigurations; 
                           emission_profile::Function, bin_size::Number=NaN, num_bins::Int=NaN,
                           start::Number=NaN, stop::Number=NaN)

Compute the binned intensity of a line emission spectrum.

# Arguments
- `initial_data`: Initial condition data.
- `output_data`: Output data from the radiation model.
- `configurations::VacuumOTEConfigurations`: Configuration parameters for the model.

# Keywords
- `emission_profile::Function`: User-defined function describing the emission profile.
- `bin_size::Union{Nothing, Number}=nothing`: Size of each bin. Either `bin_size` or `num_bins` must be specified.
- `num_bins::Union{Nothing, Number}=nothing`: Number of bins. Either `bin_size` or `num_bins` must be specified.
- `start::Union{Nothing, Number}=nothing`: Lower bound of the range to be binned. If unspecified, the minimum of the energy quotients will be used.
- `stop::Union{Nothing, Number}=nothing`: Upper bound of the range to be binned. If unspecified, the maximum of the energy quotients will be used.
- `observer_four_velocity::AbstractVector` (optional): The four-velocity of the observer. If not provided, a default static four-velocity will be used.
- `flux_direction::AbstractVector` (optional): The direction in which to measure the flu. If not provided, a default direction will be used.

# Returns
- `binned_fluxes`: Array of the binned intensity in each bin.
- `bins`: Array of the bin edges.

# Notes
`observer_four_velocity` and `flux_direction` are only accepted if `configurations.camera` is a `PinholeCamera`.
"""
line_emission_spectrum(initial_data, output_data, configurations::VacuumOTEConfigurations; kwargs...) = line_emission_spectrum(initial_data, output_data, configurations, configurations.camera; kwargs...)