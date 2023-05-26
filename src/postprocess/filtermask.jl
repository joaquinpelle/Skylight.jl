# Take the output data array and the configurations, and return a boolean array 
# that will filter out the rays which are away from the source

function get_filter_mask(output_data, configurations::VacuumOTEConfigurations)
    spacetime = configurations.spacetime
    model = configurations.radiative_model

    Nrays = number_of_initial_conditions(configurations)
    filter_mask = zeros(Bool, Nrays)
    for i in 1:Nrays
        @views begin 
            position = output_data[1:4,i]
        end

        filter_mask[i] = is_final_position_at_source(position, spacetime, model)
    end
    return filter_mask
end

# Now the same for the ETO scheme

function get_filter_mask(output_data, configurations::VacuumETOConfigurations)
    spacetime = configurations.spacetime
    model = configurations.radiative_model

    Nrays = number_of_initial_conditions(configurations)
    filter_mask = zeros(Bool, Nrays)

    for i in 1:Nrays
        @views begin 
            position = output_data[1:4,i]
        end

        filter_mask[i] = is_final_position_at_observer(position, spacetime, model)
    end
    return filter_mask
end