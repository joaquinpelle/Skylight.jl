function parametric_typename(instance)
    typeof(instance).name.wrapper
end

function getproperty_nosave(obj::AbstractSpacetime, field::Symbol)
    value = getfield(obj, field)
    if value isa NoSaveField
        return value.value
    else
        return value
    end
end

function getproperty_nosave(obj::AbstractRadiativeModel, field::Symbol)
    value = getfield(obj, field)
    if value isa NoSaveField
        return value.value
    else
        return value
    end
end

function getproperty_nosave(obj::AbstractCallbackParameters, field::Symbol)
    value = getfield(obj, field)
    if value isa NoSaveField
        return value.value
    else
        return value
    end
end