abstract type CallbackParameters end

# Neutron star hot spots

@with_kw struct NeutronStarHotSpotsCallbackParameters <: CallbackParameters

    rmax::Float64
    rmin::Float64

end

function get_cb_params(model::NeutronStarHotSpots,configurations) 

    rmin = model.star_radius
    rmax = get_rmax(configurations)

    return NeutronStarHotSpotsCallbackParameters(rmin=rmin, rmax=rmax)

end

#Black hole accretion disk

@with_kw struct BlackHoleAccretionDiskCallbackParameters <: CallbackParameters
    
    rmax::Float64
    rmin::Float64
    inner_radius::Float64
    outer_radius::Float64

end

function get_cb_params(model::BlackHoleAccretionDisk, configurations)
    
    inner_radius = model.inner_radius
    outer_radius = model.outer_radius

    rhorizon = event_horizon_radius(configurations.spacetime)
    rbound = model.rbound
    rmin = rhorizon + rbound
    rmax = get_rmax(configurations)
    
    return BlackHoleAccretionDiskCallbackParameters(inner_radius=inner_radius, outer_radius=outer_radius, rmin=rmin, rmax=rmax)

end

# Star across wormhole

@with_kw struct StarAcrossWormholeCallbackParameters <: CallbackParameters
    
    rmax::Float64
    l_center::Float64
    star_radius::Float64

end

function get_cb_params(model::StarAcrossWormhole, configurations) 

    l_center = model.l_center
    star_radius = model.star_radius
    rmax = get_rmax(configurations)

    return StarAcrossWormholeCallbackParameters(l_center=l_center, star_radius=star_radius, rmax=rmax)

end

# Dummy extended region

@with_kw struct DummyExtendedRegionCallbackParameters <: CallbackParameters
    
    τmax::Float64
    rmin::Float64
    rmax::Float64

end

function get_cb_params(model::DummyExtendedRegion, configurations) 

    τmax = configurations.τmax
    rmax = get_rmax(configurations)

    rhorizon = event_horizon_radius(configurations.spacetime)
    rbound = model.rbound
    rmin = rhorizon + rbound

    return DummyExtendedRegionCallbackParameters(τmax=τmax, rmin=rmin, rmax=rmax)

end

get_rmax(configurations::ETOConfigurations) = configurations.observer_distance

function get_rmax(configurations::OTEConfigurations) 
    
    d = configurations.image_plane.observer_distance
    hs = configurations.image_plane.horizontal_side_image_plane
    vs = configurations.image_plane.vertical_side_image_plane
    return 1.1*sqrt(d^2 + vs^2 + hs^2)

end