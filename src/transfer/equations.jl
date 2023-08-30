equations(configurations::AbstractConfigurations) = equations(isvacuum(configurations))
equations(::NonVacuum) = non_vacuum_equations
equations(::Vacuum) = geodesic_equations

function non_vacuum_equations(u::AbstractVector, p::NonVacuumCache, t)
    du = geodesic_equations(u, p, t)
    dτ, dI = transfer_equations(u, p, t)
    return vcat(du, dτ, dI)
end

function geodesic_equations(u::AbstractVector, p, t)
    spacetime = p.spacetime
    cache = p.multi_thread[Threads.threadid()]

    @inbounds begin
        @views begin
            position = u[1:4]
            momentum = u[5:8]
        end
    end

    Γ = cache.christoffel
    fill!(Γ, 0.0)
    christoffel_cache = cache.christoffel_cache
    christoffel!(Γ, position, spacetime, christoffel_cache)

    a = cache.acceleration
    fill!(a, 0.0)
    @inbounds begin
        for i in 1:4
            for j in 1:4
                for k in 1:4
                    a[i] += -Γ[i, j, k] * momentum[j] * momentum[k]
                end
            end
        end

        du1 = momentum[1]
        du2 = momentum[2]
        du3 = momentum[3]
        du4 = momentum[4]
        du5 = a[1]
        du6 = a[2]
        du7 = a[3]
        du8 = a[4]
    end
    return @SVector [du1, du2, du3, du4, du5, du6, du7, du8]
end

function transfer_equations(u::AbstractVector, p, t)
    spacetime = p.spacetime
    model = p.model
    coords_top = p.coordinates_topology
    observation_energies = p.observation_energies
    NE = p.NE
    cache = p.multi_thread[Threads.threadid()]

    @inbounds begin
        @views begin
            position = u[1:4]
            momentum = u[5:8]
            τε = u[9:(8 + NE)]
        end
    end

    ε = cache.ε
    αε = cache.αε
    jε = cache.jε
    vμ = cache.vμ
    metric = cache.metric
    spacetime_cache = cache.spacetime_cache

    metric!(metric, position, spacetime)
    rest_frame_four_velocity!(vμ, position, metric, spacetime, model, coords_top, spacetime_cache)
    rest_frame_energy = scalar_product(vμ, momentum, metric) #Without the negative sign because the momentum is past directed
    ε .= observation_energies * rest_frame_energy
    rest_frame_absorptivity!(αε, position, ε, metric, spacetime, model, coords_top)
    rest_frame_emissivity!(jε, position, ε, metric, spacetime, model, coords_top)
    return SVector{NE, Float64}(ε .* αε...),
    SVector{NE, Float64}(jε ./ (ε .^ 2) .* exp.(-τε)...)
end
