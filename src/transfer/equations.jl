equations(::NonVacuumConfigurations) = non_vacuum_equations!
equations(::VacuumConfigurations) = geodesic_equations!

function non_vacuum_equations!(du, u::Array{Float64,1}, p::NonVacuumCache, t)

    geodesic_equations!(du, u, p, t)
    transfer_equations!(du, u, p, t)

    return nothing
end

function geodesic_equations!(du, u::Array{Float64,1}, p, t)

    spacetime = p.spacetime
    cache = p.multi_thread[Threads.threadid()]

    point = cache.point
    vel = cache.velocity

    @inbounds begin 
        for i in 1:4
            point[i] = u[i]
            vel[i] = u[4+i]
        end
    end

    Γ = cache.christoffel
    fill!(Γ, 0.0)

    christoffel_cache = cache.christoffel_cache
    christoffel!(Γ, point, spacetime, christoffel_cache) 

    a = cache.acceleration
    fill!(a, 0.0)

    @inbounds begin 
        for i in 1:4
            for j in 1:4
                for k in 1:4
                    a[i] += -Γ[i,j,k]*vel[j]*vel[k]
                end
            end
        end

        for i in 1:4
            du[i] = vel[i]
            du[4+i] = a[i]
        end
    end

    return nothing
end

function transfer_equations!(du, u::Array{Float64,1}, p, t)

    model = p.model
    observation_energies = p.observation_energies
    NE = p.NE

    cache = p.multi_thread[Threads.threadid()]

    point = cache.point
    vel = cache.velocity

    ε  = cache.ε
    αε = cache.αε
    jε = cache.jε
    
    ε .= -observation_energies*vel[1]

    invariant_absorptivity!(αε, point, ε, model)
    invariant_emissivity!(jε, point, ε, model)

    @inbounds begin    
        for i in 1:NE
            du[8+i] = αε[i]
            du[8+NE+i] = jε[i]*exp(-u[8+i])
        end
    end

    return nothing
end