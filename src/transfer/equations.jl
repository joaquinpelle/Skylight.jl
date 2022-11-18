equations(configurations::VacuumConfigurations) = geodesic_equations!
equations(configurations::NonVacuumConfigurations) = non_vacuum_equations!

function non_vacuum_equations!(du, u::Array{Float64,1}, p::NonVacuumCache, t)

    geodesic_equations!(du, u, p, t)
    transfer_equations!(du, u, p, t)

end

function geodesic_equations!(du, u::Array{Float64,1}, p, t)

    spacetime = p.spacetime
    cache = p.multi_thread[Threads.threadid()]

    point = cache.point
    vel = cache.velocity

    for i in 1:4
        point[i] = u[i]
        vel[i] = u[4+i]
    end

    Γ = cache.christoffel
    fill!(Γ, 0.0)

    christoffel_cache = cache.christoffel_cache
    set_christoffel!(Γ, point, spacetime, christoffel_cache) 

    a = cache.acceleration
    fill!(a, 0.0)

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

function transfer_equations!(du, u::Array{Float64,1}, p, t)

    model = p.model
    observed_energies = p.observed_energies
    NE = p.NE

    ε  = cache.ε
    αε = cache.αε
    jε = cache.jε
    
    ε .= -observed_energies*vel[1]

    set_invariant_absorptivity!(αε, point, ε, model)
    set_invariant_emissivity!(jε, point, ε, model)

    for i in 1:NE
 
        du[8+i] = αε[i]
        du[8+NE+i] = jε[i]*exp(-u[8+i])
 
    end

end