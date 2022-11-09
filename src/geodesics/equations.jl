function geodesic_equations!(du, u::Array{Float64,1}, p::GeodesicsCache, t)

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