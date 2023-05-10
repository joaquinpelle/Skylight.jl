using Parameters, PreallocationTools

@testset "Automatic differentiation" begin
    
    M, a = 1.0, 0.5
    position = [ 2.64500435e+03, -7.59516667e+00, -7.59516667e+00,  1.35200090e+00]
    g = zeros(4,4)
    spacetime = KerrSpacetimeKerrSchildCoordinates(M=M, a=a)
    Γ = zeros(4,4,4)     
    cache = allocate_christoffel_cache(spacetime)
    set_christoffel!(Γ, position, spacetime, cache)

    #Define a Kerr spacetime to use with automatic differentiation

    struct ADKerrSpacetimeKerrSchildCoordinates

        M::Float64
        a::Float64
    
        @assert M >= 0.0
        @assert abs(a) <= M 
        
        ADKerrSpacetimeKerrSchildCoordinates(M,a) =  new(M,a)

    end

    function Skylight.set_metric!(g, position, spacetime::ADKerrSpacetimeKerrSchildCoordinates)

        M = spacetime.M
        a = spacetime.a
        
        t, x, y, z = position
        ρ2 = x^2 + y^2 + z^2
        a2 = a^2
        r2 = 0.5 * (ρ2 - a2) + sqrt(0.25 * (ρ2 - a2)^2 + a2 * z^2)
        r = sqrt(r2)
        H2 = 2. * M * r / (r2 + a2 * z^2 / r2)
        
        l1 = 1.
        l2 = (r*x + a*y)/(r2 + a2)
        l3 = (r*y - a*x)/(r2 + a2)
        l4 = z/r
        
        g[1,1]=-1. + H2 * l1*l1
        g[1,2]= 0. + H2 * l1*l2
        g[1,3]= 0. + H2 * l1*l3
        g[1,4]= 0. + H2 * l1*l4
        g[2,1]= g[1,2]
        g[2,2]= 1. + H2 * l2*l2
        g[2,3]= 0. + H2 * l2*l3
        g[2,4]= 0. + H2 * l2*l4
        g[3,1]= g[1,3]
        g[3,2]= g[2,3]
        g[3,3]= 1. + H2 * l3*l3
        g[3,4]= 0. + H2 * l3*l4
        g[4,1]= g[1,4]
        g[4,2]= g[2,4]
        g[4,3]= g[3,4]
        g[4,4]= 1. + H2 * l4*l4
        
        return nothing
        
    end
    
    spacetimeAD = ADKerrSpacetimeKerrSchildCoordinates(M, a)
    ΓAD = zeros(4,4,4) 
    #We didn't define set_christoffel! for this spacetime so it's going to use the default AutoDiff routines
    cacheAD = allocate_christoffel_cache(spacetimeAD)
    set_christoffel!(ΓAD, position, spacetimeAD, cacheAD)

    @test Γ ≈ ΓAD rtol=1e-15

end

