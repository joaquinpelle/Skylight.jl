using Skylight
using FiniteDiff
using BenchmarkTools

function finite_difference_metric_jacobian(position, spacetime::AbstractSpacetime, cache)
    return FiniteDiff.finite_difference_jacobian(Skylight.metric_field(spacetime),
        position,
        cache)
end

function finite_difference_metric_jacobian!(∂g,
    position,
    spacetime::AbstractSpacetime,
    cache)
    FiniteDiff.finite_difference_jacobian!(∂g,
        Skylight.metric_field(spacetime),
        position,
        cache)
    return nothing
end

function jacobian_finitediff(position, spacetime)
    ∂g = zeros(16, 4)
    g = zeros(4, 4)
    cache = FiniteDiff.JacobianCache(zero(position), zero(g), zero(g), Val(:central))
    finite_difference_metric_jacobian!(∂g, position, spacetime, cache)

    ∂gad = zeros(4, 4, 4)
    Skylight.metric_jacobian!(∂gad, position, spacetime, g)
    return reshape(∂g, 4, 4, 4), ∂gad
end

function full_benchmark_finitediff(position, spacetime)
    ∂g = zeros(16, 4)
    g = zeros(4, 4)
    cache = FiniteDiff.JacobianCache(zero(position), zero(g), zero(g), Val(:central))
    bfd = @benchmark finite_difference_metric_jacobian!($∂g, $position, $spacetime, $cache)

    ∂gad = zeros(4, 4, 4)
    cache = Skylight.AutoDiffChristoffelCache(spacetime)

    spacetime_metric_field = cache.spacetime_metric_field
    cfg = cache.cfg

    bad = @benchmark Skylight.metric_jacobian!($∂gad,
        $position,
        $spacetime_metric_field,
        $g,
        $cfg)
    return bfd, bad
end

function benchmark_finitediff(position, spacetime)
    ∂g = zeros(16, 4)
    g = zeros(4, 4)
    cache = FiniteDiff.JacobianCache(zero(position), zero(g), zero(g), Val(:central))
    println("Finite differences")
    @btime finite_difference_metric_jacobian!($∂g, $position, $spacetime, $cache)

    ∂gad = zeros(4, 4, 4)
    cache = Skylight.AutoDiffChristoffelCache(spacetime)

    spacetime_metric_field = cache.spacetime_metric_field
    cfg = cache.cfg

    println("Automatic differentiation")
    @btime Skylight.metric_jacobian!($∂gad, $position, $spacetime_metric_field, $g, $cfg)
    return nothing
end

println(FiniteDiff.default_relstep(Val(:central), Float64))

position = [0.0, 5.0, 0.25π, 0.0]
spacetime = KerrSpacetimeBoyerLindquistCoordinates(1.0, 0.9)
∂g, ∂gad = jacobian_finitediff(position, spacetime)
println(all(.≈(∂g, ∂gad, rtol = 1e-9)))
benchmark_finitediff(position, spacetime)

position = [0.0, 15.0, 0.0, 5.0]
spacetime = SuperposedPNSpacetime(m = (0.5, 0.5), chi = (0.9, 0.9), b = 20.0)
∂g, ∂gad = jacobian_finitediff(position, spacetime)
println(all(.≈(∂g, ∂gad, atol = 3e-10)))
benchmark_finitediff(position, spacetime)

position = [0.0, 1000.0, 0.0, 5.0]
spacetime = SuperposedPNSpacetime(m = (0.5, 0.5), chi = (0.6, 0.6), b = 20.0)
∂g, ∂gad = jacobian_finitediff(position, spacetime)
println(all(.≈(∂g, ∂gad, atol = 3e-10)))
benchmark_finitediff(position, spacetime)
