# Test suite
using Test

@testset "compute_gradient_norm tests" begin
    z = [i + j for i in 1:10, j in 1:10]
    dx_spacing = 1.0
    dy_spacing = 1.0
    grad_norm = compute_gradient_norm(z, dx_spacing, dy_spacing)
    @test size(grad_norm) == size(z)
    @test all(grad_norm .>= 0)

    z = [sin(i) + cos(j) for i in 1:0.1:10, j in 1:0.1:10]
    dx_spacing = 0.1
    dy_spacing = 0.1
    grad_norm = compute_gradient_norm(z, dx_spacing, dy_spacing)
    @test size(grad_norm) == size(z)
    @test all(grad_norm .>= 0)
end