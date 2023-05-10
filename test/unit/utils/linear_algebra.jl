using Skylight, Test

@testset "Linear algebra" begin
    
    using Skylight

    ginv = zeros(4,4)
    g = [1 2 3 4 ; 2 5 6 7 ; 3 6 8 9 ; 4 7 9 10]
    Skylight.inverse_4x4_symmetric!(ginv,g)
    @test ginv ≈ inv(g) atol=1e-14

end