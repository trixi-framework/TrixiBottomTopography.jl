module TestLaverySpline1D

using Test
using TrixiBottomTopography
using JuMP
using HiGHS

@testset "Lavery spline 2D interpolation" begin
    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    # This is a text file version of the data in `lavery_spline_2d.jl`
    # found in the `examples` folder.
    data = joinpath(root_dir, "test", "data", "lavery_2d_10.txt")

    # Define Lavery spline structure
    spline_struct = LaverySpline2D(data; lambda = 5.0)
    # Define Lavery spline interpolation function
    spline_func(x, y) = spline_interpolation(spline_struct, x, y)

    # Test function at arbitrary points.

    # Evaluation in first Sibson element
    @test spline_func(1.5, 4.1) ≈ 1.514

    # Evaluation in second Sibson element
    @test spline_func(1.75, 4.53) ≈ 2.615

    # Evaluation in third Sibson element
    @test spline_func(1.5, 4.53) ≈ 1.772473

    # Evaluation in fourth Sibson element
    @test spline_func(1.25, 4.53) ≈ 0.929946
end

end # module
