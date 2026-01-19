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

    # Test function at arbitrary point
    @test spline_func(1.5, 4.1) ≈ 1.514
end

end # module
