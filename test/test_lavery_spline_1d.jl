module TestLaverySpline1D

using Test
using TrixiBottomTopography
using JuMP
using HiGHS

@testset "Lavery spline 1D interpolation" begin
    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    data = joinpath(root_dir, "test", "data", "rhine_1d_10_x_1.txt")

    # Define Lavery-spline structure
    spline_struct = LaverySpline1D(data)
    # Define Lavery-spline interpolation function
    spline_func(x) = spline_interpolation(spline_struct, x)

    # Test function at arbitrary point
    @test spline_func(357555) ≈ 46.19
end

end # module
