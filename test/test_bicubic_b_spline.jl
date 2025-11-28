module TestBicubicBSpline

using Test
using TrixiBottomTopography

# Free end condition
@testset "Free end" begin

    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    data = joinpath(root_dir, "test", "data", "rhine_2d_10.txt")

    # Define B-spline structure
    spline_struct = BicubicBSpline(data)
    # Define B-spline interpolation function
    spline_func(x, y) = spline_interpolation(spline_struct, x, y)

    # Test function at arbitrary point
    @test 39.30464305919 < spline_func(357555, 5646555) < 39.30464305920
end

# Free end condition with smoothing
@testset "Free end + smoothing" begin

    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    data = joinpath(root_dir, "test", "data", "rhine_2d_10.txt")

    # Define B-spline structure
    spline_struct = BicubicBSpline(data; smoothing_factor = 9999)
    # Define B-spline interpolation function
    spline_func(x, y) = spline_interpolation(spline_struct, x, y)

    # Test function at arbitrary point
    @test 39.23237242 < spline_func(357555, 5646555) < 39.23237243
end

# Not-a-knot end condition
@testset "Not-a-knot end" begin

    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    data = joinpath(root_dir, "test", "data", "rhine_2d_10.txt")

    # Define B-spline structure
    spline_struct = BicubicBSpline(data; end_condition = "not-a-knot")
    # Define B-spline interpolation function
    spline_func(x, y) = spline_interpolation(spline_struct, x, y)

    # Test function at arbitrary point
    @test 39.30464305919 < spline_func(357555, 5646555) < 39.30464305920
end

# Not-a-knot end condition with smoothing
@testset "Not-a-knot end + smoothing" begin

    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    data = joinpath(root_dir, "test", "data", "rhine_2d_10.txt")

    # Define B-spline structure
    spline_struct = BicubicBSpline(data; end_condition = "not-a-knot",
                                   smoothing_factor = 9999)
    # Define B-spline interpolation function
    spline_func(x, y) = spline_interpolation(spline_struct, x, y)

    # Test function at arbitrary point
    @test 39.23237242 < spline_func(357555, 5646555) < 39.23237243
end

end # module