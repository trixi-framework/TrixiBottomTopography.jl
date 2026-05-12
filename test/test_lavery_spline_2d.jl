module TestLaverySpline2D

using Test
using TrixiBottomTopography
using MathOptInterface
using HiGHS
using Downloads: download

@testset "Lavery spline 2D interpolation" begin
    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    # This is a text file version of the data in `lavery_spline_2d.jl`
    # found in the `examples` folder.
    data = joinpath(root_dir, "test", "data", "lavery_2d_10.txt")

    # Define Lavery spline structure
    spline_struct = LaverySpline2D(data)
    # Define Lavery spline interpolation function
    spline_func(x, y) = spline_interpolation(spline_struct, x, y)

    # Test function at four arbitrary points.
    @test spline_func(1.5, 4.1) ≈ 1.514
    @test spline_func(1.75, 4.53) ≈ 2.6163978125
    @test spline_func(1.5, 4.53) ≈ 1.772473
    @test spline_func(1.25, 4.53) ≈ 0.9285481875

    # Larger dataset of the Monai topography with 300 x 200 points
    # as a stress test to help catch significant performance regressions
    # The baseline time for construction of approximately 0.886560 seconds
    # is taken from https://github.com/trixi-framework/TrixiBottomTopography.jl/pull/99
    # This value is multiplied by 3 for "safety" of different runners and other fluctuations.
    spline_bathymetry_file = download("https://gist.githubusercontent.com/andrewwinters5000/21255c980c4eda5294f91e8dfe6c7e33/raw/1afb73928892774dc3a902e0c46ffd882ef03ee3/monai_bathymetry_data.txt",
                                      joinpath(@__DIR__, "monai_bathymetry_data.txt"));
    timing = @timed LaverySpline2D(spline_bathymetry_file)
    @test timing.time < 3 # seconds
end

end # module
