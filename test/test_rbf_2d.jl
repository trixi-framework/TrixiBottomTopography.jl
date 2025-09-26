module TestRBF2D

using Test
using TrixiBottomTopography
using KernelInterpolation: ThinPlateSplineKernel

@testset "2D RBF Interpolation" begin
    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    data = joinpath(root_dir, "test", "data", "rhine_2d_10.txt")

    # Define RBF interpolation structure
    itp = RBFInterpolation2D(data, ThinPlateSplineKernel{2}())
    # Define RBF interpolation function
    itp_func(x, y) = itp([x; y])

    # Test function at arbitrary point
    @test 39.285 <= itp_func(357555, 5646555) <= 39.2686
end

end # module
