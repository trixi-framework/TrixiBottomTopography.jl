module TestRBF1D

using Test
using TrixiBottomTopography
using KernelInterpolation: ThinPlateSplineKernel

@testset "1D RBF Interpolation" begin
    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    data = joinpath(root_dir, "test", "data", "rhine_1d_10_x_1.txt")

    # Define RBF interpolation structure
    itp = RBFInterpolation1D(data, ThinPlateSplineKernel{1}())
    # Define RBF interpolation function
    itp_func(x) = itp(x)

    # Test function at arbitrary point
    @test itp_func(357555) ≈ 46.19
end

end # module
