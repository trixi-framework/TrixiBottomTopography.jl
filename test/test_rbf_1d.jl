module TestRBF1D

using Test
using TrixiBottomTopography
using KernelInterpolation: ThinPlateSplineKernel
using Random

@testset "1D RBF Interpolation" begin
    # Define data path
    root_dir = pkgdir(TrixiBottomTopography)
    data = joinpath(root_dir, "test", "data", "rhine_1d_10_x_1.txt")
    x, y = TrixiBottomTopography.parse_txt_1D(data)
    # Test scattered data by scrambling the data
    Random.seed!(24601)
    v = shuffle(collect(1:length(x)))
    xx = x[v]
    yy = y[v]
    itp = RBFInterpolation1D(xx, yy, ThinPlateSplineKernel{1}())
    # Define RBF interpolation function
    itp_func(x) = itp(x)

    # Test function at arbitrary point
    # The system is relatively ill-conditioned, so we need a bit looser tolerance here
    @test isapprox(itp_func(357555), 46.191872246844014; atol = 1e-10)
end

end # module
