# Include packages
using TrixiBottomTopography
using Test

# Start testing
@testset "TrixiBottomTopography.jl" begin

    # Test data conversion
    @testset "Conversion" begin
        include("test_convert.jl")
    end

    # Linear B-splines
    @testset "Linear B-spline interpolation" begin
        include("test_linear_b_spline.jl")
    end

    # Cubic B-splines
    @testset "Cubic B-spline interpolation" begin
        include("test_cubic_b_spline.jl")
    end

    # Bilinear B-splines
    @testset "Bilinear B-spline interpolation" begin
        include("test_bilinear_b_spline.jl")
    end

    # Bicubic B-splines
    @testset "Bicubic B-spline interpolation" begin
        include("test_bicubic_b_spline.jl")
    end

    # 1D RBF interpolation
    @testset "1D RBF interpolation" begin
        include("test_rbf_1d.jl")
    end

    # 2D RBF interpolation
    @testset "2D RBF interpolation" begin
        include("test_rbf_2d.jl")
    end

    # 1D Lavery splines
    @testset "1D Lavery spline interpolation" begin
        include("test_lavery_spline_1d.jl")
    end
end
