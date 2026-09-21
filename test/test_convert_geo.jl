module TestConvertGeoRoutines

using Test
using TrixiBottomTopography

# The `xyz` files created by `create_topography_data` contain one `x y z` triple per line,
# rounded to five digits. Check the conversion routines with some self generated data, so no data needs to be downloaded.
# The grid is deliberately rectangular (`nx != ny`) because this is what `convert_geo_1d`
# and `convert_geo_2d` add on top of `convert_dgm_1d` and `convert_dgm_2d`.
const NX = 9
const NY = 5

topography(x, y) = round(sin(x / 2) + cos(y / 2); digits = 5)

function write_synthetic_xyz(path)
    open(path, "w") do io
        for y in 0:(NY - 1), x in 0:(NX - 1)
            println(io, join((Float64(x), Float64(y), topography(x, y)), " "))
        end
    end
    return path
end

# Convenience wrappers around the parsing routines of TrixiBottomTopography which are also
# used by the B-spline constructors (make sure that the written files are readable by
# the package itself).
parse_1d(path) = TrixiBottomTopography.parse_txt_1D(path)
parse_2d(path) = TrixiBottomTopography.parse_txt_2D(path)

tmp_dir = mktempdir()
path_src_file = write_synthetic_xyz(joinpath(tmp_dir, "synthetic_geo.xyz"))

@testset "Check two dimensional conversion routine" begin
    path_out_file = joinpath(tmp_dir, "synthetic_2d_1.txt")
    convert_geo_2d(path_src_file, path_out_file; nx = NX, ny = NY)

    @test isfile(path_out_file)

    # The header must report the rectangular dimensions of the grid.
    lines = readlines(path_out_file)
    @test parse(Int, lines[2]) == NX
    @test parse(Int, lines[4]) == NY

    x, y, z = parse_2d(path_out_file)
    @test x == collect(0.0:8.0)
    @test y == collect(0.0:4.0)

    # `parse_txt_2D` returns `z` with the shape `size(y)` by `size(x)`
    @test size(z) == (NY, NX)
    @test z[1, 1] == topography(0, 0)
    @test z[1, NX] == topography(8, 0)
    @test z[NY, 1] == topography(0, 4)
    @test z[NY, NX] == topography(8, 4)
    @test z == [topography(i, j) for j in 0:(NY - 1), i in 0:(NX - 1)]

    # The converted file must be usable by the two dimensional B-splines, which reproduce
    # the data at the interpolation knots.
    spline_struct = BicubicBSpline(path_out_file)
    spline_func(x, y) = spline_interpolation(spline_struct, x, y)
    @test spline_func(0.0, 0.0)≈topography(0, 0) atol=1e-10
    @test spline_func(3.0, 2.0)≈topography(3, 2) atol=1e-10
    @test spline_func(8.0, 4.0)≈topography(8, 4) atol=1e-10
end

@testset "Check two dimensional conversion routine with stride" begin
    path_out_file = joinpath(tmp_dir, "synthetic_2d_2.txt")
    convert_geo_2d(path_src_file, path_out_file; nx = NX, ny = NY, excerpt = 2)

    # Every second value is taken in both directions
    lines = readlines(path_out_file)
    @test parse(Int, lines[2]) == 5
    @test parse(Int, lines[4]) == 3

    x, y, z = parse_2d(path_out_file)
    @test x == collect(0.0:2.0:8.0)
    @test y == collect(0.0:2.0:4.0)
    @test z == [topography(i, j) for j in 0:2:(NY - 1), i in 0:2:(NX - 1)]
end

@testset "Check one dimensional conversion routine in x direction" begin
    path_out_file = joinpath(tmp_dir, "synthetic_1d_1_x_3.txt")
    convert_geo_1d(path_src_file, path_out_file; nx = NX, ny = NY, section = 3)

    @test isfile(path_out_file)
    @test parse(Int, readlines(path_out_file)[2]) == NX

    # The `x` values are the coordinates along the `x` direction and the `y` values are the
    # corresponding elevations of the third section in `y` direction.
    x, y = parse_1d(path_out_file)
    @test x == collect(0.0:8.0)
    @test y == [topography(i, 2) for i in 0:(NX - 1)]

    spline_struct = CubicBSpline(path_out_file)
    spline_func(x) = spline_interpolation(spline_struct, x)
    @test spline_func(0.0)≈topography(0, 2) atol=1e-10
    @test spline_func(5.0)≈topography(5, 2) atol=1e-10
end

@testset "Check one dimensional conversion routine in y direction" begin
    path_out_file = joinpath(tmp_dir, "synthetic_1d_1_y_4.txt")
    convert_geo_1d(path_src_file, path_out_file; nx = NX, ny = NY, direction = "y",
                   section = 4)

    @test isfile(path_out_file)
    @test parse(Int, readlines(path_out_file)[2]) == NY

    # The `x` values must be the coordinates along the `y` direction. 
    # The `y` values are the elevations of the fourth section in `x` direction.
    x, y = parse_1d(path_out_file)
    @test x == collect(0.0:4.0)
    @test y == [topography(3, j) for j in 0:(NY - 1)]

    spline_struct = CubicBSpline(path_out_file)
    spline_func(y) = spline_interpolation(spline_struct, y)
    @test spline_func(0.0)≈topography(3, 0) atol=1e-10
    @test spline_func(4.0)≈topography(3, 4) atol=1e-10
end

@testset "Check one dimensional conversion routine with stride" begin
    path_out_file_x = joinpath(tmp_dir, "synthetic_1d_2_x_1.txt")
    convert_geo_1d(path_src_file, path_out_file_x; nx = NX, ny = NY, excerpt = 2)

    x, y = parse_1d(path_out_file_x)
    @test parse(Int, readlines(path_out_file_x)[2]) == 5
    @test x == collect(0.0:2.0:8.0)
    @test y == [topography(i, 0) for i in 0:2:(NX - 1)]

    path_out_file_y = joinpath(tmp_dir, "synthetic_1d_2_y_1.txt")
    convert_geo_1d(path_src_file, path_out_file_y; nx = NX, ny = NY, excerpt = 2,
                   direction = "y")

    x, y = parse_1d(path_out_file_y)
    @test parse(Int, readlines(path_out_file_y)[2]) == 3
    @test x == collect(0.0:2.0:4.0)
    @test y == [topography(0, j) for j in 0:2:(NY - 1)]
end

@testset "Check error handling" begin
    path_out_file = joinpath(tmp_dir, "synthetic_invalid.txt")

    # Invalid grid dimensions
    @test_throws ArgumentError convert_geo_1d(path_src_file, path_out_file; nx = 0, ny = NY)
    @test_throws ArgumentError convert_geo_1d(path_src_file, path_out_file; nx = NX, ny = 0)
    @test_throws ArgumentError convert_geo_2d(path_src_file, path_out_file; nx = 0, ny = NY)
    @test_throws ArgumentError convert_geo_2d(path_src_file, path_out_file; nx = NX, ny = 0)

    # Dimensions which do not match the number of lines in the data file
    @test_throws ArgumentError convert_geo_1d(path_src_file, path_out_file; nx = NX,
                                              ny = NY + 1)
    @test_throws ArgumentError convert_geo_2d(path_src_file, path_out_file; nx = NX,
                                              ny = NY + 1)

    # Invalid direction
    @test_throws ArgumentError convert_geo_1d(path_src_file, path_out_file; nx = NX,
                                              ny = NY, direction = "z")

    # Invalid sections
    @test_throws ArgumentError convert_geo_1d(path_src_file, path_out_file; nx = NX,
                                              ny = NY, section = 0)
    @test_throws ArgumentError convert_geo_1d(path_src_file, path_out_file; nx = NX,
                                              ny = NY, section = NY + 1)
    @test_throws ArgumentError convert_geo_1d(path_src_file, path_out_file; nx = NX,
                                              ny = NY, direction = "y", section = NX + 1)
end

end # module
