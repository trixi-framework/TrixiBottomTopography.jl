# Converting DGM data files

## Introduction

Often geographical data is provided in the form of `.xyz` files. These files organized geographical
data in three columns separated by `[space]` with an `x`, `y` and `z` component. For example:
```
1.0 1.0 1.0
2.0 1.0 2.0
3.0 1.0 3.0
1.0 2.0 4.0
2.0 2.0 5.0
3.0 2.0 6.0
1.0 3.0 7.0
2.0 3.0 8.0
3.0 3.0 9.0
```
In the example above, we have `1.0, 2.0, 3.0` as possible `x` and `y` coordinate values.
Such `.xyz` formatted files provide the corresponding `z` value for all possible `x-y` combinations.
In this case, for example, at `[2.0, 3.0]`, the corresponding `z` value is `6.0`.

## DGM data set

In the [examples folder](https://github.com/trixi-framework/TrixiBottomTopography.jl/tree/main/examples) of this repo, the underlying data has been received from
[Geobasis NRW](https://www.bezreg-koeln.nrw.de/brk_internet/geobasis/hoehenmodelle/digitale_gelaendemodelle/gelaendemodell/index.html).
They provide a [geographical data set](https://www.opengeodata.nrw.de/produkte/geobasis/hm/dgm1_xyz/dgm1_xyz/) of the whole German state of North Rhine-Westphalia called **DGM**.
This data set contains patches of $1\,km^2$ where each patch has the elevation data for $1,000,000$ data points equally distributed as a grid with a grid size of $1\,m$. The data is given as `.xyz` files:
```
357000.00 5646999.00 47.40
357001.00 5646999.00 47.43
357002.00 5646999.00 47.49
357003.00 5646999.00 47.47
357004.00 5646999.00 47.39
357005.00 5646999.00 47.30
357006.00 5646999.00 47.24
...       ...        ...
```
where the first column provides the corresponding ETRS89 East coordinates, the second column the ETRS89 North coordinates and the third column the DHHN2016 height.

## Data format of `TrixiBottomTopography.jl`

The provided `.xyz` files of DGM are not directly accepted by `TrixiBottomTopography.jl` to
define B-spline interpolation structures. To make the raw topography data work with the
package, it must be converted into `.txt` files and organized in a specific format.

For one dimensional interpolation, `TrixiBottomTopography.jl` requires the following form:
```
# Number of x values
n
# x values
x_1
...
x_n
# y values
y_1
...
y_n
```
and to interpolate two dimensional data:
```
# Number of x values
n
# Number of y values
m
# x values
x_1
...
x_n
# y values
y_1
...
y_m
# z values
z_1,1
z_1,2
...
z_1,n
z_2,1
...
z_m,n
```
The `x, y` and `z` values must be set to `Float64` format.

## Conversion functions

To make matters easier, `TrixiBottomTopography.jl` provides the functions [`convert_dgm_1d`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.convert_dgm_1d-Tuple{String,%20String})
for one dimensional interpolation and [`convert_dgm_2d`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.convert_dgm_2d-Tuple{String,%20String}) for two dimensional interpolation.

To explain these functions, we consider the example file [`convert_data.jl`](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/convert_data.jl).

First, we import our package to be able to use the functions.

```julia
# Include packages
using TrixiBottomTopography
```

Next, we define the path of the `.xyz` file `path_src_file` that we want to convert,
as well as the paths of the files where we want to save them files.

```julia
# Get root directory
dir_path = pkgdir(TrixiBottomTopography)

# Define file paths
path_src_file = string(dir_path, "/examples/data/dgm1_32_357_5646_1_nw.xyz")
path_out_file_1d_x = string(dir_path, "/examples/data/rhine_data_1d_20_x.txt")
path_out_file_1d_y = string(dir_path, "/examples/data/rhine_data_1d_20_y.txt")
path_out_file_2d = string(dir_path, "/examples/data/rhine_data_2d_20.txt")
```

The source data from `path_src_file` looks as follows:

```
357000.00 5646999.00 47.40
357001.00 5646999.00 47.43
357002.00 5646999.00 47.49
357003.00 5646999.00 47.47
357004.00 5646999.00 47.39
357005.00 5646999.00 47.30
357006.00 5646999.00 47.24
...       ...        ...
```

[Here](https://gist.github.com/maxbertrand1996/c6917dcf80aef1704c633ec643a531d5), you can see the whole file.

Now the data can be converted.

```julia
# Convert data
convert_dgm_1d(path_src_file, path_out_file_1d_x; excerpt = 20, section = 100)
```

Calling this expression tells [`convert_dgm_1d`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.convert_dgm_1d-Tuple{String,%20String}) that the source file is `path_src_file` and the converted file will be saved in the file `path_out_file_1d_x`.
The optional attribute `excerpt` tells the function that only every `20`th point in the `x` direction (in this case, in the ETRS89 East coordinate) will be considered. Setting `section` to `100` tells the function that the corresponding `z` values (DHHN2016 in this case) from the `100`th `y` coordinate (ETRS89 North) will be chosen. The entire converted file produced by
this routine is available [here](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/data/rhine_data_1d_20_x.txt).

```julia
convert_dgm_1d(path_src_file, path_out_file_1d_y; excerpt = 20, direction = "y", section = 100)
```

Similar to the previous expression, this one has the additional attribute `directon = "y"` which tells [`convert_dgm_1d`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.convert_dgm_1d-Tuple{String,%20String}) that the data will be read from the `y` direction. (Click [here](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/data/rhine_data_1d_20_y.txt) to view the resulting file)

```julia
convert_dgm_2d(path_src_file, path_out_file_2d; excerpt = 20)
```
The two dimensional version [`convert_dgm_2d`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.convert_dgm_2d-Tuple{String,%20String}) works similar to the one dimensional case except that the optional attributes
`direction` and `section` do not exist, but only `excerpt`. Setting e.g. to `20` tells the
function that only every `20`th value in the `x` and `y` direction of the source file
`path_src_file` will be stored in `path_out_file_2d`. The resulting file after the two
dimensional conversion is available [here](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/data/rhine_data_2d_20.txt).