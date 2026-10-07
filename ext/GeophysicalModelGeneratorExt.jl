module GeophysicalModelGeneratorExt

using GeophysicalModelGenerator
using GMT
using TrixiBottomTopography

import TrixiBottomTopography: geo_topo_impression, create_topography_data

"""
    geo_topo_impression(; resolution::String, lon_min::Float64, lon_max::Float64,
                        lat_min::Float64, lat_max::Float64)

Function to load topographical data of the region enclosed by the given longitude and
latitude bounds with
[GeophysicalModelGenerator.jl](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl)
and convert it into Cartesian coordinates.

# Arguments
- `resolution`: String specifying the resolution of the topography data, see the table below
- `lon_min`: Minimum longitude value for the area of interest
- `lon_max`: Maximum longitude value for the area of interest
- `lat_min`: Minimum latitude value for the area of interest
- `lat_max`: Maximum latitude value for the area of interest

# Returns
- `Topo`: Original topography data from GeophysicalModelGenerator
- `p`: Projection point used for the coordinate transformation, placed in the middle of the
  chosen area
- `Topo_Cart`: Converted topography data in Cartesian coordinates

# Notes
The following topography data sets are available.

| Dataset | Resolution | Description |
|:--------|:----------:|:------------|
| `"@earth_relief_01s"` | 1 arc sec | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS] |
| `"@earth_relief_03s"` | 3 arc sec | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS] |
| `"@earth_relief_15s"` | 15 arc sec | SRTM15+ [David Sandwell, SIO/UCSD] |
| `"@earth_relief_30s"` | 30 arc sec | SRTM30+ [Becker et al., 2009, SIO/UCSD] |
| `"@earth_relief_01m"` | 1 arc min | ETOPO1 Ice surface [NEIC/NOAA] |
| `"@earth_relief_02m"` | 2 arc min | ETOPO2v2 Ice surface [NEIC/NOAA] |
| `"@earth_relief_03m"` | 3 arc min | ETOPO1 after Gaussian spherical filtering (5.6 km fullwidth) |
| `"@earth_relief_04m"` | 4 arc min | ETOPO1 after Gaussian spherical filtering (7.5 km fullwidth) |
| `"@earth_relief_05m"` | 5 arc min | ETOPO1 after Gaussian spherical filtering (9 km fullwidth) |
| `"@earth_relief_06m"` | 6 arc min | ETOPO1 after Gaussian spherical filtering (10 km fullwidth) |
| `"@earth_relief_10m"` | 10 arc min | ETOPO1 after Gaussian spherical filtering (18 km fullwidth) |
| `"@earth_relief_15m"` | 15 arc min | ETOPO1 after Gaussian spherical filtering (28 km fullwidth) |
| `"@earth_relief_20m"` | 20 arc min | ETOPO1 after Gaussian spherical filtering (37 km fullwidth) |
| `"@earth_relief_30m"` | 30 arc min | ETOPO1 after Gaussian spherical filtering (55 km fullwidth) |
| `"@earth_relief_60m"` | 60 arc min | ETOPO1 after Gaussian spherical filtering (111 km fullwidth) |
"""
function geo_topo_impression(;
                             resolution::String,
                             lon_min::Float64,
                             lon_max::Float64,
                             lat_min::Float64,
                             lat_max::Float64,)

    # Specify the limits of the topography data based on the coordinates
    lon_mean = (lon_max + lon_min) / 2
    lat_mean = (lat_min + lat_max) / 2

    # Loading the topography data
    Topo = import_topo(lon = [lon_min, lon_max], lat = [lat_min, lat_max],
                       file = resolution) # here we load the topography data

    p = ProjectionPoint(Lon = lon_mean, Lat = lat_mean) # to use Cartesian coordinates we choose a projection point here

    Topo_Cart = convert2CartData(Topo, p) # here we use the projection point to convert the topography data to Cartesian coordinates

    return Topo, p, Topo_Cart
end

"""
    create_topography_data(; low_x::Float64, high_x::Float64, gridsize_x::Float64,
                           low_y::Float64, high_y::Float64, gridsize_y::Float64,
                           write_path::String, dataname::String, Topo, p)

Function to project topography data onto a structured Cartesian grid and write it to an
`xyz` file. It takes the output of [`geo_topo_impression`](@ref) and creates a regular grid
which can afterwards be converted with [`convert_geo_1d`](@ref) or [`convert_geo_2d`](@ref).

# Arguments
- `low_x`: Minimum x-coordinate for the grid [km]
- `high_x`: Maximum x-coordinate for the grid [km]
- `gridsize_x`: Grid spacing in x-direction [km]
- `low_y`: Minimum y-coordinate for the grid [km]
- `high_y`: Maximum y-coordinate for the grid [km]
- `gridsize_y`: Grid spacing in y-direction [km]
- `write_path`: Directory path where the output file should be saved
- `dataname`: Name of the output file
- `Topo`: Topography data from GeophysicalModelGenerator (output from `geo_topo_impression`)
- `p`: Projection point for the coordinate transformation (output from `geo_topo_impression`)

# Returns
- `xyz`: `NamedTuple` with the vectors `x`, `y` and `z` containing the projected coordinates
  and elevations (in meters)
- `Topo_Cart_orth`: CartData object containing the projected topography data

# Notes
- All distance inputs (`low_x`, `high_x`, etc.) are given in kilometers
- The returned `xyz` and the written file contain coordinates in meters
- The output file is space separated with the columns `x y z`
"""
function create_topography_data(;
                                low_x::Float64,
                                high_x::Float64,
                                gridsize_x::Float64,
                                low_y::Float64,
                                high_y::Float64,
                                gridsize_y::Float64,
                                write_path::String,
                                dataname::String,
                                Topo,
                                p,)
    values_x = collect(low_x:gridsize_x:high_x)
    values_y = collect(low_y:gridsize_y:high_y)

    Topo_Cart_orth = CartData(xyz_grid(values_x, values_y, 0)) # create a grid with the desired gridpoints

    Topo_Cart_orth = project_CartData(Topo_Cart_orth, Topo, p) # project the topography data on the Cartesian grid

    # Extract x, y, z as vectors. Note: you have to scale the values to meters
    x = convert.(Float64, vec(Topo_Cart_orth.x.val[:, :, 1])) .* 1000
    y = convert.(Float64, vec(Topo_Cart_orth.y.val[:, :, 1])) .* 1000
    z = convert.(Float64, vec(Topo_Cart_orth.z.val[:, :, 1])) .* 1000

    # Write the data to a file and save it in the specified directory with the specified name
    output_file = joinpath(write_path, dataname)
    open(output_file, "w") do file
        for i in eachindex(x, y, z)
            println(file, round(x[i], digits = 5), " ", round(y[i], digits = 5), " ",
                    round(z[i], digits = 5))
        end
    end
    return (; x, y, z), Topo_Cart_orth
end
end
