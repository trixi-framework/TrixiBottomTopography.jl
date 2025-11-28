module TrixiBottomTopographyGeophysicalModelGeneratorExt

using GeophysicalModelGenerator
using DataFrames
using TrixiBottomTopography

# Import the functions so they can be extended
import TrixiBottomTopography: geo_topo_impression, create_topography_data

"""
    geo_topo_impression(; resolution::String, lon_min::Float64, lon_max::Float64, 
                         lat_min::Float64, lat_max::Float64)

Function to load and convert topographical data from GeophysicalModelGenerator into cartesian 
coordinates.

# Arguments
- `resolution`: String specifying the resolution of topography data (e.g. "@earth_relief_01s")
- `lon_min`: Minimum longitude value for the area of interest
- `lon_max`: Maximum longitude value for the area of interest
- `lat_min`: Minimum latitude value for the area of interest
- `lat_max`: Maximum latitude value for the area of interest

# Returns
- `Topo`: Original topography data from GeophysicalModelGenerator
- `p`: Projection point used for coordinate transformation
- `Topo_Cart`: Converted topography data in cartesian coordinates

# Notes
Available resolutions for topography data:
# # | Dataset                 |   Resolution |   Description                                               |
# # |:----------------        | ------------ | ----------------------------------------------------------- |
# # | "@earth\\_relief\\_01s" |	1 arc sec 	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
# # | "@earth\\_relief\\_03s"	|   3 arc sec	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
# # | "@earth\\_relief\\_15s"	|  15 arc sec	 | SRTM15+ [David Sandwell, SIO/UCSD]                          |
# # | "@earth\\_relief\\_30s"	|  30 arc sec	 | SRTM30+ [Becker et al., 2009, SIO/UCSD]                     |
# # | "@earth\\_relief\\_01m"	|   1 arc min	 | ETOPO1 Ice surface [NEIC/NOAA]                              |
# # | "@earth\\_relief\\_02m"	|   2 arc min	 | ETOPO2v2 Ice surface [NEIC/NOAA]                            |
# # | "@earth\\_relief\\_03m"	|   3 arc min	 | ETOPO1 after Gaussian spherical filtering (5.6 km fullwidth)|
# # | "@earth\\_relief\\_04m"	|   4 arc min	 | ETOPO1 after Gaussian spherical filtering (7.5 km fullwidth)|
# # | "@earth\\_relief\\_05m"	|   5 arc min	 | ETOPO1 after Gaussian spherical filtering (9 km fullwidth)  |
# # | "@earth\\_relief\\_06m"	|   6 arc min	 | ETOPO1 after Gaussia30n spherical filtering (10 km fullwidth)|
# # | "@earth\\_relief\\_10m"	|  10 arc min	 | ETOPO1 after Gaussian spherical filtering (18 km fullwidth) |
# # | "@earth\\_relief\\_15m"	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (28 km fullwidth) |
# # | "@earth\\_relief\\_20m"	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (37 km fullwidth) |
# # | "@earth\\_relief\\_30m"	|  30 arc min	 | ETOPO1 after Gaussian spherical filtering (55 km fullwidth) |
# # | "@earth\\_relief\\_60m"	|  60 arc min	 | ETOPO1 after Gaussian spherical filtering (111 km fullwidth)|

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

    p = ProjectionPoint(Lon = lon_mean, Lat = lat_mean) # to use cartesian coordinates we choose a projection point here

    Topo_Cart = convert2CartData(Topo, p) # here we use the projection point to convert the topography data to cartesian coordinates

    return Topo, p, Topo_Cart
end

"""
    create_topography_data(; low_x::Float64, high_x::Float64, gridsize_x::Float64,
                            low_y::Float64, high_y::Float64, gridzize_y::Float64,
                            write_path::String, dataname::String, Topo, p)

Function to create a structured grid of topography data and project it onto cartesian coordinates.
Takes the output from `geo_topo_impression()` and creates a regular grid for further processing.

# Arguments
- `low_x`: Minimum x-coordinate for the grid [km]
- `high_x`: Maximum x-coordinate for the grid [km]
- `gridsize_x`: Grid spacing in x-direction [km]
- `low_y`: Minimum y-coordinate for the grid [km]
- `high_y`: Maximum y-coordinate for the grid [km]
- `gridzize_y`: Grid spacing in y-direction [km]
- `write_path`: Directory path where the output file should be saved
- `dataname`: Name of the output file 
- `Topo`: Topography data from GeophysicalModelGenerator (output from geo_topo_impression)
- `p`: Projection point for coordinate transformation (output from geo_topo_impression)

# Returns
- `df_xyz`: DataFrame containing the projected coordinates and elevation (in meters)
- `Topo_Cart_orth`: CartData object containing the projected topography data

# Notes
- All distance inputs (low_x, high_x, etc.) should be in kilometers
- The output DataFrame and saved file will contain coordinates in meters
- The output file format is space-separated with columns: x y z
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

    Topo_Cart_orth = project_CartData(Topo_Cart_orth, Topo, p) # project the topography data on the cartasian grid

    # combine x, y, z into a DataFrame, warning: you have to scale the values to meters
    df_xyz = DataFrame(x = convert.(Float64, vec(Topo_Cart_orth.x.val[:, :, 1])) .* 1000,
                       y = convert.(Float64, vec(Topo_Cart_orth.y.val[:, :, 1])) .* 1000,
                       z = convert.(Float64, vec(Topo_Cart_orth.z.val[:, :, 1])) .* 1000)

    data_dir = write_path

    # Write the data to a file and save it in the specified directory with the specified name
    output_file = joinpath(data_dir, dataname)
    open(output_file, "w") do file
        for row in eachrow(df_xyz)
            rounded_row = [round(value, digits = 5) for value in row]
            println(file, join(rounded_row, " "))
        end
    end
    return df_xyz, Topo_Cart_orth
end
end