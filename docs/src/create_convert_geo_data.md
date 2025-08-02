# Create and convert real topography data
In this section, it is explained how to get real topography data using GeophysicalModelGenerator.jl and make it accessible to TrixiBottomTopography. The necessary structure of the data for TrixiBottomTopography is explained in the `docs/src/conversion.md` file. In the file `docs/src/trixishallowwater_jl_geo_data_examples`, it is then explained how to properly run simulations with this data. In the following, we will consider the example file `create_convert_data_geo.jl`.

## Getting an impression of the topography data

To get a first impression of the desired topography, TrixiBottomTopography.jl provides the function `geo_topo_impression`. First, we import our packages to be able to use all the functions.

```@example create_convert_geo
# Include packages
using GeophysicalModelGenerator
using TrixiBottomTopography
using GMT # needs to be installed for the functions geo_topo_impression and create_topography_data
```

Next, we can run `geo_topo_impression` to load and preview the topography data:

```@example create_convert_geo
# Get a first impression of the topography data
Topo, p, Topo_Cart = geo_topo_impression(
    resolution="@earth_relief_01s", 
    lon_min=6.963880, 
    lon_max=6.978499, 
    lat_min=50.947861, 
    lat_max=50.957095
)
```

This function requires you to:

- **Define the geographic region**: Specify the longitude and latitude coordinates for your area of interest. The values used here correspond to the example from `docs/src/conversion.md`.

- **Choose the data resolution**: Select from the available topography datasets. Higher resolution means more detailed data but larger file sizes.

The function returns three important objects:
- `Topo`: The raw topography data from GeophysicalModelGenerator.
- `p`: The projection point used for coordinate transformation. The projection point is choosen to be in the middle of the chosen area.
- `Topo_Cart`: The topography data converted to Cartesian coordinates.

<table style="border-collapse: collapse; width: 100%; border: 2px solid #666;">
    <thead>
        <tr style="border: 2px solid #666; background-color: #f3f3f3;">
            <th style="border: 2px solid #666; padding: 8px; text-align: left;">Dataset</th>
            <th style="border: 2px solid #666; padding: 8px; text-align: center;">Resolution</th>
            <th style="border: 2px solid #666; padding: 8px; text-align: left;">Description</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_01s"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">1 arc sec</td>
            <td style="border: 1px solid #666; padding: 8px;">SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_03s"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">3 arc sec</td>
            <td style="border: 1px solid #666; padding: 8px;">SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_15s"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">15 arc sec</td>
            <td style="border: 1px solid #666; padding: 8px;">SRTM15+ [David Sandwell, SIO/UCSD]</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_30s"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">30 arc sec</td>
            <td style="border: 1px solid #666; padding: 8px;">SRTM30+ [Becker et al., 2009, SIO/UCSD]</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_01m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">1 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 Ice surface [NEIC/NOAA]</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_02m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">2 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO2v2 Ice surface [NEIC/NOAA]</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_03m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">3 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 after Gaussian spherical filtering (5.6 km fullwidth)</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_04m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">4 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 after Gaussian spherical filtering (7.5 km fullwidth)</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_05m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">5 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 after Gaussian spherical filtering (9 km fullwidth)</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_06m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">6 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 after Gaussian spherical filtering (10 km fullwidth)</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_10m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">10 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 after Gaussian spherical filtering (18 km fullwidth)</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_15m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">15 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 after Gaussian spherical filtering (28 km fullwidth)</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_20m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">20 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 after Gaussian spherical filtering (37 km fullwidth)</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_30m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">30 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 after Gaussian spherical filtering (55 km fullwidth)</td>
        </tr>
        <tr>
            <td style="border: 1px solid #666; padding: 8px;">"@earth_relief_60m"</td>
            <td style="border: 1px solid #666; padding: 8px; text-align: center;">60 arc min</td>
            <td style="border: 1px solid #666; padding: 8px;">ETOPO1 after Gaussian spherical filtering (111 km fullwidth)</td>
        </tr>
    </tbody>
</table>

## Creating structured topography data

Based on the topography data from `geo_topo_impression`, we can now create a structured grid onto which the topography data is then interpolated by using the `create_topography_data` function.

```@example create_convert_geo
# Create a structured grid of topography data
df_xyz, Topo_Cart_orth = create_topography_data(
    low_x = -0.5, 
    high_x = 0.499, 
    gridsize_x = 0.001, 
    low_y = -0.5, 
    high_y = 0.499, 
    gridsize_y = 0.002, 
    write_path = joinpath(@__DIR__, "TrixiBottomTopography.jl/examples/data"), 
    dataname = "geo.xyz",
    Topo = Topo, 
    p = p
)
```

This function requires you to specify:

- **Spatial domain**: Define the boundaries of your simulation area using `low_x`, `high_x`, `low_y`, and `high_y` (in kilometers).
- **Grid resolution**: Set the spacing between grid points using `gridsize_x` and `gridsize_y` (in kilometers). This defines the resolution of the grid onto which the topography data is interpolated.
- **Output settings**: Choose the directory (`write_path`) and filename (`dataname`) for the resulting data file.
- **Input data**: Pass the `Topo` and `p` objects from the previous `geo_topo_impression` call.

The function returns:
- `df_xyz`: A DataFrame containing the projected coordinates and elevation data.
- `Topo_Cart_orth`: A CartData object with the topography data interpolated onto regular grid points.

**Important notes:**
- All input distance parameters are in **kilometers**.
- The output file contains coordinates and elevations in **meters**.
- The output file format is space-separated with columns: x y z.



For more detailed information about the coordinate projections and transformations used in these functions, refer to the [projection section of the GeophysicalModelGenerator.jl user guide](https://juliageodynamics.github.io/GeophysicalModelGenerator.jl/dev/man/projection/).


## Conversion functions

After creating the structured topography data with `create_topography_data`, we need to convert the `geo.xyz` file into the correct format required by TrixiBottomTopography. The required format is explained in `docs/src/conversion.md`. TrixiBottomTopography provides two specialized functions for this conversion using not only quadratic but also rectangular domains:

- `convert_geo_1d`: For one-dimensional interpolation along either x or y direction.
- `convert_geo_2d`: For two-dimensional interpolation preserving the full grid structure.

### Converting to 1D format

```@example create_convert_geo
# Define file paths
data_dir = joinpath(@__DIR__, "TrixiBottomTopography.jl/examples/data")
path_src_file = joinpath(data_dir, "geo.xyz")

path_out_file_1d_x = joinpath(data_dir, "rhine_data_1d_20_x_geo.txt")
path_out_file_1d_y = joinpath(data_dir, "rhine_data_1d_20_y_geo.txt")
path_out_file_2d = joinpath(data_dir, "rhine_data_2d_20_geo.txt")

# Get grid dimensions from the CartData object
nx = size(Topo_Cart_orth.x.val, 1)
ny = size(Topo_Cart_orth.y.val, 2)

# Convert to 1D format in x-direction (taking section 10 in y-direction)
convert_geo_1d(path_src_file, path_out_file_1d_x, nx=nx, ny=ny; 
               excerpt=20, section=10)

# Convert to 1D format in y-direction (taking section 100 in x-direction)
convert_geo_1d(path_src_file, path_out_file_1d_y, nx=nx, ny=ny; 
               excerpt=20, direction="y", section=100)
```

The `convert_geo_1d` function parameters:
- **`nx`, `ny`**: Grid dimensions from the original structured data.
- **`excerpt`**: Stride parameter - extracts every nth data point (20 means every 20th point).
- **`direction`**: Either "x" or "y" to specify the extraction direction.
- **`section`**: Which cross-section to extract from the perpendicular direction.

### Converting to 2D format

```@example create_convert_geo
# Convert to 2D format
convert_geo_2d(path_src_file, path_out_file_2d, nx=nx, ny=ny; excerpt=20)
```

The `convert_geo_2d` function preserves the two-dimensional structure while applying the stride parameter to reduce data density.
