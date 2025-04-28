#the goal is to use the GeophysicalModelGenerator.jl to load real topography data 

# wie lade ich beliebeige topographische daten in GeophysicalModelGenerator.jl:
# es gibt import topo, gmg laden und gmt, dann gibt es eine funktion import_topo, denn gebe ich die longitude und lattitude ein und dann läd man das runter
# 
import Pkg
Pkg.add(["GeophysicalModelGenerator", "GMT", "Plots"])
using GMT
using Plots
using GeophysicalModelGenerator


lon_min, lon_max = -10.0, 10.0  # Longitude range from -10° to 10°
lat_min, lat_max = 45.0, 55.0   # Latitude range from 45° to 55°

# Load topography data using GMT.grdtrack
region = [lon_min, lon_max, lat_min, lat_max]
spacing = "0.1/0.1"  # Grid spacing in degrees
topo_data = GMT.grdtrack(region=region, spacing=spacing)

# Reshape the data for plotting
z = reshape(topo_data[:, 3], length(lon_min:0.1:lon_max), length(lat_min:0.1:lat_max))

# Visualize the data
heatmap(lon_min:0.1:lon_max, lat_min:0.1:lat_max, z', 
        xlabel="Longitude", 
        ylabel="Latitude",
        title="Topography")