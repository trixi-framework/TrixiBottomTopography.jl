#the goal is to use the GeophysicalModelGenerator.jl to load real topography data 

# wie lade ich beliebeige topographische daten in GeophysicalModelGenerator.jl:
# es gibt import topo, gmg laden und gmt, dann gibt es eine funktion import_topo, denn gebe ich die longitude und lattitude ein und dann läd man das runter
# 
import Pkg
Pkg.add(["GeophysicalModelGenerator", "GMT", "Plots","CSV","DataFrames","TrixiBottomTopography", "Downloads"])
using GMT
using Plots
using GeophysicalModelGenerator
using CSV
using DataFrames
using TrixiBottomTopography
using Downloads: download

##########################
# some topography data from the Rhine close to the Theodor-Heuss-Brücke in Mainz
#50.010910 8.275880
#50.009480 8.270871
#50.006184 8.274388
#50.007914 8.278718



###########################
# Specify the limits of the topography data based on the coordinates
lon_min = 8.270871
lon_max = 8.278718
lat_min = 50.006184 
lat_max = 50.010910
limits=[lon_min, lon_max, lat_min, lat_max]

lon_mean = (lon_max+lon_min)/2
lat_mean=(lat_min+lat_max)/2

#rausfinden wie man die richtigen topographischen daten bekommt, anscheinend ist die gewählte form: 
#the first column provides the corresponding ETRS89 East coordinates, the second column the 
#ETRS89 North coordinates and the third column the DHHN2016 height.

###
#Lon=17.3, Lat=37.5
#Topo = import_topo(limits, file="@earth_relief_20m")
Topo = import_topo(lon = [8.270871, 8.278718], lat=[50.006184, 50.010910], file="@earth_relief_01m")
p=ProjectionPoint(Lon=lon_mean, Lat=lat_mean)
Topo_Cart = convert2CartData(Topo,p)

#es muss eine datei erstellt werden, sodass die wurzel der länge der datei ein integer ist



Topo_Cart_orth  = CartData(xyz_grid(-0.5:0.1:0.5,-0.5:0.1:0.5,0))

Topo_Cart_orth  = project_CartData(Topo_Cart_orth, Topo, p)

df_x = DataFrame(Topo_Cart_orth.x.val[:,:,1], :auto)

df_y = DataFrame(Topo_Cart_orth.y.val[:,:,1], :auto);

df_z = DataFrame(Topo_Cart_orth.z.val[:,:,1], :auto);

# Kombiniere die DataFrames für x, y und z in einen DataFrame
df_xyz = DataFrame(
                   x = convert.(Float64,vec(Topo_Cart_orth.x.val[:,:,1])), 
                   y = convert.(Float64,vec(Topo_Cart_orth.y.val[:,:,1])), 
                   z = convert.(Float64,vec(Topo_Cart_orth.z.val[:,:,1])),
                   )


 #alohomora, spieglein spieglein an der wand bitte löse meine probleme wie von zauberhand

# Schreibe die Datei ohne Header und mit Leerzeichen als Delimiter
# Make sure the data directory exists
data_dir = joinpath(@__DIR__, "data")
mkpath(data_dir)  # Create the directory if it doesn't exist

# Write directly to the data directory
output_file = joinpath(data_dir, "test.xyz")
open(output_file, "w") do file
    for row in eachrow(df_xyz)
        # Round each value in the row to 2 decimal places
        rounded_row = [round(value, digits=2) for value in row]
        println(file, join(rounded_row, " "))
    end
end

data_dir = joinpath(@__DIR__, "data")
data_dir
# Download the raw bottom topography data
path_src_file = joinpath(@__DIR__, "test.xyz")

path_out_file_1d_x = joinpath(data_dir, "rhine_data_1d_x_theodor.txt")
path_out_file_1d_y = joinpath(data_dir, "rhine_data_1d_20_y_theodor.txt")
path_out_file_2d = joinpath(data_dir, "rhine_data_2d_20_theodor.txt")

# Convert data
convert_dgm_1d(path_src_file, path_out_file_1d_x; excerpt = 1, section = 1)
convert_dgm_1d(path_src_file, path_out_file_1d_y; excerpt = 1, direction = "y", section = 1)
convert_dgm_2d(path_src_file, path_out_file_2d; excerpt = 1)

Int(sqrt(122))
round(Int, sqrt(14.45))
####################

#####################
# Define file paths
root_dir = pkgdir(TrixiBottomTopography)

# Download the raw bottom topography data
path_src_file = download("https://gist.githubusercontent.com/maxbertrand1996/c6917dcf80aef1704c633ec643a531d5/raw/f09b43f604adf9e2cfb45a7d998418f1e72f251d/dgm1_32_357_5646_1_nw.xyz",
                         joinpath(root_dir, "examples", "data",
                                  "dgm1_32_357_5646_1_nw.xyz"))

file_path = joinpath(root_dir, "examples", "data", "dgm1_32_357_5646_1_nw.xyz")
line_count = countlines(file_path)

path_out_file_1d_x = joinpath(root_dir, "examples", "data", "rhine_data_1d_20_x.txt")
path_out_file_1d_y = joinpath(root_dir, "examples", "rhine_data_1d_20_y.txt")
path_out_file_2d = joinpath(root_dir, "examples", "rhine_data_2d_20.txt")

# Convert data
convert_dgm_1d(path_src_file, path_out_file_1d_x; excerpt = 1, section = 1)
convert_dgm_1d(path_src_file, path_out_file_1d_y; excerpt = 20, direction = "y",
               section = 100)
convert_dgm_2d(path_src_file, path_out_file_2d; excerpt = 20)


print("lol")


############################
#Loading the topography data
#important: if you have a small area you might use hiher resolution: 
# Note: 
# ====
# - latitude values in the southern hemisphere should have a minus sign (e.g., -2.8)
# - longitude values that are "west" should *either* come with a minus sign *or* are defined by values >180

# | Dataset                 |   Resolution |   Description                                               |
# |:----------------        | ------------ | ----------------------------------------------------------- |
# | "@earth\\_relief\\_01s" |	1 arc sec 	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
# | "@earth\\_relief\\_03s"	|   3 arc sec	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
# | "@earth\\_relief\\_15s"	|  15 arc sec	 | SRTM15+ [David Sandwell, SIO/UCSD]                          |
# | "@earth\\_relief\\_30s"	|  30 arc sec	 | SRTM30+ [Becker et al., 2009, SIO/UCSD]                     |
# | "@earth\\_relief\\_01m"	|   1 arc min	 | ETOPO1 Ice surface [NEIC/NOAA]                              |
# | "@earth\\_relief\\_02m"	|   2 arc min	 | ETOPO2v2 Ice surface [NEIC/NOAA]                            |
# | "@earth\\_relief\\_03m"	|   3 arc min	 | ETOPO1 after Gaussian spherical filtering (5.6 km fullwidth)|
# | "@earth\\_relief\\_04m"	|   4 arc min	 | ETOPO1 after Gaussian spherical filtering (7.5 km fullwidth)|
# | "@earth\\_relief\\_05m"	|   5 arc min	 | ETOPO1 after Gaussian spherical filtering (9 km fullwidth)  |
# | "@earth\\_relief\\_06m"	|   6 arc min	 | ETOPO1 after Gaussia30n spherical filtering (10 km fullwidth) |
# | "@earth\\_relief\\_10m"	|  10 arc min	 | ETOPO1 after Gaussian spherical filtering (18 km fullwidth) |
# | "@earth\\_relief\\_15m"	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (28 km fullwidth) |
# | "@earth\\_relief\\_20m"	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (37 km fullwidth) |
# | "@earth\\_relief\\_30m"	|  30 arc min	 | ETOPO1 after Gaussian spherical filtering (55 km fullwidth) |
# | "@earth\\_relief\\_60m"	|  60 arc min	 | ETOPO1 after Gaussian spherical filtering (111 km fullwidth)|


# Topo = import_topo(lon = [-10, 45], lat=[25, 50], file="@earth_relief_20m")
# p=ProjectionPoint(Lon=17.3, Lat=37.5)
# convert2UTMzone(Topo,p)
# Topo_Cart = convert2CartData(Topo,p)
# Topo_Cart_orth  = CartData(xyz_grid(-2000:20:2000,-1000:20:1000,0))
# Topo_Cart_orth  = project_CartData(Topo_Cart_orth, Topo, p)
# write_paraview(Topo_Cart,"Topo_Cart_orth");


# #Topo = import_topo(lon = [-10, 45], lat=[25, 50], file="@earth_relief_20m")
# p=ProjectionPoint(Lon=lon_mean, Lat=lat_mean)
# Topo_Cart = convert2CartData(Topo,p)#Topo_Cart = convert2UTMzone(Topo,p) #convert2UTMzone(Topo,p),convert2CartData(Topo,p)
# write_paraview(Topo_Cart,"Rhein_highest");

# #load the csv file for the topography data from the paraview file
# # now we load the csv file for the topography data from the paraview file
# #in vs code you might get asked if you like to install Rainbow CSV
# # ...existing code...
# #df = CSV.read("/Users/Vincent_1/Library/CloudStorage/OneDrive-JGU(2)/HiWi/TrixiBottomTopography_geforkt/TrixiBottomTopography.jl/examples/rh.csv", DataFrame)
# df = CSV.read("/Users/Vincent_1/Library/CloudStorage/OneDrive-JGU(2)/HiWi/TrixiBottomTopography_geforkt/TrixiBottomTopography.jl/examples/modelgenerator3.csv", DataFrame)
# #try to load the data as an orthogonal cartesian grid, vlt geht das aber nicht 

# ############################
# #convert the csv file to a xyz file to be able to use it in TrixiBottomTopography

# function csv_to_xyz(csv_path, xyz_path)
#     # CSV-Datei einlesen
#     data = CSV.read(csv_path, DataFrame)
    
#     # Spaltenüberschriften bereinigen
#     column_names = names(data)
#     clean_names = [replace(String(name), "Points:" => "") for name in column_names]
#     rename!(data, Dict(zip(column_names, clean_names)))
    
#     # XYZ-Datei erstellen
#     open(xyz_path, "w") do file
#         for i in 1:size(data, 1)
#             x = data[i, "0"]
#             y = data[i, "1"]
#             z = data[i, "2"]
#             println(file, "$x $y $z")
#         end
#     end
    
#     println("Konvertierung abgeschlossen. XYZ-Datei wurde erstellt: $xyz_path")
# end

# ## Convert the CSV from the topo data to XYZ format:
# csv_to_xyz("/Users/Vincent_1/Library/CloudStorage/OneDrive-JGU(2)/HiWi/TrixiBottomTopography_geforkt/TrixiBottomTopography.jl/examples/modelgenerator3.csv", "/Users/Vincent_1/Library/CloudStorage/OneDrive-JGU(2)/HiWi/TrixiBottomTopography_geforkt/TrixiBottomTopography.jl/examples/modelgenerator3.xyz")

# # ############################################################################
# #convert this data to a format which can be used by TrixiBottomTopography, namely a special .txt file format

# # Define file paths
# data_dir = joinpath(@__DIR__, "data")
# data_dir
# # Download the raw bottom topography data
# path_src_file = joinpath(@__DIR__, "modelgenerator3.xyz")

# path_out_file_1d_x = joinpath(data_dir, "rhine_data_1d_x_theodor.txt")
# path_out_file_1d_y = joinpath(data_dir, "rhine_data_1d_20_y_theodor.txt")
# path_out_file_2d = joinpath(data_dir, "rhine_data_2d_20_theodor.txt")

# # Convert data
# convert_dgm_1d(path_src_file, path_out_file_1d_x; excerpt = 20, section = 10)
# convert_dgm_1d(path_src_file, path_out_file_1d_y; excerpt = 20, direction = "y", section = 10)
# convert_dgm_2d(path_src_file, path_out_file_2d; excerpt = 20)

# # Download one dimensional Rhine bottom data from gist
# #Rhine_data = download("https://gist.githubusercontent.com/maxbertrand1996/19c33682b99bfb1cc3116f31dd49bdb9/raw/d96499a1ffe250bc8e4cca8622779bae61543fd8/Rhine_data_1D_40_x_841.txt")

# #