#the goal is to use the GeophysicalModelGenerator.jl to load real topography data 

# wie lade ich beliebeige topographische daten in GeophysicalModelGenerator.jl:
# es gibt import topo, gmg laden und gmt, dann gibt es eine funktion import_topo, denn gebe ich die longitude und lattitude ein und dann läd man das runter
# 
import Pkg
Pkg.add(["GeophysicalModelGenerator", "GMT", "Plots"])
using GMT
using Plots
using GeophysicalModelGenerator


##########################
#some topography data from the Rhine close to the Theodor-Heuss-Brücke in Mainz
#50.008306, 8.276694
#50.008200, 8.277027
#50.007888, 8.276223
#50.008217, 8.275613
###########################
# Specify the limits of the topography data based on the coordinates
lon_min = 8.276223
lon_max = 8.277027
lat_min = 50.007888
lat_max = 50.008306
limits=[lon_min, lon_max, lat_min, lat_max]

lon_mean = (lon_max+lon_min)/2
lat_mean=(lat_min+lat_max)/2
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

# *Note*: this routine is only available once the GMT.jl package is loaded in the REPL

Topo = import_topo(limits; file="@earth_relief_01s", maxattempts=5) 


#Topo = import_topo(lon = [-10, 45], lat=[25, 50], file="@earth_relief_20m")
p=ProjectionPoint(Lon=lon_mean, Lat=lat_mean)
Topo_Cart = convert2CartData(Topo,p)#Topo_Cart = convert2UTMzone(Topo,p) #convert2UTMzone(Topo,p),convert2CartData(Topo,p)
write_paraview(Topo_Cart,"Rhein_highest");

############################
#building a .vts file to first get an impression of the topography data in Paraview


    
#Uses `GMT` to download the topography of a certain region, specified with limits=[lon_min, lon_max, lat_min, lat_max].
#Sometimes download fails because of the internet connection. We do `maxattempts` to download it.

# using GeophysicalModelGenerator, GMT
# Topo = import_topo([4,20,37,49], file="@earth_relief_01m")

Topo = import_topo(lon = [-10, 45], lat=[25, 50], file="@earth_relief_20m")

convert(UTMData, Topo)

p=ProjectionPoint(Lon=17.3, Lat=37.5)

convert2UTMzone(Topo,p)

Topo_Cart = convert2CartData(Topo,p)

Topo_Cart_orth  = CartData(xyz_grid(-2000:20:2000,-1000:20:1000,0))

Topo_Cart_orth  = project_CartData(Topo_Cart_orth, Topo, p)