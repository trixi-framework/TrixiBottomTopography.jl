##########################################################################
# This script converts the data from                                     #
# https://www.opengeodata.nrw.de/produkte/geobasis/hm/dgm1_xyz/dgm1_xyz/ #
# into data files which can be used by TrixiBottomTopography             #
##########################################################################

# Include packages
using TrixiBottomTopography
using Downloads: download

# Get root directory
dir_path = pkgdir(TrixiBottomTopography)

# Define file paths
path_src_file = download("https://gist.githubusercontent.com/maxbertrand1996/c6917dcf80aef1704c633ec643a531d5/raw/f09b43f604adf9e2cfb45a7d998418f1e72f251d/dgm1_32_357_5646_1_nw.xyz")
path_out_file_1d_x = string(dir_path, "/examples/data/rhine_data_1d_20_x.txt")
path_out_file_1d_y = string(dir_path, "/examples/data/rhine_data_1d_20_y.txt")
path_out_file_2d = string(dir_path, "/examples/data/rhine_data_2d_20.txt")

# Convert data
convert_dgm_1d(path_src_file, path_out_file_1d_x; excerpt = 20, section = 100)
convert_dgm_1d(path_src_file, path_out_file_1d_y; excerpt = 20, direction = "y", section = 100)
convert_dgm_2d(path_src_file, path_out_file_2d; excerpt = 20)