using TrixiBottomTopography
using Downloads: download

# Define file paths
path_src_file = download("https://gist.githubusercontent.com/maxbertrand1996/c6917dcf80aef1704c633ec643a531d5/raw/f09b43f604adf9e2cfb45a7d998418f1e72f251d/dgm1_32_357_5646_1_nw.xyz")
path_out_file_1d_1_x_1 = "data/rhine_1d_1_x_1.txt"
path_out_file_1d_1_y_1 = "data/rhine_1d_1_y_1.txt"
path_out_file_1d_1_x_100 = "data/rhine_1d_1_x_100.txt"
path_out_file_1d_1_y_100 = "data/rhine_1d_1_y_100.txt"
path_out_file_1d_10_x_1 = "data/rhine_1d_10_x_1.txt"
path_out_file_1d_10_y_1 = "data/rhine_1d_10_y_1.txt"
path_out_file_1d_10_x_100 = "data/rhine_1d_10_x_100.txt"
path_out_file_1d_10_y_100 = "data/rhine_1d_10_y_100.txt"
path_out_file_2d = "data/rhine_2d_1.txt"
path_out_file_2d_10 = "data/rhine_2d_10.txt"

# Convert data
convert_dgm_1d(path_src_file, path_out_file_1d_1_x_1)
convert_dgm_1d(path_src_file, path_out_file_1d_1_y_1; direction = "y")
convert_dgm_1d(path_src_file, path_out_file_1d_1_x_100; section = 100)
convert_dgm_1d(path_src_file, path_out_file_1d_1_y_100; direction = "y", section = 100)
convert_dgm_1d(path_src_file, path_out_file_1d_10_x_1; excerpt = 10)
convert_dgm_1d(path_src_file, path_out_file_1d_10_y_1; excerpt = 10, direction = "y")
convert_dgm_1d(path_src_file, path_out_file_1d_10_x_100; excerpt = 10, section = 100)
convert_dgm_1d(path_src_file, path_out_file_1d_10_y_100; excerpt = 10, direction = "y", section = 100)
convert_dgm_2d(path_src_file, path_out_file_2d)
convert_dgm_2d(path_src_file, path_out_file_2d_10; excerpt = 10)

# Check if files exist
@test isfile(path_out_file_1d_1_x_1)
@test isfile(path_out_file_1d_1_y_1)
@test isfile(path_out_file_1d_1_x_100)
@test isfile(path_out_file_1d_1_y_100)
@test isfile(path_out_file_1d_10_x_1)
@test isfile(path_out_file_1d_10_y_1)
@test isfile(path_out_file_1d_10_x_100)
@test isfile(path_out_file_1d_10_y_100)
@test isfile(path_out_file_2d)
@test isfile(path_out_file_2d_10)