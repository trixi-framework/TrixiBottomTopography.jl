using TrixiBottomTopography

# Define data path
data = "data/rhine_1d_10_x_1.txt"

# Define B-spline structure
spline_struct = linear_b_spline(data)
# Define B-spline interpolation function
spline_func(x) = spline_interpolation(spline_struct, x)

# Get y values from file
file = open(data)
lines = readlines(file)
close(file)

num_elements = parse(Int64,lines[2])
y = [parse(Float64, val) for val in lines[5+num_elements:end]]

# Check if x values put into spline function
# correspond to actual y values
@test spline_func.(spline_struct.x) == y