###########################################################################
# This example file visualises the topography data and shows a simulation #
# created with Trixi2Vtk of a two dimensional dam break where the bottom  #
# topography is generated with GeophysicalModelGenerator.                 #
###########################################################################
# Include packages
using TrixiBottomTopography
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using CairoMakie: scatter
using Trixi2Vtk

# Load the created two dimensional Rhine bottom data topography data
data_dir = joinpath(@__DIR__, "data")
data = joinpath(data_dir, "rhine_data_2d_20_geo.txt")

# B-spline interpolation of the underlying data
spline_struct = BicubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)

# Define B-spline interpolation function
spline_func(x, y) = spline_interpolation(spline_struct, x, y)

# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))

# Get interpolated matrix

z_int_pts = evaluate_bicubicspline_interpolant(spline_func, x_int_pts, y_int_pts)

# Plot the topography
plot_topography(
    x_int_pts,
    y_int_pts,
    z_int_pts;
    xlabel = "x\n [m]",
    ylabel = "y\n [m]",
    zlabel = "z\n [m]",
    azimuth_angle = 54 * pi / 180,
    elevation_angle = 27 * pi / 180,
)

equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 70.0) # 55, 65

function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)
    inicenter = SVector(-10.0, -20.0) # center of the domain
    x_norm = x - inicenter
    r = sqrt(x_norm[1]^2 + x_norm[2]^2)

    # Calculate primitive variables
    H = r < 50 ? 80.0 : 70.0 # 65, 70
    v1 = 0.0
    v2 = 0.0

    x1, x2 = x
    b = spline_func(x1, x2)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_wave

# Setting the boundary to be a free-slip wall for a P4estMesh

boundary_condition = Dict(
    :x_neg => boundary_condition_slip_wall,
    :y_neg => boundary_condition_slip_wall,
    :y_pos => boundary_condition_slip_wall,
    :x_pos => boundary_condition_slip_wall,
)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(
    polydeg = 3,
    surface_flux = (flux_fjordholm_etal, flux_nonconservative_fjordholm_etal),
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
)

###############################################################################
# create the mesh and semidiscretization
coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])

mesh_1 = P4estMesh(
    (1, 1);
    polydeg = 1,
    coordinates_min = coordinates_min,
    coordinates_max = coordinates_max,
    initial_refinement_level = 6,
    periodicity = false,
)

semi = SemidiscretizationHyperbolic(
    mesh_1,
    equations,
    initial_condition,
    solver,
    boundary_conditions = boundary_condition,
)

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Clear the output directory if it exists and create it new for saving the output later

output_dir = "out"
if isdir(output_dir)
    rm(output_dir, recursive = true)
end
mkpath(output_dir)

###############################################################################
# Callbacks for the ODE solver

stepsize_callback = StepsizeCallback(cfl = 0.6)

save_solution = SaveSolutionCallback(
    interval = 2,
    save_initial_solution = true,
    save_final_solution = true,
    output_directory = output_dir,
    solution_variables = cons2prim,
)

callbacks = CallbackSet(stepsize_callback, save_solution);#amr_callback,

###############################################################################

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,)) #stage limiter for the water height

sol = solve(
    ode,
    RDPK3SpFSAL49(stage_limiter!),
    dt = 1.0,
    adaptive = false,
    callback = callbacks,
) # the stepsize_callback will adjust the time step size automatically

# To visualize the solution and bathymetry we post-processing the Trixi.jl output file(s)
# with the Trixi2Vtk.jl functionality and plot them with ParaView.

Trixi.save_mesh_file(mesh_1, output_dir)
trixi2vtk(joinpath(output_dir, "solution_*.h5"), output_directory = output_dir)

# It is possible to open the created solution_00000.pvd file with ParaView and create a video of the simulation.

# In ParaView, after opening the solution_00000.pvd file, one can apply two instances
# of the Warp By Scalar filter to visualize the water height and bathymetry in three dimensions.
# Many additional customizations, e.g., color scaling, fonts, etc. are available in ParaView.

# A created video of the simulation can be found here: https://jgumainz-my.sharepoint.com/:v:/g/personal/vimarks_uni-mainz_de/EVy3I6lXCXZFkbtQpdEbFX4BL-3bJ4-ueNPHsqYd__0pdA?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=whBBO4
# In addition there is also a top view of the solution available here: https://jgumainz-my.sharepoint.com/:v:/g/personal/vimarks_uni-mainz_de/EQ9GSqVOilJKjmWUElKHq5UBylGwgLlPLY1oJsRx3r4s6Q?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=pJIz3Z
