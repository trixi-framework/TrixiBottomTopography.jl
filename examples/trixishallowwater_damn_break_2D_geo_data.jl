###################################################################################
# This example file visualises the topography data and creats a .gif file of a two #
#dimensional dam break where the bottom topography is is generated with            #
#GeophysicalModelGenerator.                                                        #
####################################################################################

# Include packages

using TrixiBottomTopography
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using CairoMakie

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
plot_topography(x_int_pts, y_int_pts, z_int_pts;
                xlabel = "x\n [m]",
                ylabel = "y\n [m]",
                zlabel = "z\n [m]",
                azimuth_angle = 54 * pi / 180,
                elevation_angle = 27 * pi / 180)

equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 65.0)

function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)
    inicenter = SVector(0.0, 0.0)
    x_norm = x - inicenter
    r = sqrt(x_norm[1]^2 + x_norm[2]^2)

    # Calculate primitive variables
    H = r < 50 ? 75.0 : 65.0
    v1 = 0.0
    v2 = 0.0

    x1, x2 = x
    b = spline_func(x1, x2)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_wave

# Setting the boundary to be a free-slip wall for a P4estMesh

boundary_condition = Dict(:x_neg => boundary_condition_slip_wall,
                          :y_neg => boundary_condition_slip_wall,
                          :y_pos => boundary_condition_slip_wall,
                          :x_pos => boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(polydeg = 3,
               surface_flux = (flux_fjordholm_etal, flux_nonconservative_fjordholm_etal),
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the P4estMesh and setup a periodic mesh. Here we use a P4estMesh instead of a TreeMesh
# because the mesh is not cubic.

coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])

mesh = P4estMesh((1, 1);
                 polydeg = 1,
                 coordinates_min = coordinates_min,
                 coordinates_max = coordinates_max,
                 initial_refinement_level = 3,
                 periodicity = false)

####
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 70.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# run the simulation

# define equidistant nodes in time for visualization of an animation
visnodes = range(tspan[1], tspan[2], length = 175)

# use a Runge-Kutta method with error based time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8,
            saveat = visnodes)

##################################################################################
# Create an animation of the solution

j = Observable(1)
time = Observable(0.0)

f = Figure()
title_string = lift(t -> "Dam break simulation - time t = $(round(t, digits=3)) s", time)
ax = Axis3(f[1, 1], xlabel = "x[m]", ylabel = "y[m]", zlabel = "z[m]",
           title = title_string, azimuth = 130 * pi / 180, elevation = 18 * pi / 180)

x_coords = vec(pd_list[1].x)
y_coords = vec(pd_list[1].y)

record(f, "animation_2d.gif", 1:length(pd_list); framerate = 10) do tt
    # actualisation of the plot data
    height_data = [pd_list[tt].data[i][1] for i in 1:length(pd_list[tt].data)]
    bottom_data = [pd_list[tt].data[i][4] for i in 1:length(pd_list[tt].data)]

    #get the data for the current time step
    empty!(ax)
    scatter!(ax, x_coords, y_coords, bottom_data,
             color = bottom_data, colormap = :terrain, markersize = 4)
    scatter!(ax, x_coords, y_coords, height_data,
             color = height_data, colormap = :blues, markersize = 3)

    time[] = sol.t[tt]
    zlims!(ax, 35, 70)
end
