###################################################################################################
# # This example file visualises the topography data and creates a .gif file of a one dimensional # 
# dam break where the bottom topography is generated with GeophysicalModelGenerator.              #
###################################################################################################

using TrixiBottomTopography
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using CairoMakie

#get the data
data_dir = joinpath(@__DIR__, "data")
data = joinpath(data_dir, "rhine_data_1d_20_x_geo.txt")

# Define B-spline structure
spline_struct = CubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)

spline_func(x) = spline_interpolation(spline_struct, x)

# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))

# Get interpolated values
y_int_pts = spline_func.(x_int_pts)

# Plot the topography
plot_topography(x_int_pts, y_int_pts; xlabel = "x[m]", ylabel = "z[m]")

# Get the original interpolation knots
x_knots = spline_struct.x
y_knots = spline_func.(x_knots)

plot_topography_with_interpolation_knots(
    x_int_pts,
    y_int_pts,
    x_knots,
    y_knots;
    xlabel = "x[m]",
    ylabel = "z[m]",
)

equations = ShallowWaterEquations1D(gravity = 1.0, H0 = 55.0)
# Defining initial condition for the dam break problem
function initial_condition_dam_break(x, t, equations::ShallowWaterEquations1D)
    inicenter = SVector(0.0)
    x_norm = x[1] - inicenter[1]
    r = abs(x_norm)

    # Calculate primitive variables
    H = r < 50 ? 70.0 : 60.0
    v = 0.0
    b = spline_func(x[1])

    return prim2cons(SVector(H, v, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_dam_break

# Setting the boundary to be a reflective wall
boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(
    polydeg = 3,
    surface_flux = (flux_hll, flux_nonconservative_fjordholm_etal),
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
)

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = spline_struct.x[1]
coordinates_max = spline_struct.x[end]
mesh = TreeMesh(
    coordinates_min,
    coordinates_max,
    initial_refinement_level = 3,
    n_cells_max = 10_000,
    periodicity = false,
)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition,
    solver,
    boundary_conditions = boundary_condition,
)

###############################################################################
# ODE solvers

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# run the simulation

# define equidistant nodes in time for visualization of an animation
visnodes = range(tspan[1], tspan[2], length = 90)

# use a Runge-Kutta method with error-based time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8, saveat = visnodes)

# Create animation of the solution
j = Observable(1)
time = Observable(0.0)

pd_list = [PlotData1D(sol.u[i], semi) for i = 1:length(sol.t)]
f = Figure()
title_string = lift(t -> "time t = $(round(t, digits=3))", time)
ax = Axis(f[1, 1], xlabel = "x[m]", ylabel = "z[m]", title = title_string)

height = lift(i -> pd_list[i].data[:, 1], j)
bottom = lift(i -> pd_list[i].data[:, 3], j)
CairoMakie.lines!(ax, pd_list[1].x, height)
CairoMakie.lines!(ax, pd_list[1].x, bottom)
ylims!(ax, 30, 90)

record(f, "animation_geo.gif", 1:length(pd_list)) do tt
    j[] = tt
    time[] = sol.t[tt]
end
