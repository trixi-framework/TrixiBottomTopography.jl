
# Examples with real topography data from GeophysicalModelGenerator.jl

This section demonstrates how to use the topography data created with GeophysicalModelGenerator.jl functions (`geo_topo_impression` and `create_topography_data`) in combination with [TrixiShallowWater.jl](https://github.com/trixi-framework/TrixiShallowWater.jl) to simulate shallow water flow over real-world topography.

The workflow follows these steps:
1. Load the converted topography data files (as shown in `docs/src/create_convert_geo_data.md`).
2. Create B-spline interpolations of the topography.
3. Set up and run shallow water simulations.
4. Visualize the results.

The example files used in this section can be found here:
-  1D example:  `examples/trixishallowwater_damn_break_1D_geo_data.jl`
- 2D example: `examples/trixishallowwater_damn_break_2D_geo_data.jl`

## One dimensional dam break with real topography

This example demonstrates a 1D dam break simulation using real Rhine River topography data processed through the GeophysicalModelGenerator workflow.

### Setup and data loading

First, we include the necessary packages and load the topography data:

```@example geo1d
# Include packages
using TrixiBottomTopography
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using CairoMakie
using GeophysicalModelGenerator
using GMT


data_dir = joinpath(@__DIR__, "examples/data")
mkpath(data_dir)

data_file = joinpath(data_dir, "rhine_data_1d_20_x_geo.txt")

data = data_file

```

### B-spline interpolation

The topography data is interpolated using cubic B-splines to create a smooth bottom function:

```@example geo1d
# Define B-spline structure
spline_struct = CubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)

spline_func(x) = spline_interpolation(spline_struct, x)
```

### Visualization of the topography

Before running the simulation, we can visualize the interpolated topography:

```@example geo1d
# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))

# Get interpolated values
y_int_pts = spline_func.(x_int_pts)

# Plot the topography
plot_topography(x_int_pts, y_int_pts; xlabel = "x[m]", ylabel = "z[m]")
```

### Shallow water equations setup

We define the 1D shallow water equations with appropriate physical parameters:

```@example geo1d
equations = ShallowWaterEquations1D(gravity = 1.0, H0 = 55.0)
```

### Initial condition for dam break

The initial condition creates a dam break scenario where water is initially higher in a central region:

```@example geo1d
# Defining initial condition for the dam break problem
function initial_condition_dam_break(x, t, equations::ShallowWaterEquations1D)
    inicenter = SVector(0.0)
    x_norm = x[1] - inicenter[1]
    r = abs(x_norm)

    # Calculate primitive variables
    H = r < 50 ? 70.0 : 60.0  # Higher water in center region
    v = 0.0                   # Initial velocity is zero
    b = spline_func(x[1])     # Bottom topography from B-spline

    return prim2cons(SVector(H, v, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_dam_break

# Setting the boundary to be a reflective wall
boundary_condition = boundary_condition_slip_wall
```

### Numerical solver setup

Get the DG approximation in space:

```@example geo1d
###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(polydeg = 3, surface_flux = (flux_hll, flux_nonconservative_fjordholm_etal),
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
```

### Mesh and semidiscretization

Here we use a TreeMesh:

```@example geo1d
###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = spline_struct.x[1]
coordinates_max = spline_struct.x[end]
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
                n_cells_max = 10_000,
                periodicity = false)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)
```

### Time integration and solution

Finally solving the PDE:

```@example geo1d
###############################################################################
# ODE solvers

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# run the simulation

# define equidistant nodes in time for visualization of an animation
visnodes = range(tspan[1], tspan[2], length = 90)

# use a Runge-Kutta method with error-based time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8,
            saveat = visnodes)
```

### Animation creation

The solution is visualized as an animation showing the evolution of water height over the real topography:

```@example geo1d
# Create animation of the solution
j = Observable(1)
time = Observable(0.0)

pd_list = [PlotData1D(sol.u[i], semi) for i in 1:length(sol.t)]
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
```

## Two dimensional dam break with real topography

This example extends the simulation to 2D, using the real Rhine River topography data processed through GeophysicalModelGenerator.

### Setup and data loading

```@example geo2d
# Include packages
using TrixiBottomTopography
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using CairoMakie
using Trixi2Vtk
using GeophysicalModelGenerator
using GMT


data_dir = joinpath(@__DIR__, "examples/data")
mkpath(data_dir)

data_file = joinpath(data_dir, "rhine_data_2d_20_geo.txt")

data = data_file

```

### Bicubic B-spline interpolation

For 2D topography, we use bicubic B-splines:

```@example geo2d
# B-spline interpolation of the underlying data
spline_struct = BicubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)

# Define B-spline interpolation function
spline_func(x, y) = spline_interpolation(spline_struct, x, y)
```

### Visualization of 2D topography

```@example geo2d
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
```

### 2D shallow water equations

```@example geo2d
equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 70.0)

function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)
    inicenter = SVector(-10.0, -20.0) # center of the domain
    x_norm = x - inicenter
    r = sqrt(x_norm[1]^2 + x_norm[2]^2)

    # Calculate primitive variables
    H = r < 50 ? 80.0 : 70.0 # Higher water in center region
    v1 = 0.0                 # Initial x-velocity
    v2 = 0.0                 # Initial y-velocity

    x1, x2 = x
    b = spline_func(x1, x2)  # Bottom topography from bicubic B-spline

    return prim2cons(SVector(H, v1, v2, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_wave

# Setting boundary conditions for all sides
boundary_condition = Dict(:x_neg => boundary_condition_slip_wall,
                          :y_neg => boundary_condition_slip_wall,
                          :y_pos => boundary_condition_slip_wall,
                          :x_pos => boundary_condition_slip_wall)
```

### 2D solver and mesh setup

```@example geo2d
###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(polydeg = 3,
               surface_flux = (flux_fjordholm_etal, flux_nonconservative_fjordholm_etal),
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# create the mesh and semidiscretization
coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])

mesh = P4estMesh((1, 1);
                 polydeg = 1,
                 coordinates_min = coordinates_min,
                 coordinates_max = coordinates_max,
                 initial_refinement_level = 6,
                 periodicity = false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)
```
**Note on mesh choice**: For this 2D simulation, we use `P4estMesh` instead of a `TreeMesh` (as used in the 1D case) because the `TreeMesh` requires cubic domains which we do not have here.

### Solution with callbacks and output

For the 2D case, we set up callbacks to save the solution for post-processing with ParaView:

```@example geo2d
tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Clear the output directory if it exists and create it new for saving the output

output_dir = "examples/data/out"
if isdir(output_dir)
    rm(output_dir, recursive = true)
end
mkpath(output_dir)

###############################################################################
# Callbacks for the ODE solver

stepsize_callback = StepsizeCallback(cfl = 0.6)

save_solution = SaveSolutionCallback(interval = 2,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = output_dir,
                                     solution_variables = cons2prim)

callbacks = CallbackSet(stepsize_callback, save_solution)

###############################################################################
# Solve with positivity preserving limiter

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

sol = solve(ode, RDPK3SpFSAL49(stage_limiter!), dt = 1.0, adaptive = false,
            callback = callbacks)
```

### Post-processing for ParaView

The 2D results are converted to VTK format for visualization in ParaView:

```@example geo2d
# Save mesh and convert to VTK format
Trixi.save_mesh_file(mesh, output_dir)
trixi2vtk(joinpath(output_dir, "solution_*.h5"), output_directory = output_dir)
```

## Visualization in ParaView

For 2D simulations, the generated VTK files can be opened in ParaView for advanced visualization:

1. Open the `solution_00000.pvd` file in ParaView.
2. Apply "Warp By Scalar" filter twice: once for water height, once for bathymetry.
3. Customize colors, lighting, and camera angles.
4. Create animations of the time evolution.

Example videos of such simulations are available here:
- [3D perspective view](https://jgumainz-my.sharepoint.com/:v:/g/personal/vimarks_uni-mainz_de/EVy3I6lXCXZFkbtQpdEbFX4BL-3bJ4-ueNPHsqYd__0pdA)
- [Top-down view](https://jgumainz-my.sharepoint.com/:v:/g/personal/vimarks_uni-mainz_de/EQ9GSqVOilJKjmWUElKHq5UBylGwgLlPLY1oJsRx3r4s6Q)