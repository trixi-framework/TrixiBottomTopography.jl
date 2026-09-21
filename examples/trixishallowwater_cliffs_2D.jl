###############################################################################
# This example file simulates a two dimensional wave breaking against the     #
# Cliffs of Moher, Ireland. The bottom topography is obtained with            #
# GeophysicalModelGenerator.jl, see `examples/create_convert_data_geo.jl`.    #
###############################################################################

# Include packages
using TrixiBottomTopography
using OrdinaryDiffEqSSPRK
using Trixi
using Trixi2Vtk
using TrixiShallowWater

# Load the two dimensional Cliffs of Moher topography. The domain covers the open sea
# and the cliffs.
root_dir = pkgdir(TrixiBottomTopography)
cliffs_data = joinpath(root_dir, "examples", "data", "cliffs_data_2d_10.txt")

# B-spline interpolation of the underlying data. No smoothing is applied here because it
# would flatten the steep cliff.
spline_struct = BicubicBSpline(cliffs_data; end_condition = "not-a-knot")
spline_func(x, y) = spline_interpolation(spline_struct, x, y)

###############################################################################
# Visualization of the topography

if isdefined(Main, :Makie)
    # Define interpolation points
    n = 200
    x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
    y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))

    # Get interpolated matrix
    z_int_pts = evaluate_two_dimensional_interpolant(spline_func, x_int_pts, y_int_pts)

    # Plot the topography
    plot_topography(x_int_pts,
                    y_int_pts,
                    permutedims(z_int_pts);
                    xlabel = "x\n [m]",
                    ylabel = "y\n [m]",
                    zlabel = "z\n [m]",
                    azimuth_angle = -120 * pi / 180,
                    elevation_angle = 20 * pi / 180)
end

###############################################################################
# Defining two dimensional shallow water equations

# The topography is given with respect to sea level. A positive background water height
# H_0 = 30 therefore floods the shallow shelf in front of the cliffs while the cliff face
# itself stays dry.
equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 30.0)

# Amplitude, crest position and width of the initial wave. The crest is placed on the flat
# shelf in front of the cliffs where the undisturbed water depth is constant.
wave_amplitude = 15.0
wave_center = -100.0
wave_width = 120.0

# Defining initial condition of a plane wave which travels towards the cliffs
function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)
    b = spline_func(x[1], x[2])

    # Undisturbed water depth
    h0 = max(equations.H0 - b, equations.threshold_limiter)

    # Smooth elevation of the water surface on top of it
    h = h0 + wave_amplitude * exp(-((x[1] - wave_center) / wave_width)^2)

    # Velocity of a right running simple wave
    v1 = 2 * (sqrt(equations.gravity * h) - sqrt(equations.gravity * h0))
    v2 = 0.0

    return SVector(h, h * v1, h * v2, b)
end

# Setting initial condition
initial_condition = initial_condition_wave

# Setting the boundary to be a free-slip wall
boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

basis = LobattoLegendreBasis(7)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = waterheight_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the mesh and semidiscretization

coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])
mesh = P4estMesh((1, 1);
                 polydeg = 1,
                 coordinates_min = coordinates_min,
                 coordinates_max = coordinates_max,
                 initial_refinement_level = 3,
                 periodicity = false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 45.0)
ode = semidiscretize(semi, tspan)

# Clear the output directory if it exists and create it anew for saving the output later
output_dir = "out"
if isdir(output_dir)
    rm(output_dir, recursive = true)
end
mkpath(output_dir)

# positivity limiter for the water height
stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

# adaptive mesh refinement to resolve the moving wave front
amr_indicator = IndicatorLöhner(semi, variable = first)
amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level = 3,
                                      med_level = 4, med_threshold = 0.1,
                                      max_level = 5, max_threshold = 0.5)
amr_callback = AMRCallback(semi, amr_controller,
                           interval = 1,
                           adapt_initial_condition = true,
                           adapt_initial_condition_only_refine = true,
                           limiter! = stage_limiter!)

stepsize_callback = StepsizeCallback(cfl = 0.2)

save_solution = SaveSolutionCallback(dt = 0.5,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = output_dir,
                                     solution_variables = cons2prim)

callbacks = CallbackSet(amr_callback, stepsize_callback, save_solution)

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(; stage_limiter!); dt = 1.0, adaptive = false,
            callback = callbacks)

# To visualize the solution and bathymetry we post-process the Trixi.jl output file(s)
# with the Trixi2Vtk.jl functionality and plot them with ParaView.
trixi2vtk(joinpath(output_dir, "solution_*.h5"), output_directory = output_dir)

# It is possible to open the created solution_00000.pvd file with ParaView and create a
# video of the simulation.
#
# In ParaView, after opening the solution_00000.pvd file, one can apply two instances
# of the Warp By Scalar filter to visualize the water height and bathymetry.
#
# For example: https://jgumainz-my.sharepoint.com/:v:/g/personal/vimarks_uni-mainz_de/IQAKqXmSRmVDSavwyq4GUSzwAdsBEl7-8c_S_O8agqehVZk?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=cfnoyx
