###############################################################################
# This example file creates a .gif file of a one dimensional wave running up  #
# the Cliffs of Moher, Ireland. The bottom topography is obtained with        #
# GeophysicalModelGenerator.jl, see `examples/create_convert_data_geo.jl`.    #
###############################################################################

# Include packages
using TrixiBottomTopography
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater

# Load the one dimensional cut through the Cliffs of Moher topography.
root_dir = pkgdir(TrixiBottomTopography)
cliffs_data = joinpath(root_dir, "examples", "data", "cliffs_data_1d_10_x.txt")

# B-spline interpolation of the underlying data. No smoothing is applied here because it
# would flatten the steep cliff.
spline_struct = CubicBSpline(cliffs_data; end_condition = "not-a-knot")
spline_func(x) = spline_interpolation(spline_struct, x)

# Defining one dimensional shallow water equations. The topography is given with respect to
# sea level, so a positive `H0` floods the shallow shelf in front of the cliffs while the
# cliff face itself stays dry.
equations = ShallowWaterEquations1D(gravity = 9.81, H0 = 30.0)

# Amplitude, crest position and width of the initial wave. The crest is placed on the flat
# shelf in front of the cliffs where the undisturbed water depth is constant.
wave_amplitude = 15.0
wave_center = -100.0
wave_width = 120.0

# Defining initial condition of a wave which travels towards the cliffs
function initial_condition_wave(x, t, equations::ShallowWaterEquations1D)
    b = spline_func(x[1])

    # Undisturbed water depth
    # Note, it is mandatory to shift the water level at dry areas to ensure that the water height
    # stays positive. The system would not be stable for `h` set to a hard `0` due to division in
    # the computation of velocity, e.g., (h v) / h. The small dry state threshold has
    # a default value of 500*eps() ≈ 1e-13 in double precision.
    h0 = max(equations.H0 - b, equations.threshold_limiter)

    # Smooth elevation of the water surface on top of it
    h = h0 + wave_amplitude * exp(-((x[1] - wave_center) / wave_width)^2)

    # Velocity of a right running simple wave
    v = 2 * (sqrt(equations.gravity * h) - sqrt(equations.gravity * h0))

    return SVector(h, h * v, b)
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

coordinates_min = spline_struct.x[1]
coordinates_max = spline_struct.x[end]
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                periodicity = false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers

tspan = (0.0, 45.0)
ode = semidiscretize(semi, tspan)

# positivity limiter for the water height
stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

###############################################################################
# run the simulation

# equidistant nodes in time for visualization
visnodes = range(tspan[1], tspan[2], length = 90)

sol = solve(ode, SSPRK43(; stage_limiter!), abstol = 1.0e-6, reltol = 1.0e-6,
            saveat = visnodes);

# Create an animation of the solution
if isdefined(Main, :Makie)
    j = Makie.Observable(1)
    time = Makie.Observable(0.0)

    pd_list = [PlotData1D(sol.u[i], semi) for i in 1:length(sol.t)]
    f = Makie.Figure()
    title_string = Makie.lift(t -> "time t = $(round(t, digits=3))", time)
    ax = Makie.Axis(f[1, 1], xlabel = "x [m]", ylabel = "z [m]", title = title_string)

    height = Makie.lift(i -> pd_list[i].data[:, 1], j)
    bottom = Makie.lift(i -> pd_list[i].data[:, 3], j)
    Makie.lines!(ax, pd_list[1].x, height)
    Makie.lines!(ax, pd_list[1].x, bottom)

    Makie.record(f, "animation_cliffs_1d.gif", 1:length(pd_list)) do tt
        j[] = tt
        time[] = sol.t[tt]
    end
end
