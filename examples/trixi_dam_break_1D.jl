###############################################################################
# This example file creates a .gif file of a one dimensional dam break        #
# where the bottom topography is chosen from the DGM data set.                #
###############################################################################

# Include packages
using TrixiBottomTopography
using Plots
using OrdinaryDiffEqLowStorageRK
using Trixi

# Download one dimensional Rhine bottom data from gist
Rhine_data = Trixi.download("https://gist.githubusercontent.com/maxbertrand1996/19c33682b99bfb1cc3116f31dd49bdb9/raw/d96499a1ffe250bc8e4cca8622779bae61543fd8/Rhine_data_1D_40_x_841.txt")

# B-spline interpolation of the underlying data
spline_struct = CubicBSpline(Rhine_data)
spline_func(x) = spline_interpolation(spline_struct, x)

# Defining one dimensional shallow water equations
equations = ShallowWaterEquations1D(gravity_constant=1.0, H0=55.0)

# Defining initial condition for the dam break problem
function initial_condition_dam_break(x, t, equations::ShallowWaterEquations1D)

  inicenter = SVector(357490.0)
  x_norm = x[1] - inicenter[1]
  r = abs(x_norm)

  # Calculate primitive variables
  H =  r < 50 ? 60.0 : 55.0
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
solver = DGSEM(polydeg=3, surface_flux=(flux_hll, flux_nonconservative_fjordholm_etal),
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh

coordinates_min = spline_struct.x[1]
coordinates_max = spline_struct.x[end]
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=3,
                n_cells_max=10_000,
                periodicity = false)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# run the simulation

# use a Runge-Kutta method with error based time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol=1.0e-8, reltol=1.0e-8,
            save_everystep=true);

# Create .gif animation of the solution
pyplot()
animation = @animate for k= 1:2:length(sol.t)
  pd = PlotData1D(sol.u[k], semi)
  plot(pd["H"])
  plot!(pd["b"], ylim=(38,65), title="t=$(sol.t[k])", xlabel="ETRS89 East", ylabel="DHHN2016")
end

gif(animation, "examples\\plots\\dam_break_1d.gif", fps=15)