###############################################################################
# This example file creates a .gif file of a two dimensional dam break        #
# where the bottom topography is chosen from the DGM data set.                #
###############################################################################

# Include packages
using TrixiBottomTopography
# using Plots
using LinearAlgebra
using OrdinaryDiffEq
using Trixi

Rhine_data = download("https://gist.githubusercontent.com/maxbertrand1996/a30db4dc9f5427c78160321d75a08166/raw/fa53ceb39ac82a6966cbb14e1220656cf7f97c1b/Rhine_data_2D_40.txt")

spline_struct = BicubicBSpline(Rhine_data)
spline_func(x,y) = spline_interpolation(spline_struct, x, y)

equations = ShallowWaterEquations2D(gravity_constant=3.0, H0=55.0)

function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)

  inicenter = SVector(357490.0, 5646519.0)
  x_norm = x - inicenter
  r = norm(x_norm)

  # Calculate primitive variables
  H =  r < 50 ? 65.0 : 55.0
  v1 = 0.0
  v2 = 0.0

  x1, x2 = x
  b = spline_func(x1, x2)

  return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_wave

# boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(polydeg=3, surface_flux=(flux_fjordholm_etal, flux_nonconservative_fjordholm_etal),
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=3,
                n_cells_max=10_000)#,
                # periodicity = false)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)#,
                                    # boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# run the simulation

# use a Runge-Kutta method with automatic (error based) time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol=1.0e-8, reltol=1.0e-8, save_everystep=true);

# display(length(sol.t))

# Create .gif animation of the solution
pyplot()
animation = @animate for k= 1:6:length(sol.t)
  pd = PlotData2D(sol.u[k], semi)
  wireframe(pd["H"])
  surface!(pd["b"], zlim=(38,65), camera = (30,20), title="t=$(sol.t[k])",
            xlabel="E", ylabel="N", zlabel="H")
end

gif(animation, "examples\\plots\\dam_break_2d.gif", fps=15)