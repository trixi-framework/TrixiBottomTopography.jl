# One dimensional dam brack

As mentioned in the [Home](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/)
section of this documentation, `TrixiBottomTopography.jl` was initially developed as a
supplementary package for the numerical solver [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)
to enable the user to use real life geographical data for the bottom topography
function of the shallow water equations.

In this section a one dimensional example is presented which uses the functionalities of
`TrixiBottomTopography.jl` with [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)
to simulate a dam break problem.

The underlying example file can be found [here](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/blob/main/examples/trixi_dam_break_1D.jl).

First, all the necessary packages must be included at the beginning of the file.
```julia
# Include packages
using TrixiBottomTopography
using Plots
using OrdinaryDiffEq
using Trixi
```
- `Plots` is responsible for visualizing the approximate solution of the dam break problem.
- `OrdinaryDiffEq` always has to be added when working with `Trixi`.

Next, the underlying bottom topography data is downloaded from a gist.
```julia
# Download one dimensional Rhine bottom data from gist
Rhine_data = download("https://gist.githubusercontent.com/maxbertrand1996/19c33682b99bfb1cc3116f31dd49bdb9/raw/d96499a1ffe250bc8e4cca8622779bae61543fd8/Rhine_data_1D_40_x_841.txt")
```
The downloaded data is then used to define the B-spline interpolation function as described in
[B-spline structure]("https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/structure/")
and [B-spline function]("https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/function/").
In this case a cubic B-spline interpolation function with free end condition is chosen.
```julia
# B-spline interpolation of the underlying data
spline_struct = CubicBSpline(Rhine_data)
spline_func(x) = spline_interpolation(spline_struct, x)
```
Now that the B-spline interpolation function is determined, the one dimensional shallow water equations implemented in `Trixi.jl` can be defined by calling:
```julia
# Defining one dimensional shallow water equations
equations = ShallowWaterEquations1D(gravity_constant=1.0, H0=55.0)
```
Here the gravity constant has been chosen to be $1.0$ and the background
total water height $H_0$ has been set to $55.0$.

Next the initial condition for the dam break problem can be defined.
At time $t=0$, a part of the water hight in the center of the domain with a diameter of $100$
is set to $60.0$ while the rest of the domain stays at the background water height $55.0$.
Additionally we can see that the bottom topography `b` is defined by the
B-spline interpolation function `spline_func` and is set in the initial condition.
```julia
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
```
After the initial condition we can set the boundary conditions.
In this case a reflective wall condition is chosen, which is already implemented
in `Trixi.jl` for the one dimensional shallow water equations.
```julia
# Setting initaial condition
initial_condition = initial_condition_dam_break

# Setting the boundary to be a reflective wall
boundary_condition = boundary_condition_slip_wall
```
The upcoming code parts will not be covered in full detail.
To get a more profound understanding of the routines, please see the
[Trixi.jl documentation](https://trixi-framework.github.io/Trixi.jl/stable/).

The following code snippet sets up the discontinous Galerking spectral element method (DGSEM).
In this solver-type, we can specify which flux functions for the surface and volume fluxes
will be taken, as well as the polynomial degree (`polydeg`) of the polynomials used
in the approximation space.
```julia
###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(polydeg=3, surface_flux=(flux_hll, flux_nonconservative_fjordholm_etal),
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))
```
After the solver comes the specification of the mesh in the approximation.
In this case a [`TreeMesh`](https://trixi-framework.github.io/Trixi.jl/stable/meshes/tree_mesh/) is chosen, which is a Cartesian mesh.
Here the domain borders must be defined, as well as the number of initial elements
($2$ to the power of `inital_refinement_level`).
Also, we have to indicate if the domain is periodic.
In this example boundary conditions were defined, thus the periodicity is set to `false`.

Once the underlying mesh is constructed, a semidiscretization object can be created
by calling `SemiDiscretizationHyperbolic`. This collects all the building blocks needed to set up the semi discretization:
- The underlying mesh.
- The set of equations.
- The initial condition.
- The solver (in this case DGSEM).
- The boundary conditions.
```julia
###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = spline_struct.x[1]
coordinates_max = spline_struct.x[end]
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=3,
                n_cells_max=10_000,
                periodicity = false)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)
```
Now an ordinary differential equations object is set up using a specified time range `tspan` and the semidiscretization object `semi`.
```julia
###############################################################################
# ODE solvers

tspan = (0.0, 0.0)
ode = semidiscretize(semi, tspan)
```
The ordinary differential equations object `ode` is solved by the function `sol`
which is part of the `OrdinaryDiffEq` package. Here the time stepping method can
be specified (in this case `RDPK3SpFSAL49()`) as well as some tolerances which
are responsible for an error based time step control.
```julia
###############################################################################
# run the simulation

# use a Runge-Kutta method with automatic (error based) time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol=1.0e-8, reltol=1.0e-8,
            save_everystep=false);
```
At this point the calculations would normally be finished.
But to have a nice visualization of the dam break problem, we want to create a .gif
file of the solution. To do so, a rather unorthodox approach is chosen and outlined below.

Creating a `PlotData1D` object of the solution enables the user to plot the solution
in a nice interpolated way. The only drawback is that we can only plot the final solution and
not the interim steps as the solution evolves. To get around this, a for loop is started which
increases `tspan` in every time step and performs the ODE calculation over and over again.
The solutions are saved in the vector `sol_vec`.

Although this seems a bit overkill, the calculations in one spatial dimension are actually
so fast that this does not have a huge impact on the overall calculation time.
In fact, the longest part is taken up by creating the .gif file from the solutions vector.

```julia
pd = PlotData1D(sol)

sol_vec = [pd]

# Run for t = 1,...,100
for i = 1:100

  ###############################################################################
  # ODE solvers

  local tspan = (0.0, convert(Float64, i))
  local ode = semidiscretize(semi, tspan)

  ###############################################################################
  # run the simulation

  # use a Runge-Kutta method with automatic (error based) time step size control
  local sol = solve(ode, RDPK3SpFSAL49(), abstol=1.0e-8, reltol=1.0e-8,
              save_everystep=false);

  local pd = PlotData1D(sol)

  append!(sol_vec, [pd])
end
```
From the variable `sol_vec`, the .gif file can be created using the macro `@animate`.
```julia
# Create .gif animation of the solution
pyplot()
animation = @animate for k= 1:101
    plot(sol_vec[k]["H"])
    plot!(sol_vec[k]["b"], ylim=(38,65), title="t=$(k-1)", xlabel="ETRS89 East", ylabel="DHHN2016")
end

gif(animation, "examples\\plots\\dam_break_1d.gif", fps=10)
```
This is the resulting .gif.

![gif](https://user-images.githubusercontent.com/101979498/203507054-2faca609-2628-4fea-9a4c-5788d02a237b.gif)