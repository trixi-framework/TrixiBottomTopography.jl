# Examples with Trixi.jl

As mentioned in the [Home](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/)
section of this documentation, `TrixiBottomTopography.jl` was initially developed as a
supplementary package for the numerical solver [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)
to enable the user to use real world geographical data for the bottom topography
function of the shallow water equations.

In this section, a one dimensional example is presented which uses the functionalities of
`TrixiBottomTopography.jl` with [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)
to simulate a dam break problem.

The underlying example file can be found [here](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/trixi_dam_break_1D.jl).

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
[B-spline structure](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/structure/)
and [B-spline function](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/function/).
In this case, a cubic B-spline interpolation function with free end condition is chosen.
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
Here the gravity constant has been chosen to be $1.0$, and the background
total water height $H_0$ has been set to $55.0$.

Next, the initial condition for the dam break problem can be defined.
At time $t=0$, a part of the water height in the center of the domain with a diameter of $100$
is set to $60.0$ while the rest of the domain stays at the background water height $55.0$.
Additionally, we can see that the bottom topography `b` is defined by the
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
After the initial condition, we can set the boundary conditions.
In this case, a reflective wall condition is chosen, which is already implemented
in `Trixi.jl` for the one dimensional shallow water equations.
```julia
# Setting initaial condition
initial_condition = initial_condition_dam_break

# Setting the boundary to be a reflective wall
boundary_condition = boundary_condition_slip_wall
```
The upcoming code parts will **not** be covered in full detail.
To get a more profound understanding of the routines, please see the
[Trixi.jl documentation](https://trixi-framework.github.io/Trixi.jl/stable/).

The following code snippet sets up the discontinuous Galerkin spectral element method (DGSEM).
In this solver type, we can specify which flux functions for the surface and volume fluxes
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
In this case, a [`TreeMesh`](https://trixi-framework.github.io/Trixi.jl/stable/meshes/tree_mesh/) is chosen, which is a Cartesian mesh.
Here the domain borders must be defined, as well as the number of initial elements
($2$ to the power of `inital_refinement_level`).
Also, we have to indicate if the domain is periodic.
In this example, boundary conditions were defined. Thus the periodicity is set to `false`.

Once the underlying mesh is constructed, a semidiscretization object can be created
by calling `SemiDiscretizationHyperbolic`. This collects all the building blocks needed to set up the semidiscretization:
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
An ordinary differential equations object is set up using a specified time range, `tspan`, and the semidiscretization object, `semi`.
```julia
###############################################################################
# ODE solvers

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)
```
The ordinary differential equations object `ode` is solved by the function `sol`
which is part of the `OrdinaryDiffEq` package. Here the time stepping method can
be specified (in this case, `RDPK3SpFSAL49()`) as well as some tolerances responsible for an error-based time step control.
```julia
###############################################################################
# run the simulation

# use a Runge-Kutta method with automatic (error-based) time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol=1.0e-8, reltol=1.0e-8,
            save_everystep=true);
```
At this point, the calculations would generally be finished. However, to visualize the dam break problem, we want to create a .gif file of the solution. To do so, we have set the `save_everystep` attribute to `true`. This means that the solution for every time step will be callable afterwards.

First of all, a plotting backend is chosen. Here we use `pyplot()` as the resulting plots look very clear. Then we define an `animation` loop using the macro `@animate` over every second of the interim solutions. Inside the loop, the `PlotData2D` functionality from `Trixi.jl` is called to create a plotting object. Afterwards, this plotting object can be plotted using the known `plot` command.

The `gif` function uses `animation` to create a .gif from the plots for every second-time step and saves it in the specified location. Additionally, the frames per second rate can be set in the `fps` attribute.
```julia
# Create .gif animation of the solution
pyplot()
animation = @animate for k= 1:2:length(sol.t)
  pd = PlotData1D(sol.u[k], semi)
  plot(pd["H"])
  plot!(pd["b"], ylim=(38,65), title="t=$(sol.t[k])", xlabel="ETRS89 East", ylabel="DHHN2016")
end

gif(animation, "examples\\plots\\dam_break_1d.gif", fps=15)
```
This is the resulting .gif animation.

![gif](https://user-images.githubusercontent.com/101979498/203507054-2faca609-2628-4fea-9a4c-5788d02a237b.gif)

## Two dimensional dam break

The underlying example file can be found [here](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/trixi_dam_break_2D.jl).

The two dimensional example is very similar to the one dimensional case.

First, all the necessary packages and the underlying bottom topography data are loaded.

```julia
# Include packages
using TrixiBottomTopography
using Plots
using LinearAlgebra
using OrdinaryDiffEq
using Trixi

Rhine_data = download("https://gist.githubusercontent.com/maxbertrand1996/a30db4dc9f5427c78160321d75a08166/raw/fa53ceb39ac82a6966cbb14e1220656cf7f97c1b/Rhine_data_2D_40.txt")
```

Using the data, a bicubic B-spline interpolation is performed on the data to define a bottom topography function.

```julia
# B-spline interpolation of the underlying data
spline_struct = BicubicBSpline(Rhine_data)
spline_func(x,y) = spline_interpolation(spline_struct, x, y)
```

Then the two dimensional shallow water equations are defined, where the gravitational constant has been chosen to be `3.0` and the initial water height `55.0`. Afterwards, the initial condition is defined. Similar to the one dimensional case, in the center of the domain, a circular part with a diameter of `100.0` is chosen where the initial water height is chosen to be `10.0` units higher.

```julia
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
```

As we can see, there is no boundary condition specified. This is because, at this stage, `boundary_condition_slip_wall` has not been implemented into `Trixi.jl` for the two dimensional shallow water equations.

The DGSEM solver is set up as in the one dimensional case. 

```julia
###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(polydeg=3, surface_flux=(flux_fjordholm_etal, flux_nonconservative_fjordholm_etal),
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))
```

Now the mesh has to be specified. Because we do not have any boundary conditions defined, we can only assume to have a periodic domain. Therefore `periodicity` does not have to be specified as it is set to `true` by default.

```julia
###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=3,
                n_cells_max=10_000)
```

When calling the semidiscretization object again, `boundary_conditions` does not have to be specified.

```julia
# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)
```

The solution of the PDE and the.gif animation is analogous to the one dimensional case except that we chose `PlotData2D` to create the plotting object instead of `PlotData1D` as we are in the two dimensional case now.

```julia
###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# run the simulation

# use a Runge-Kutta method with automatic (error based) time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol=1.0e-8, reltol=1.0e-8, save_everystep=true);

# Create .gif animation of the solution
pyplot()
animation = @animate for k= 1:6:length(sol.t)
  pd = PlotData2D(sol.u[k], semi)
  wireframe(pd["H"])
  surface!(pd["b"], zlim=(38,65), camera = (30,20), title="t=$(sol.t[k])",
            xlabel="E", ylabel="N", zlabel="H")
end

gif(animation, "examples\\plots\\dam_break_2d.gif", fps=15)
```

This is the resulting .gif animation.

![gif](https://user-images.githubusercontent.com/101979498/203507057-f4fa5ef2-e852-493d-8df6-497c1e2a9a51.gif)

For the bottom topography, the domain's boundaries look weird. The reason for that is a bug in `PlotData2D` of `Trixi.jl`. Once this has been addressed, the plotted bottom topography will look similar to the one in [the previous section](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/function/#Two-dimensional-case).