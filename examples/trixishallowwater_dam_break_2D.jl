###############################################################################
# This example file creates a .gif file of a two dimensional dam break        #
# where the bottom topography is chosen from the DGM data set.                #
###############################################################################

# Include packages
using TrixiBottomTopography
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

# Download two dimensional Rhine bottom data from gist
root_dir = pkgdir(TrixiBottomTopography)
Rhine_data = Trixi.download(
    "https://gist.githubusercontent.com/maxbertrand1996/a30db4dc9f5427c78160321d75a08166/raw/fa53ceb39ac82a6966cbb14e1220656cf7f97c1b/Rhine_data_2D_40.txt",
    joinpath(root_dir, "examples", "Rhine_data_2D_40.txt"),
)

# B-spline interpolation of the underlying data
spline_struct = BicubicBSpline(Rhine_data)
spline_func(x, y) = spline_interpolation(spline_struct, x, y)

equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 55.0)

function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)
    inicenter = SVector(357490.0, 5646519.0)
    x_norm = x - inicenter
    r = sqrt(x_norm[1]^2 + x_norm[2]^2)

    # Calculate primitive variables
    H = r < 50 ? 65.0 : 55.0
    v1 = 0.0
    v2 = 0.0

    x1, x2 = x
    b = spline_func(x1, x2)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_wave

# Setting the boundary to be a free-slip wall
boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(
    polydeg = 3,
    surface_flux = (flux_fjordholm_etal, flux_nonconservative_fjordholm_etal),
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
)

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])
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
# ODE solvers, callbacks etc.

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# run the simulation

# define equidistant nodes in time for visualization of an animation
visnodes = range(tspan[1], tspan[2], length = 175)

# use a Runge-Kutta method with error based time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8, saveat = visnodes);

# Create an animation of the solution
if isdefined(Main, :Makie)
    j = Makie.Observable(1)
    time = Makie.Observable(0.0)

    pd_list = [PlotData2D(sol.u[i], semi) for i = 1:length(sol.t)]
    f = Makie.Figure()

    title_string = lift(t -> "time t = $(round(t, digits=3))", time)
    az = 130 * pi / 180
    el = 18 * pi / 180
    ax = Makie.Axis3(
        f[1, 1],
        xlabel = "E",
        ylabel = "N",
        zlabel = "H",
        title = title_string,
        azimuth = az,
        elevation = el,
    )

    height = lift(i -> pd_list[i].data[1], j)
    bottom = lift(i -> pd_list[i].data[4], j)
    Makie.surface!(ax, pd_list[1].x, pd_list[1].y, bottom; colormap = :greenbrownterrain)
    Makie.wireframe!(
        ax,
        pd_list[1].x,
        pd_list[1].y,
        height;
        color = Makie.RGBA(0, 0.5, 1, 0.4),
    )
    Makie.zlims!(ax, 35, 70)

    Makie.record(f, "animation_2d.gif", 1:length(pd_list)) do tt
        j[] = tt
        time[] = sol.t[tt]
    end
end
