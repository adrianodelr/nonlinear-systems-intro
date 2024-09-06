# used packages 
using LinearAlgebra
using CairoMakie
using LaTeXStrings
using Colors, ColorSchemes
using ForwardDiff

"""
    rk4(x,u,h,k)

Fourth order Runge Kutta scheme with ZOH on controls. 

# Arguments
- `x::Vector{Float64}`: System state 
- `u::Vector{Float64}`: Applied controls
- `h::Float64`: Step size for discretization 
- `f::Function`: Continuous time forward dynamic model as first order system   
"""            
function rk4(x,u,h,f)
    k1 = f(x,u); 
    k2 = f(x + (h/2)*k1,u);         
    k3 = f(x + (h/2)*k2,u);        
    k4 = f(x + h*k3,u);        
    xnew = x + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    return xnew
end

"""
    rk4(x,h,k)

Fourth order Runge Kutta scheme for autonomous system. 

# Arguments
- `x::Vector{Float64}`: System state 
- `h::Float64`: Step size for discretization 
- `f::Function`: Continuous time forward dynamic model as first order system   
"""            
function rk4(x,h,f)
    k1 = f(x); 
    k2 = f(x + (h/2)*k1);         
    k3 = f(x + (h/2)*k2);        
    k4 = f(x + h*k3);        
    xnew = x + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    return xnew
end

"""
    simple_harmonic_oscillator(x)

State space representation of simple harmonic oscillator in continuous time. It is artificially brought in first order form by a substitution, so we can use numerical tools as rk4.  
# Arguments
- `x::Vector{Float64}`: System state (position, velocity)
"""
function simple_harmonic_oscillator(x)
    x,v = x[1],x[2]
    k = 0.5     # spring constant 
    m = 1.0     # mass 
    return [v,-(k/m)*x]
end

"""
    coordinate_grid(xlims,ylims,xh,yh)

Returns a grid of coordinates (reshaped in vector form) for plotting. 
# Arguments
- `xlims::Tuple{Float64,Float64}`: Lower and upper x-limit 
- `ylims::Tuple{Float64,Float64}`: Lower and upper y-limit 
- `xh::Float64`: Grid step size in x-direction 
- `yh::Float64`: Grid step size in y-direction

# Returns 
- `xxs::Vector{Float64}`: X coordinates of a 2D grid
- `yys::Vector{Float64}`: Y coordinates of a 2D grid
- `xs::Vector{Float64}`: X coordinates range
- `ys::Vector{Float64}`: Y coordinates range
"""
function coordinate_grid(xlims,ylims,xh,yh)
    xs = xlims[1]:xh:xlims[2]
    ys = ylims[1]:yh:ylims[2]
    xxs = [x for x in xs for y in ys]
    yys = [y for x in xs for y in ys]
    return xxs,yys,xs,ys
end 


"""
    harmonic_oscillator_vector_field(;xlims=(-1,1),ylims=(-1,1),xh=0.15, yh=0.15)

Plot the vector field of a linear harmonic oscillator 
# Arguments
- `xlims::Tuple{Float64,Float64}`: Lower and upper x-limit 
- `ylims::Tuple{Float64,Float64}`: Lower and upper y-limit 
- `xh::Float64`: Grid step size in x-direction 
- `yh::Float64`: Grid step size in y-direction
"""
function linear_harmonic_oscillator_vector_field(;xlims=(-1,1),ylims=(-1,1),xh=0.15, yh=0.15)
    xxs,yys,xs,ys = coordinate_grid(xlims,ylims,xh,yh)
    f = Figure(size = (700, 400),fontsize = 24)
    ax = Axis(f[1, 1], xlabel=L"$x$ (m)",ylabel=L"$v$ (m/s)", limits = (xlims,ylims))
    arrows!(ax, xs, ys, simple_harmonic_oscillator , arrowsize = 5, lengthscale = 0.05, arrowcolor = "black", linecolor = "black", linewidth = 1)
    f
    return f 
end     

# linear_harmonic_oscillator_vector_field(xlims=(-2,2),ylims=(-1,1),xh=0.15, yh=0.15)


"""
    solve_trajectories_2D_sys(num_phasepoints, initial_conditions, system, h, tf)

Numerically solve the 2D differential equation system, given multiple initial condition. Stack resulting trajectories in a matrix   
# Arguments
- `num_traj::Int`: number of trajectories  
- `initial_conditions::Vector{Float}`: Vector of initial conditions of the form [x01,y01, x02,y02, ...]  
- `system::Function`: 2D continuous time system (autonomous) that takes in a 2D state  
- `h::Float64`: Integration step size 
- `tf::Float64`: Final integration 

# Returns 
- `traj::Matrix{Float64}`: Matrix where 2 rows represent the time evolution of the state 
"""
function solve_trajectories_2D_sys(num_traj, initial_conditions, system, h, tf)
    # Solve numerically for some trajectories, given random initial conditions 
    traj = zeros(2num_traj,Int(tf/h)+1)
    traj[:,1] = initial_conditions
    for k in eachindex(h:h:tf)
        for j in 1:num_traj
            is = 2j-1
            ie = 2j
            traj[is:ie,k+1] = rk4(traj[is:ie,k],h,system)
        end 
    end
    return traj 
end

"""
    linear_harmonic_oscillator_phase_portrait(;num_phasepoints=5, tf=10, h=0.01,xlims=(-1,1),ylims=(-1,1),xh=0.15, yh=0.15)

Plot the phase portrait of a linear harmonic oscillator, with randomly initialized trajectories   
# Arguments
- `num_traj::Int`: number of trajectories  
- `tf::Float64`: Final time of the trajectories 
- `h::Float64`: Integration step size 
- `xlims::Tuple{Float64,Float64}`: Lower and upper x-limit 
- `ylims::Tuple{Float64,Float64}`: Lower and upper y-limit 
- `xh::Float64`: Grid step size in x-direction 
- `yh::Float64`: Grid step size in y-direction
"""
function linear_harmonic_oscillator_phase_portrait(;num_traj=5, tf=10, h=0.01,xlims=(-1,1),ylims=(-1,1),xh=0.15, yh=0.15)
    # Solve numerically for some trajectories 
    initial_conditions = rand(2num_traj)*1.25 # random initial conditions
    traj = solve_trajectories_2D_sys(num_traj, initial_conditions, simple_harmonic_oscillator, h, tf)

    # plot trajectories along with streamline plots 
    f = Figure(size = (700, 400),fontsize = 24)
    ax = Axis(f[1, 1], xlabel=L"$x$ (m)",ylabel=L"$v$ (m/s)", limits = (xlims,ylims))
    xxs,yys,xs,ys = coordinate_grid(xlims,ylims,xh,yh)
    sho(x) = Point2f(simple_harmonic_oscillator(x))
    strength = vec(sqrt.(xxs .^ 2 .+ yys .^ 2))
    heatmap!(ax, xxs, yys, strength, colormap = :Spectral, alpha=0.5)
    streamplot!(ax, sho, xs, ys, color=(p)->RGBf(0/255, 0/255, 0/255), gridsize= (20,20), arrow_size = 5, linewidth=0.5)
    for j in 1:num_traj
        is = 2j-1
        ie = 2j
        lines!(ax, traj[is,:],traj[ie,:], color="black")
        scatter!(ax,[traj[is,1]],[traj[ie,1]], color="black")
    end
    f
    return f
end  

# linear_harmonic_oscillator_phase_portrait(num_phasepoints=5, tf=10, h=0.01,xlims=(-2,2),ylims=(-1,1),xh=0.10, yh=0.10)

"""
    linear_system_solution(t,x0;a=1)

Solution to the uncoupled linear first order system [ẋ,ẏ]^T = [a 0; 0 -1][x,y]^T

# Arguments
- `t::Float64`: time 
- `x0::Vector{Float64}`: initial condition 
- `a::Float64`: constant coefficient 
"""
function linear_system_solution(t,x0;a=1.0)
    return [x0[1]*exp(a*t);x0[2]*exp(-t)]
end 


"""
    uncoupled_linear_system_phase_portraits(;a = [-1.5, -1, -0.25, 0, 0.5], tr = 0:0.001:4, titles::Vector{String} = ["a)","b)","c)","d)","e)"])

Plot solutions to uncoupled linear first order systems [ẋ,ẏ]^T = [a 0; 0 -1][x,y]^T. Trajectories are determined from random initial conditins. 
Trajectories are mirrored at horizontal and vertical axis, such that trajectories 'start' in all quadrants. 

# Arguments
- `a::Vector{Float64}`: different coefficients for the linear system (number determines how many subplots) 
- `tr:StepRangeLen`: timerange for the trajectories  
- `titles::Vector{String}`: titles for the systems detemined by the coefficients a 
"""
function uncoupled_linear_system_phase_portraits(;a = [-1.5, -1, -0.25, 0, 0.5], tr = 0:0.001:4, titles::Vector{String} = ["a)","b)","c)","d)","e)"])
    if length(a) != length(titles)
        throw(error("number of titles not equal to number of systems (determined by lenght(a))"))
    end 
    f = Figure(size = (1500, 300),fontsize = 24)
    cs = ColorScheme([colorant"yellow",colorant"fuchsia",colorant"deepskyblue", colorant"seagreen1"]);
    origcolors = get(cs,0:0.01:1)
    l = 0.8 # linewidth
    num_traj = 20
    for j in eachindex(a)
        ax = Axis(f[1, j], title=titles[j])
        c = 1
        for k in 1:num_traj
            x0 = rand(2)*0.5
            xi = map(t-> linear_system_solution(t,[-x0[1],x0[2]],a=a[j]),tr) |> stack
            lines!(ax, xi[1,:],xi[2,:],color=origcolors[c], linewidth=l)
            c+=1
            xi = map(t-> linear_system_solution(t,[x0[1],-x0[2]],a=a[j]),tr) |> stack
            lines!(ax, xi[1,:],xi[2,:],color=origcolors[c], linewidth=l)
            c+=1
            xi = map(t-> linear_system_solution(t,[-x0[1],-x0[2]],a=a[j]),tr) |> stack
            lines!(ax, xi[1,:],xi[2,:],color=origcolors[c], linewidth=l)
            c+=1
            xi = map(t-> linear_system_solution(t,[x0[1],x0[2]],a=a[j]),tr) |> stack
            lines!(ax, xi[1,:],xi[2,:],color=origcolors[c], linewidth=l)
            c+=1
            scatter!(ax,[x0[1]],[x0[2]], markersize=3, color = "black")
            scatter!(ax,[-x0[1]],[x0[2]], markersize=3, color = "black")
            scatter!(ax,[x0[1]],[-x0[2]], markersize=3, color = "black")
            scatter!(ax,[-x0[1]],[-x0[2]], markersize=3, color = "black")
        end
    end
    f
    return f
end  

# uncoupled_linear_system_phase_portraits(a = [-1.5, -1, -0.25, 0, 0.5], tr = 0:0.001:4, titles = ["a)","b)","c)","d)","e)"]) 

"""
    nonlinear_system(x)

Some nonlinear 2D sample system

# Arguments
- `x::Vector{Float64}`: state vector 
"""
function nonlinear_system(x)
    return [x[1]+exp(-x[2]),-x[2]]
end  


"""
    nonlinear_system_direction_field(;initial_conditions = [-0.5,1, -0.9,1,-2,-1, -2,-0.5 ,-1,0 ,-1.3,0 ,-0.7,0 ,-0.63215,1, -1.5,-1], h = 0.01, tf = 6, xlims=(-3,3),ylims=(-1.5,1.5),xh=0.15, yh=0.15)

Plot the direction field of a nonlinear system, starting from given initial conditions
# Arguments
- `initial_conditions::Vector{Float}`: Vector of initial conditions of the form [x01,y01, x02,y02, ...]
- `tf::Float64`: Final time of the trajectories 
- `h::Float64`: Integration step size 
- `xlims::Tuple{Float64,Float64}`: Lower and upper x-limit 
- `ylims::Tuple{Float64,Float64}`: Lower and upper y-limit 
- `xh::Float64`: Grid step size in x-direction 
- `yh::Float64`: Grid step size in y-direction
"""
function nonlinear_system_direction_field(;initial_conditions = [-0.5,1, -0.9,1,-2,-1, -2,-0.5 ,-1,0 ,-1.3,0 ,-0.7,0 ,-0.63215,1, -1.5,-1], h = 0.01, tf = 6, xlims=(-3,3),ylims=(-1.5,1.5),xh=0.15, yh=0.15)    
    # solve numerically for trajectories 
    num_traj = Int(round(length(initial_conditions)/2))
    traj = solve_trajectories_2D_sys(num_traj, initial_conditions, nonlinear_system, h, tf)

    xxs,yys,xs,ys = coordinate_grid(xlims,ylims,xh,yh)
    ns(x) = Point2f(nonlinear_system(x))
    f = Figure(size = (700, 400),fontsize = 24)
    ax = Axis(f[1, 1], xlabel=L"$x_1$",ylabel=L"$x_2$", limits = (xlims,ylims))
    streamplot!(ax, ns, xs, ys, colormap = :Spectral, gridsize= (40,30), arrow_size = 0, linewidth=1)
    for j in 1:num_traj
        is = 2*j-1
        ie = 2*j
        lines!(ax, traj[is,:],traj[ie,:], color="black", linewidth=2)
        scatter!(ax,[traj[is,1]],[traj[ie,1]], color="black")
    end
    f
    return f
end  

# nonlinear_system_direction_field() 

"""
    double_well_potential(x)

Double well potential ODE

# Arguments
- `x::Vector{Float64}`: state vector 
"""
double_well_potential(x) = Point2f([x[2],x[1]-x[1]^3])

"""
    double_well_potential_energy(x,y)

Double well potential energy function

# Arguments
- `x::Vector{Float64}`: state vector 
"""
function double_well_potential_energy(x)
    x,y = x[1],x[2] 
    return  1/2*y^2 - 1/2*x^2 + 1/4*x^4
end  


"""
    energy_function_and_phase_plane_double_well(;initial_conditions = [0,0.01, -1,0.0, 1,0], h = 0.01, tf = 30, xlims=(-1.5,1.5),ylims=(-1.5,1.5),xh=0.05, yh=0.05)

Plot the energy function of a double well potential system together with its phase portrait and fixed points/ trajectories  
# Arguments
- `initial_conditions::Vector{Float}`: Vector of initial conditions of the form [x01,y01, x02,y02, ...]
- `tf::Float64`: Final time of the trajectories 
- `h::Float64`: Integration step size 
- `xlims::Tuple{Float64,Float64}`: Lower and upper x-limit 
- `ylims::Tuple{Float64,Float64}`: Lower and upper y-limit 
- `xh::Float64`: Grid step size in x-direction 
- `yh::Float64`: Grid step size in y-direction
"""
function energy_function_and_phase_plane_double_well(;initial_conditions = [0,0.01, -1,0.0, 1,0], h = 0.01, tf = 30, xlims=(-1.5,1.5),ylims=(-1.5,1.5),xh=0.05, yh=0.05)
    # bring energy function in format needed by makie and build gradient function 
    Ewrapper(x,y) = [double_well_potential_energy([x,y])] 
    gradE(xx, yy) = Point2f(ForwardDiff.jacobian(x -> Ewrapper(x[1], x[2]), [xx, yy])...) # gradient of energy funvtion 
    # discrete grid of energy values (in vector form)    
    x = xlims[1]:xh:xlims[2]
    y = ylims[1]:yh:ylims[2]
    z = [double_well_potential_energy([i, j]) for i in x, j in y];

    # solving numerically for trajectories    
    num_traj = Int(round(length(initial_conditions)/2))
    traj = solve_trajectories_2D_sys(num_traj, initial_conditions, double_well_potential, h, tf)
    
    # plotting
    zmin, zmax = minimum(z), maximum(z)
    cmap = :Spectral
    fig = Figure(fontsize = 24, size=(1200,500))
    ga = fig[1:2, 1]= GridLayout() 
    gb = fig[1:2, 2]= GridLayout() 
    # Energy function plot 
    ax1 = Axis3(ga[1,1],  perspectiveness = 0.3, elevation = 0.2, azimuth = 1.2,
    xlabel = L"$x$ (m)", ylabel = L"y (m)",zlabel = L"$E(x,y)$ (Joule)", 
    title=L"\textbf{a) Energy surface}", width=500, height=400)
    # 3d surface 
    surface!(ax1, x, y, z; colormap = cmap,colorrange = (zmin, zmax))
    contour3d!(ax1, x, y, z .+ 0.005; levels = 10, linewidth = 1, color = :white)
    wireframe!(ax1, x, y, z; color = (:black, 0.05))
    # projection on the floor
    heatmap!(ax1, x, y, z, colormap = cmap, alpha=1, transformation = (:xy, -0.3))
    contour!(ax1, x, y, z; levels = 10, linewidth = 1, color = :white, transformation = (:xy, -0.3))
    streamplot!(ax1, double_well_potential, x, y; color=(p)->RGBf(0/255, 0/255, 0/255), gridsize = (20, 20),arrow_size = 1, linewidth = 0.5, transformation = (:xy, -0.15))

    # phase portrait 
    ax2 = Axis(gb[1, 1], xlabel = L"$x$ (y)", ylabel = L"$y$ (m)", title=L"\textbf{b) Phase plane}", titlegap = 30, limits=(xlims,ylims))
    heatmap!(ax2, x, y, z, colormap = cmap, alpha=1.0)
    streamplot!(ax2, double_well_potential, x, y; color=(p)->RGBf(0/255, 0/255, 0/255), gridsize = (30, 30),arrow_size = 8, linewidth = 1)
    # fixed points 
    scatter!(ax2, [0.0], [0.0], color="white")
    scatter!(ax2, [-π,π], [0.0,0.0], color="blue")
    # trajectories 
    for j in 1:num_traj
        is = 2*j-1
        ie = 2*j
        lines!(ax2, traj[is,:],traj[ie,:], color=:blue, linewidth=2)
        scatter!(ax2,[traj[is,1]],[traj[ie,1]], color=:white)
    end

    colgap!(ga, 0)
    rowgap!(ga, 0)
    resize_to_layout!(fig)
    fig
    return fig 
end  

# energy_function_and_phase_plane_double_well()

"""
    single_pendulum(x)

single_pendulum ODE (assuming mass, link length etc = 1)

# Arguments
- `x::Vector{Float64}`: state vector 
"""
single_pendulum(x) = Point2f([x[2],-sin(x[1])])

"""
    single_pendulum_energy(x,y)

Single pendulum energy function (assuming mass, link length etc = 1)

# Arguments
- `x::Vector{Float64}`: state vector 
"""
function single_pendulum_energy(x)
    θ,v = x[1],x[2] 
    return  1/2*v^2 -cos(θ)
end  



"""
    energy_function_and_phase_plane_single_pendulum(;initial_conditions = [0,0.01, -1,0.0, 1,0], h = 0.01, tf = 30, xlims=(-1.5,1.5),ylims=(-1.5,1.5),xh=0.05, yh=0.05)

Plot the energy function of a double well potential system together with its phase portrait and fixed points/ trajectories  
# Arguments
- `initial_conditions::Vector{Float}`: Vector of initial conditions of the form [x01,y01, x02,y02, ...]
- `tf::Float64`: Final time of the trajectories 
- `h::Float64`: Integration step size 
- `xlims::Tuple{Float64,Float64}`: Lower and upper x-limit 
- `ylims::Tuple{Float64,Float64}`: Lower and upper y-limit 
- `xh::Float64`: Grid step size in x-direction 
- `yh::Float64`: Grid step size in y-direction
"""
function energy_function_and_phase_plane_single_pendulum(;initial_conditions = [0,0.01, -1,0.0, 1,0], h = 0.01, tf = 30, xlims=(-2π,2π),ylims=(-3,3),xh=0.1, yh=0.1)
    # bring energy function in format needed by makie and build gradient function 
    Ewrapper(x,y) = [single_pendulum_energy([x,y])] 
    gradE(xx, yy) = Point2f(ForwardDiff.jacobian(x -> Ewrapper(x[1], x[2]), [xx, yy])...) # gradient of energy funvtion 
    # discrete grid of energy values (in vector form)    
    x = xlims[1]:xh:xlims[2]
    y = ylims[1]:yh:ylims[2]
    z = [single_pendulum_energy([i, j]) for i in x, j in y];

    # solving numerically for trajectories    
    num_traj = Int(round(length(initial_conditions)/2))
    traj = solve_trajectories_2D_sys(num_traj, initial_conditions, single_pendulum, h, tf)
    
    # plotting
    zmin, zmax = minimum(z), maximum(z)
    cmap = :Spectral
    fig = Figure(fontsize = 24, size=(1200,500))
    ga = fig[1:2, 1]= GridLayout() 
    gb = fig[1:2, 2]= GridLayout() 
    # Energy function plot 
    ax1 = Axis3(ga[1,1],  perspectiveness = 0.3, elevation = 0.2, azimuth = 1.2,
    xlabel = L"$x$ (m)", ylabel = L"y (m)",zlabel = L"$E(x,y)$ (Joule)", 
    title=L"\textbf{a) Energy surface}", width=500, height=400)
    # 3d surface 
    surface!(ax1, x, y, z; colormap = cmap,colorrange = (zmin, zmax))
    contour3d!(ax1, x, y, z .+ 0.005; levels = 10, linewidth = 1, color = :white)
    wireframe!(ax1, x, y, z; color = (:black, 0.05))
    # projection on the floor
    heatmap!(ax1, x, y, z, colormap = cmap, alpha=1, transformation = (:xy, -1.4))
    contour!(ax1, x, y, z; levels = 10, linewidth = 1, color = :white, transformation = (:xy, -1.4))
    streamplot!(ax1, single_pendulum, x, y; color=(p)->RGBf(0/255, 0/255, 0/255), gridsize = (20, 20),arrow_size = 1, linewidth = 0.5, transformation = (:xy, -0.7))

    # phase portrait 
    ax2 = Axis(gb[1, 1], xlabel = L"$x$ (y)", ylabel = L"$y$ (m)", title=L"\textbf{b) Phase plane}", titlegap = 30, limits=(xlims,ylims))
    heatmap!(ax2, x, y, z, colormap = cmap, alpha=1.0)
    streamplot!(ax2, single_pendulum, x, y; color=(p)->RGBf(0/255, 0/255, 0/255), gridsize = (30, 30),arrow_size = 8, linewidth = 1)
    # fixed points 
    scatter!(ax2, [0.0], [0.0], color="white")
    scatter!(ax2, [-π,π], [0.0,0.0], color="blue")

    colgap!(ga, 0)
    rowgap!(ga, 0)
    resize_to_layout!(fig)
    fig
    return fig 
end  

energy_function_and_phase_plane_single_pendulum()