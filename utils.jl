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
    linear_harmonic_oscillator_phase_portrait(;num_phasepoints=5, tf=10, h=0.01,xlims=(-1,1),ylims=(-1,1),xh=0.15, yh=0.15)

Plot the phase portrait of a linear harmonic oscillator, with randomly initialized trajectories   
# Arguments
- `num_phasepoints::Int`: number of trajectories  
- `tf::Float64`: Final time of the trajectories 
- `h::Float64`: Integration step size 
- `xlims::Tuple{Float64,Float64}`: Lower and upper x-limit 
- `ylims::Tuple{Float64,Float64}`: Lower and upper y-limit 
- `xh::Float64`: Grid step size in x-direction 
- `yh::Float64`: Grid step size in y-direction
"""
function linear_harmonic_oscillator_phase_portrait(;num_phasepoints=5, tf=10, h=0.01,xlims=(-1,1),ylims=(-1,1),xh=0.15, yh=0.15)
    # Solve numerically for some trajectories, given random initial conditions 
    traj = zeros(2num_phasepoints,Int(tf/h)+1)
    traj[:,1] = rand(2num_phasepoints)*1.25 # random initial conditions
    for k in eachindex(h:h:tf)
        for j in 1:num_phasepoints
            is = 2j-1
            ie = 2j
            traj[is:ie,k+1] = rk4(traj[is:ie,k],h,simple_harmonic_oscillator)
        end 
    end 
    # plot trajectories along with streamline plots 
    f = Figure(size = (700, 400),fontsize = 24)
    ax = Axis(f[1, 1], xlabel=L"$x$ (m)",ylabel=L"$v$ (m/s)", limits = (xlims,ylims))
    xxs,yys,xs,ys = coordinate_grid(xlims,ylims,xh,yh)
    sho(x) = Point2f(simple_harmonic_oscillator(x))
    strength = vec(sqrt.(xxs .^ 2 .+ yys .^ 2))
    heatmap!(ax, xxs, yys, strength, colormap = :Spectral, alpha=0.5)
    streamplot!(ax, sho, xs, ys, color=(p)->RGBf(0/255, 0/255, 0/255), gridsize= (20,20), arrow_size = 5, linewidth=0.5)
    for j in 1:num_phasepoints
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
    f = Figure(size = (1500, 300))
    cs = ColorScheme([colorant"yellow",colorant"fuchsia",colorant"deepskyblue", colorant"seagreen1"]);
    origcolors = get(cs,0:0.01:1)
    l = 0.8 # linewidth
    num_phasepoints = 20
    for j in eachindex(a)
        ax = Axis(f[1, j], title=titles[j])
        c = 1
        for k in 1:num_phasepoints
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

# uncoupled_linear_system_phase_portraits(a = [-1.5, -1, -0.25, 0, 0.5], tr = 0:0.001:4, titles::Vector{String} = ["a)","b)","c)","d)","e)"]) 