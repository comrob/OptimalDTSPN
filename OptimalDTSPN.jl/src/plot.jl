##################################################
# DTSPN Plot Functions
##################################################

function BNBSolver_plot_map(solver::BNBSolver)  
    plt.clf()
    plt.axis("equal")
    g_points = [x.center for x in solver.goals]
    plot_points(g_points, specs="rx", a = 0.8)
    for i in 1:length(solver.goals)
        c = solver.goals[i].center.coords
        plt.text(c[1]+0.05, c[2]-0.45, "\$R_{$i}\$")
    end
    for g in solver.goals
        plot_circle(g.center, g.radius)
    end  

end

function BNBSolver_plot_tour(solver, color)
    BNBSolver_plot_tour(solver, solver.feasible.ub_states, color)
end

function BNBSolver_plot_tour(solver, solution::Vector{Point}, color)
    n = length(solution)
    for a in 1:n
        b = (a % n) + 1
        s1 = solution[a]
        s2 = solution[b]
        gdip_init_dubins_maneuver(toArray(s1), toArray(s2), solver.radius)
        step_size = 0.01 * solver.radius
        if solver.radius == 0.0
            step_size = gdip_get_length() / 3.0
        end
        configuration = gdip_sample_many(step_size)
        plot_points(configuration, specs=color)
    end
end

function BNBSolver_plot_tour(solver, configurations::Vector{Vector{Point}}, color)
    if length(configurations) == 0
        return 
    end
    for configuration in configurations      
        plot_points(configuration, specs = color)
    end
end

function BNBSolver_plot_lb(solver, color)
    BNBSolver_plot_lb(solver, solver.lb.lb_samples, color)
end

function BNBSolver_plot_lb(solver, solution::Vector{Sample}, color)
    n = length(solution)
    if n == 0 
        return 0.0 
    end
    len = 0
    for a in 1:n
        b = (a % n) + 1
        s1 = solution[a]
        s2 = solution[b]
        lowerPathGDIP(s1, s2, solver.radius)
        step_size = solver.radius == 0 ? 0.1 : 0.01 * solver.radius
        configuration = gdip_sample_many(step_size)
        plot_points(configuration, specs=color)
    end
    return len
end

function BNBSolver_plot_samples(solver)
    if !isnothing(solver.lb.samples)
        for s in solver.lb.samples.samples
            for sample in s
                ax = BNB.plt.gca()
                Circle = matplotlib.patches.Circle
                circle = Circle(sample.center, sample.radius, facecolor=nothing ,edgecolor="green", linewidth=1, alpha=0.2)
                ax.add_patch(circle)
                # if sample.alphaResolution >= 8
                #     p1 = sample.center .+ solver.radius .* [cos(sample.alpha1), sin(sample.alpha1)]
                #     p2 = sample.center .+ solver.radius .* [cos(sample.alpha2), sin(sample.alpha2)]
                #     plot_points([sample.center, p1, p2, sample.center]; specs="r", a=0.5)
                # end
            end
        end
    end
end

function BNBSolver_plot_lb(solver, configurations::Vector{Vector{Point}}, color)
    for configuration in configurations      
        plot_points(configuration, specs = color)
    end
end

function BNBSolver_plot_save(solver, problem, filename)
    BNBSolver_plot_map(solver)
    BNBSolver_plot_tour(solver, "b-")
    BNBSolver_plot_lb(solver, "r--")
    BNBSolver_plot_samples(solver)

    gap = ( solver.feasible.ub - solver.lb.lb) /  solver.feasible.ub * 100.0
    s = @sprintf("Upper: %6.2f  Lower: %6.2f  GAP: %2.2f", solver.feasible.ub, solver.lb.lb, gap)
    plt.title("DTSPN $(problem): turning radius $(solver.radius), max_res:$(MAX_RES)\n$(s)")
    savefig(filename)
end

function plot_points(points::Vector{Point}; specs="b", a=1.0)
    x_val = [x.coords[1] for x in points]
    y_val = [x.coords[2] for x in points]
    p = plt.plot(x_val, y_val, specs, alpha=a, ms = 3, solid_capstyle="butt")	 	
end 

function plot_points(points; specs="b", a=1.0)
    x_val = [x[1] for x in points]
    y_val = [x[2] for x in points]
    p = plt.plot(x_val, y_val, specs, alpha=a, ms = 6, solid_capstyle="butt")	 	
end 

function plot_circle(xy, radius)
    ax = plt.gca()
    Circle = matplotlib.patches.Circle
    circle = Circle([xy.coords[1], xy.coords[2]], radius, facecolor="yellow", edgecolor="orange", linewidth=1, alpha=0.3)
    p = ax.add_patch(circle)
end
