#!/usr/bin/env julia

module Run

# Imports and includes
using OptimalDTSPN

using PyPlot
using Printf
using CSV
using TimerOutputs

# Use latex for texts
plt.rc("text", usetex=true)
plt.rc("font", family="serif")

# Settings
turning_radius = 1.0    # ρ in the paper
branching_factor = 1.0  # α in the paper
filename = "instances/Behdani/dtspn-20-02-100-0.txt"
max_time = 60           # Timelimit in seconds
name = string(split(filename, "/")[end][1:end-4])

# Save results into a directory "results"
results_dir="results-example"
mkpath(results_dir)

# Solve functions
function plot(solver)
    plt.figure(1)
    BNBSolver_plot_map(solver)
    BNBSolver_plot_tour(solver, "b-")
    BNBSolver_plot_lb(solver, "r--")
    #BNBSolver_plot_samples(solver)
    plt.axis("off")

    plt.figure(2)
    plt.clf()
    plt.plot(solver.times_eval, solver.lbs_eval, label = "LB", "--r")
    plt.plot(solver.times_eval, solver.ubs_eval, label = "UB", "-b")
    plt.xlabel("time [s]")
    plt.ylabel("Distance")
    plt.xscale("log")
    plt.legend()
    plt.title("Quality of the solution")
end

function solve(filename::String, problem::String, max_time::Int64, turning_radius::Float64, branching_factor)
    global to = TimerOutput()
    start = time()
    solver = solveBNB(filename, max_time, 0.0, turning_radius, branching_factor, to) do solver
        gap = (solver.feasible.ub - solver.lb.lb) /  solver.feasible.ub * 100.0
        printstyled(@sprintf("%.3f: ", time()-start), color=:blue) 
        print(@sprintf("%.2f <= %.2f GAP: %.2f%% O:%d G:%d\n", 
            solver.lb.lb,  solver.feasible.ub, gap, length(solver.O)+1, length(solver.dtrp_solvers)
        ))
        plot(solver)
    end
    plot(solver)
    to, solver
end

@info string("Solving $(name) in  '", filename, ", T = ", max_time, " [s]")
time_eval, solver = solve(filename, name, max_time, turning_radius, branching_factor)
    
gap = ( solver.feasible.ub - solver.lb.lb) /  solver.feasible.ub * 100.0
@info @sprintf("Problem SOLVED: Upper: %6.2f  Lower: %6.2f  GAP: %2.2f", solver.feasible.ub, solver.lb.lb, gap)
 
end
