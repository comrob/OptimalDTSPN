##################################################
# Structures - Branch and Bound
##################################################

mutable struct PartialSolution
    Σ::Vector{Int64}                 # Σ - sequence
    lb::Float64                      # L(Tlb) - lower bound
    ub::Float64                      # L(Tub) - upper bound
    lb_samples::Vector{Sample}       # Tlb - lower bound path 
    ub_states::Vector{Point}         # Tub - upper bound path
    samples::Union{Sampling,Nothing}
end

PartialSolution() = PartialSolution([], 0, Inf, [],[], nothing)

Base.isless(a::PartialSolution, b::PartialSolution) = a.lb < b.lb

# Node in the search tree
mutable struct Node
    partial::PartialSolution
    feasible::PartialSolution
    ω::Int64
    parent::Union{Node, Nothing}
    childs::Vector{Node}
end

Node() = Node(PartialSolution(),PartialSolution())
Node(partial, ub) = Node(partial, ub, RES, nothing, [])

struct NotCovered
    label::Int64
    feasib::Bool
    dist::Float64
end

Base.:<(a::NotCovered, b::NotCovered) = a.dist < b.dist

##################################################
# The main solver class
##################################################
mutable struct BNBSolver{T}
    goals::Vector{T}
    radius::Float64

    feasible::PartialSolution
    lb::PartialSolution

    O::PriorityQueue{Node, Float64}

    dtrp_solvers::Dict{Vector{Int}, GDIPSolver{T}}

    times_eval::Vector{Float64}
    lbs_eval::Vector{Float64}
    ubs_eval::Vector{Float64}

    init_time::Float64
end

function BNBSolver(filename::String, sensing_radius::Float64, radius::Float64)
    targets = load(filename, sensing_radius)
    n = length(targets)
    T = typeof(first(targets))
    BNBSolver(
        targets, radius, 
        PartialSolution(), PartialSolution(), 
        PriorityQueue{Node, Float64}(), 
        Dict{Vector{Int}, GDIPSolver{T}}(), 
        Float64[], Float64[], Float64[],
        0.0
    )
end 

function log_outputs(solver, update_fce)
    now_time = (time()  - solver.init_time)
    push!(solver.times_eval, now_time)
    push!(solver.lbs_eval, solver.lb.lb)
    push!(solver.ubs_eval, solver.feasible.ub)
    update_fce(solver)
end

##################################################
# Nice printing of the search tree
##################################################

function AbstractTrees.children(x::Node) 
    ret = [a for a in x.childs]
    # Direct edge with only one child
    while length(ret) == 1 && length(ret[1].childs) > 0
        ret = [a for a in ret[1].childs]
    end
    # Direct edge for only one child
    for i in 1:length(ret)
        while length(ret[i].childs) == 1
            ret[i] = first(ret[i].childs)
        end
    end
    return ret
end

##################################################
# DTSPN Methods
##################################################

# Returns Σroot
function SelectRoot(solver)
    n = length(solver.goals)
    d = distance_matrix(solver.goals)
    argmax([i,j,k] for i=1:n, j=1:n, k=1:n) do (i,j,k)
        (j <= i || k <= j) && return -Inf
        d[i,j] + d[j,k] + d[k,i]
    end
end

function ComputeNotCovered(solver, solution::PartialSolution)
    targets = solver.goals
    n = length(targets)
    dists = [NotCovered(i, false, Inf) for i in 1:n]
    for s in solution.Σ
        dists[s] = NotCovered(s, true, -Inf)
    end

    lb_samples = solution.lb_samples
    m = length(lb_samples)
    for act in 1:n
        target = targets[act]
        best_dist = Inf
        for a=1:m
            b = mod1(a+1, m)
            lowerPathGDIP(lb_samples[a], lb_samples[b], solver.radius)
            state = gdip_dubins_closest(target.center.coords)
            act_dist = dist(state[1:2], target.center.coords) - target.radius
            if act_dist < best_dist
                best_dist = act_dist
            end
        end
        if best_dist <= dists[act].dist
            dists[act] = NotCovered(act, best_dist <= 0.0, best_dist)
        end
    end

    return filter(x -> !x.feasib, dists)
end

function optimize_LIO(ub_points, targets::Vector{TargetRegion}, radius)
    n = length(ub_points)
    @assert n == length(targets)

    d(a,b) = Dubins(a,b,radius)

    TOL = 1e-7

    for i = 1:3n
        target = targets[mod1(i, n)]

        prev = ub_points[mod1(i-1, n)]
        act  = ub_points[mod1(i, n)]
        next = ub_points[mod1(i+1, n)]

        step1 = pi
        step2 = 0.5
        act_dist = d(prev, act) + d(act, next)
        while abs(step1) > TOL || abs(step2) > TOL
            # Part 1 - rotation
            new_act = Point(act.coords, act.theta + step1)
            new_dist = d(prev, new_act) + d(new_act, next)
            if new_dist < act_dist
                act = new_act
                act_dist = new_dist
                step1 *=2
            else
                step1 *= -0.1
            end
            # Part 2 - position
            vector = act.coords - target.center.coords
            #@show vector
            vector = [[cos(step2),  -sin(step2)] [sin(step2), cos(step2)]] * vector
            if dist(vector) > target.radius
                vector /= (dist(vector) / target.radius)
            end
            new_act = Point(target.center.coords + vector, act.theta)
            new_dist = d(prev, new_act) + d(new_act, next)
            if new_dist < act_dist
                act = new_act
                act_dist = new_dist
                step2 *=2
            else
                step2 *= -0.1
            end

            if abs(step1) < TOL
                step1 = sign(step1) * TOL
            end
            if abs(step2) < TOL
                step2 = sign(step2) * TOL
            end
        end

        ub_points[mod1(i, n)] = act

        err = dist(act, target.center) - target.radius
        @assert err < TOL
    end
end

function compute_upper_bound(solver, targets::Vector{TargetRegion}, node::Node)     
    seq = copy(node.partial.Σ)
    ub_points = copy(node.partial.ub_states)

    updated = true
    while updated
        updated = false
        notCovered = [x for x in 1:length(targets) if !(x in seq)]
        isempty(notCovered) && break

        dists = fill(Inf, length(notCovered))
        states = fill(Point(), length(notCovered))
        idx = fill(-1, length(notCovered))

        n = length(seq)
        for a=1:n
            b = mod1(a+1, n)
            Dubins(ub_points[a], ub_points[b], solver.radius)
            for i in 1:length(notCovered)
                target = targets[notCovered[i]]
                state = gdip_dubins_closest(target.center.coords)
                d = dist(state[1:2], target.center.coords) - target.radius
                if d < dists[i]
                    dists[i] = d
                    states[i] = PointFromDubins(state)
                    idx[i] = a
                end
            end
        end

        selected_id = argmax(dists)            
        selected_dist = dists[selected_id]
        if selected_dist > 0.
            # @info("Target $(act) will be included at $(best_idx)")
            best_state = states[selected_id]
            target = targets[notCovered[selected_id]]
            best_idx = idx[selected_id]
            d = dist(best_state.coords[1:2], target.center.coords)
            vec = (best_state.coords - target.center.coords)/d*target.radius
            new_coords = target.center.coords .+ vec
            new_state = Point(new_coords, best_state.theta)

            insert!(seq, best_idx+1, notCovered[selected_id])
            insert!(ub_points, best_idx+1, new_state)
            updated = true
        end
        if updated
            # LIO optimazation
            optimize_LIO(ub_points, targets[seq], solver.radius)
        end
    end

    ub = tour_length(ub_points, solver.radius)
    return PartialSolution(seq, node.partial.lb, ub, [], ub_points, nothing)
end

function print_Σ_tree(solver, root)
    function clean_tree(node)
        for n in node.childs
            clean_tree(n)
        end
        empty!(node.childs)
    end
    clean_tree(root)

    function add_to_tree(node)
        if !isnothing(node.parent)
            push!(node.parent.childs, node)
            add_to_tree(node.parent)
        end
    end
    for o in solver.O
        add_to_tree(o[1])
    end

    function custom_print(io, x::Node)
        color = isempty(x.childs) ? :green : :gray
        printstyled(stdout, x.ω, ", ", x.partial.lb, color=color)
        print(io, x.ω, ", ", x.partial.lb, ", ", x.partial.Σ)
    end
    print_tree(custom_print, stdout, root; maxdepth = 20)
end

function compute_bounds(solver, Σ, parent, fce, in_time_limit; ω::Int64 = RES, new_id = -1, new_goal = nothing)
    node::Node = Node(PartialSolution(copy(Σ), 0.0, Inf, [], [], nothing), 
                PartialSolution(copy(Σ), 0.0, Inf, [], [], nothing),
                ω, parent, Node[])

    !isnothing(parent) && push!(parent.childs, node)

    dtrp_solver = get(solver.dtrp_solvers, node.partial.Σ, nothing)
    if isnothing(dtrp_solver)
        dtrp_solver_parent = nothing
        if !isnothing(parent)
            dtrp_solver_parent = get(solver.dtrp_solvers, parent.partial.Σ, nothing)
        end
        if isnothing(dtrp_solver_parent) || new_id == -1
            dtrp_solver = GDIPSolver_init(solver.radius, solver.goals[node.partial.Σ])
        else
            dtrp_solver = GDIPSolver_insert_target(dtrp_solver_parent, new_id, new_goal, solver.radius)
        end
    end

    solver.dtrp_solvers[node.partial.Σ] = dtrp_solver

    refined = refine_dtrp(dtrp_solver, node.ω, in_time_limit)
    isnothing(refined) && return nothing

    @timeit to "Refine" (lb, ub, lb_samples, ub_samples, ub_states, sampling) = refined
                        

    node.partial.lb = lb
    node.partial.ub = ub
    node.partial.lb_samples = lb_samples
    node.partial.ub_states = PointFromDubins.(ub_states) 
    # TODO
    optimize_LIO(node.partial.ub_states, solver.goals[node.partial.Σ], solver.radius)
    new_ub = tour_length(node.partial.ub_states, solver.radius)
    #@show 1-new_ub/node.partial.ub
    node.partial.ub = new_ub
    node.partial.samples = sampling  

    @timeit to "ComputeUB" node.feasible = compute_upper_bound(solver, solver.goals, node)

    if node.partial.lb > node.feasible.ub
        node.partial.lb = node.feasible.ub
        @warn "Bounds intersect"
    end
    #@assert node.partial.lb <= node.feasible.ub + 10 * 1e-5 "Very bad error!"

    if node.feasible.ub < solver.feasible.ub  
        solver.feasible = node.feasible
        log_outputs(solver, fce)
    end

    return node
end

##################################################
# Main solver
##################################################
function solveBNB(fce, filename::String, max_time::Int64, 
    sensing_radius::Float64, radius::Float64,
    branching_factor::Float64, tio = TimerOutput();
    Σroot = nothing
)
    global to = tio  

    h(sol) = (sol.ub * branching_factor + sol.lb * (1-branching_factor))

    solver = BNBSolver(filename, sensing_radius, radius)
    solver.init_time = time()
    root::Node = Node()

    @timeit to "Root" begin
        if isnothing(Σroot)
            @notimeit to Σ  = SelectRoot(solver)
        else
            Σ = Σroot
        end
        printstyled(stdout, "Root selected ", Σ, "\n", color=:red, bold = true)
        @notimeit to root = compute_bounds(solver, Σ, nothing, fce, nothing)
        solver.lb = root.partial; 
        solver.feasible = root.feasible
        log_outputs(solver, fce)
    end

    in_time_limit() = (time() - solver.init_time) <= max_time

    solver.O[root] = root.partial.lb
    @timeit to "BNB" while in_time_limit()
        @assert !isempty(solver.O) "O list cannot be empty!"

        cur::Node = dequeue!(solver.O)   
        @timeit to "NotCovered" notCovered = ComputeNotCovered(solver, cur.partial)
        if cur.partial.lb > solver.lb.lb
            solver.lb = cur.partial
            #log_outputs(solver, fce)
        end
        
        if isempty(notCovered) && cur.ω >= MAX_RES
            log_outputs(solver, fce)
            return solver # The computation is over
        end

        resolution_branching = h(cur.partial) > solver.feasible.ub
        if isempty(notCovered) || (resolution_branching && cur.ω < MAX_RES)
            #log_outputs(solver, fce)
            @timeit to "Bound1" node = compute_bounds(solver, cur.partial.Σ, cur, fce, in_time_limit; ω = 2*cur.ω)
            if isnothing(node)
                @warn "Terminated during Bound1"                
                return solver
            end
            solver.O[node] = node.partial.lb
        else
            r::Int = notCovered[findmax([x.dist for x in notCovered])[2]].label
            @assert r ∉ cur.partial.Σ
            added = 0
            for i in 1:length(cur.partial.Σ)
                Σ = insert!(copy(cur.partial.Σ), i, r)
                @timeit to "Bound2" child = compute_bounds(solver, Σ, cur, fce, in_time_limit; new_id = i, new_goal = solver.goals[r])     
                if isnothing(child)
                    @warn "Terminated during Bound2"                
                    return solver
                end
                if child.partial.lb < solver.feasible.ub
                    solver.O[child] = child.partial.lb  
                    added += 1
                end
            end
            printstyled(stdout, "Branching $(added)/$(length(cur.partial.Σ)) child nodes added (+$(added-1))\n", color=:yellow)
        end
        l1 = length(solver.O)
        filter!(o->o[1].partial.lb <= solver.feasible.ub, solver.O)
        #print_Σ_tree(solver, root)
        l2 = length(solver.O)
        if l1 != l2
            printstyled(stdout, "Filter out $(l1-l2)/$(l2) nodes (-$(l1-l2))\n", color=:green)
        end
        log_outputs(solver, fce) # Just for sure, safe results one more time
    end
    log_outputs(solver, fce) # Just for sure, safe results one more time
    return solver
end
