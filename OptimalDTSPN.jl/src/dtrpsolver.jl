#!/usr/bin/env julia

using GeneralizedDubinsIntervalProblem

INIT_RESOLUTION = 8
OPTIMIZED_DTRP = true

##################################################
# Sample on the given target region
##################################################
# T - type of target, P - problem 
mutable struct Sample{T, P}
    target::T
    center::Array{Float64}
    radius::Float64
    # heading interval [alpha1, alpha2]
    # alpha2 - alpha1 = 2π / alphaResolution
    alpha1::Float64
    alpha2::Float64
    alphaResolution::Int64
    # position interval on the boundary [beta1, beta2]
    # beta2 - beta1 = 2π / betaResolution
    beta1::Float64
    beta2::Float64
    betaResolution::Int64
end

Base.copy(s::Sample{T, P}) where T where P = Sample{T, P}(
        s.target, s.center, s.radius, 
        s.alpha1, s.alpha2, s.alphaResolution, 
        s.beta1, s.beta2, s.betaResolution
    )

function Sample_init(targetRegion::T, P) where T
    alpha1 = 0.
    alpha2 = 2 * pi
    alphaResolution = 1
    beta1 = 0.
    beta2 = 2 * pi
    betaResolution = 1
    return Sample{T, P}(
        targetRegion, targetRegion.center.coords, targetRegion.radius, 
        alpha1, alpha2, alphaResolution, 
        beta1, beta2, betaResolution
    )
end

"""Split the actual sample into two new ones.
    The first is stored directly in the actual sample, and the second is returned.
    If the required resolution is already met, then nothing is done and None is returned.
    
    Parameters:
        resolution: the requred resolution
    Returns:
        Sample - the second new sample
"""
function Sample_split(sample::Sample{A,B}, resolution::Int64) where {A,B}
    # test wheather the sample may be heading into the circle region
    function test_incoming(sam::Sample)
        if sam.alphaResolution < 4
            return true
        end
        diff = 2π/sam.betaResolution
        lim1 = sam.beta1 + pi/2
        TOL = 1e-5
        a1 = mod2pi(sam.alpha1 - lim1 + TOL)
        a2 = mod2pi(sam.alpha2 - lim1 + TOL)
        #return true
        return a1 < π + diff + 2TOL || a2 < π + diff + 2TOL
    end

    # prefer splitting according position resolution
    if sample.betaResolution < resolution && sample.target.radius > 0.0
        sam1 = copy(sample)
        sam2 = copy(sample)
        sam1.betaResolution = sam2.betaResolution = 2 * sample.betaResolution
        sam1.beta2 = sam2.beta1 = (sample.beta1 + sample.beta2) / 2
        update_center_radius(sam1)
        update_center_radius(sam2)
        return filter(test_incoming, Sample[sam1, sam2])
    end
    if sample.alphaResolution < resolution && B == :DTSPN
        sam1 = copy(sample)
        sam2 = copy(sample)
        sam1.alphaResolution = sam2.alphaResolution = 2 * sample.alphaResolution
        sam1.alpha2 = sam2.alpha1 = (sample.alpha1 + sample.alpha2) / 2
        return filter(test_incoming, Sample[sam1, sam2])
    end
    return Sample[]
end

function update_center_radius(self)
    p1 = get_position_at_boundary(self.target, self.beta1)
    p2 = get_position_at_boundary(self.target, self.beta2)
    self.center = ((p1 + p2) / 2).coords
    self.radius = dist(p1.coords, p2.coords) / 2
end

function getFeasibleState(self)
    pos = get_position_at_boundary(self.target, self.beta1)
    q = zeros(3)
    q[1:2] = pos.coords
    q[3] = self.alpha1
    return q
end

##################################################
# Sampling structure which holds all the used samples
##################################################
mutable struct Sampling
    goals
    samples
end

function Sampling_init(goals, P)
    samples = [[Sample_init(t, P)] for t in goals]
    return Sampling(goals, samples)
end

"""Refine the seleted samples if the required resolution is not met.

    Parameters:
        slected: indexes of the selected samples (vector 1 x n)
        resolution: the requred resolution
    Returns:
        boolean - true if any sample is refined
    """
function refine_samples(self, selected, resolution)
    n = length(self.samples)
    refined = false
    for i in 1:n
        to_split = selected[i]
        samp = self.samples[i][to_split]
        res = Sample_split(samp, resolution)
        if !isempty(res)
            self.samples[i][to_split] = res[1]
            if length(res) > 1
                push!(self.samples[i], res[2])
            end
            refined = true 
        end
    end
    return refined
end

##################################################
# Functions
##################################################

"""Compute lower-bound path using GDIP between two configurations

    Arguments:  
        s1 - start; s2 - end; turning_radius
    Returns:    
        Dubins maneuver (DubinsWrapper)
"""
function lowerPathGDIP(s1::Sample{A,B}, s2::Sample{A,B}, turning_radius) where {A,B}
    interval1 = [s1.alpha1, s1.alpha2 - s1.alpha1]
    interval2 = [s2.alpha1, s2.alpha2 - s2.alpha1]
    if B == :DTSPN
        # @timeit to "gdip" begin
            return gdip_init_gdip(s1.center, interval1, s1.radius, s2.center, interval2, s2.radius, turning_radius)
        # end
    else
        return @timeit to "lb euclid" dist(s1.center - s2.center) - s1.radius - s1.radius
    end
end

"""Compute feasible Dubins path two configurations

    Arguments:  
        s1 - start; s2 - end; turning_radius
    Returns:    
        (Dubins maneuver 'DubinsWrapper', length)
"""
function upperPathGDIP(s1::Sample{A,B}, s2::Sample{A,B}, turning_radius) where {A,B}
    q1 = getFeasibleState(s1)
    q2 = getFeasibleState(s2)
    if B == :DTSPN
        ub = gdip_init_dubins_maneuver(q1, q2, turning_radius)
    else
        ub = dist(q1, q2)
    end
    return ub
end

function compute_distance_matrix(samples, idx, dst_fce, turning_radius)
    n = length(samples)
    samples1 = samples[mod1(idx,n)]
    samples2 = samples[mod1(idx+1,n)]
    n1 = length(samples1)
    n2 = length(samples2)
    return [dst_fce(samples1[i1], samples2[i2], turning_radius) for i1=1:n1, i2=1:n2]
end

function compute_distances(samples, dst_fce; turning_radius = 0)
    n = length(samples)
    distances = Vector{Matrix{Float64}}()
    resize!(distances, n)
    for i in 1:n
        distances[i] = compute_distance_matrix(samples, i, dst_fce, turning_radius)
    end
    return distances
end

function update_distances(samples, distances, selected_samples, dst_fce; turning_radius = 0)
    n = length(samples)    
    for i in 1:n
        ss1 = samples[i]
        ss2 = samples[(i % n) + 1]

        n1 = length(ss1)
        n2 = length(ss2)

        sh_old::Matrix{Float64} = distances[i]
        siz_old = size(sh_old) 
        sh = fill(NaN, (n1, n2))
        sh[1:siz_old[1],1:siz_old[2]] = sh_old

        sh[selected_samples[i], :] .= NaN 
        sh[:,selected_samples[(i%n)+1]] .= NaN

        for i1 = [selected_samples[i], siz_old[1]+1:n1...], i2 in 1:n2
            @timeit to "dist_fce" begin
                sh[i1,i2] = dst_fce(ss1[i1], ss2[i2], turning_radius)
            end
        end
        for i1 in 1:n1, i2 in [selected_samples[(i%n)+1], siz_old[2]+1:n2...]
            dist = sh[i1,i2]
            if isnan(dist)
                @timeit to "dist_fce" begin
                    sh[i1,i2] = dst_fce(ss1[i1], ss2[i2], turning_radius)
                end
            end
        end
        distances[i] = sh
    end
    return distances
end

function find_shortest_tour(distances)
    n = size(distances)[1]
    best_len = Inf
    best_tour = Int64[]

    # maximal number of samples
    k_max = maximum([size(x)[1] for x in distances])
    no_start = size(distances[1])[1]
    #@show no_start, k_max, n
    for start in 1:no_start
        # shortest sh[region_idx][sample_idx]
        sh = fill(Inf, (n+1, k_max))
        # used edge
        prev = fill(-1, (n+1, k_max))

        sh[1,start] = 0.0

        @inbounds for region_idx in 1:n
            (n1, n2) = size(distances[region_idx])
            for idx2 in 1:n2
                min_sh = Inf
                min_prev = -1
                for idx3 = 1:n1
                    dst = sh[region_idx,idx3] + distances[region_idx][idx3,idx2]
                    if min_sh > dst
                        min_sh = dst
                        min_prev = idx3
                    end
                end
                sh[region_idx+1,idx2] = min_sh
                prev[region_idx+1,idx2]= min_prev
            end
        end

        act_sol = sh[n+1,start]
        if act_sol < best_len
            best_len = act_sol
            tour = Int64[]
            act = start
            for i in 1:n
                act = prev[n-i+2,act]
                push!(tour, act)
            end
            best_tour = reverse(tour)
        end
    end

    return best_tour::Vector{Int64}
end

function find_shortest_tour_simplified(distances, extra_cycles = 1)
    n = size(distances)[1]
    best_len = Inf
    best_tour = Int64[]

    # maximal number of samples
    k_max = maximum([size(x)[1] for x in distances])

    M = (1 + 2*extra_cycles)n

    # shortest sh[region_idx][sample_idx]
    sh = fill(Inf, (M+1, k_max))
    # used edge
    prev = fill(-1, (M+1, k_max))

    sh[1,:] .= 0.0

    @inbounds for region_idx in 1:M
        (n1, n2) = size(distances[mod1(region_idx,n)])
        for idx2 in 1:n2
            min_sh = Inf
            min_prev = -1
            for idx3 = 1:n1
                dst = sh[region_idx,idx3] + distances[mod1(region_idx,n)][idx3,idx2]
                if min_sh > dst
                    min_sh = dst
                    min_prev = idx3
                end
            end
            sh[region_idx+1,idx2] = min_sh
            prev[region_idx+1,idx2]= min_prev
        end
    end

    start = findmin(sh[M+1,1])[2]
    act_sol = sh[M+1,start]
    if act_sol < best_len
        best_len = act_sol
        tour = Int64[]
        act = start
        for i in 1:M
            act = prev[M-i+2,act]
            push!(tour, act)
        end
        best_tour = reverse(tour)
    end

    optimal =  best_tour[n*extra_cycles] == best_tour[n*(1+extra_cycles)]
    best_tour = best_tour[(n*extra_cycles+1):(n*(1+extra_cycles))]

    # if !optimal
    #     @warn "Not optimal $extra_cycles"
    # end
    
    return best_tour::Vector{Int64}, optimal
end

function retrieve_path(samples, dst_fce, turning_radius, selected_samples)
    n = length(samples)
    path = []
    for a in 1:n
        g1 = samples[a][selected_samples[a]]
        g2 = samples[(a+1) % n][selected_samples[(a+1) % n]]
        path.append(dst_fce(g1, g2, turning_radius))
    end
    return path
end

path_len(path) = sum_list([dub.get_length() for dub in path])
get_position_at_boundary(region, beta) = region.center + region.radius * Point([cos(beta), sin(beta)])

##################################################
# The main solver class
##################################################
mutable struct GDIPSolver{T}
    turning_radius::Float64
    sampling::Sampling
    goals::Vector{T}
    resolution::Int

    lower_path::Vector{Integer}
    upper_path::Vector{Integer}

    lower_bound::Float64
    upper_bound::Float64

    lb_dist
end

function GDIPSolver_init(turning_radius, goals)
    type = turning_radius == 0 ? :CETSP : :DTSPN
    return GDIPSolver(
        turning_radius, Sampling_init(goals, type), goals, INIT_RESOLUTION, 
        Integer[], Integer[], 
        0., Inf,
        nothing
    )
end

function GDIPSolver_insert_target(solver::GDIPSolver, new_id, new_goal, turning_radius)
    goals = copy(solver.goals)
    insert!(goals, new_id, new_goal)
    type = turning_radius == 0 ? :CETSP : :DTSPN
    sampling = deepcopy(solver.sampling)
    sampling.goals = goals
    insert!(sampling.samples, new_id, [Sample_init(new_goal, type)])
    lb_dist = nothing
    dst_fce = lowerPathGDIP
    if !isnothing(solver.lb_dist)
        lb_dist = deepcopy(solver.lb_dist)
        insert!(lb_dist, new_id, compute_distance_matrix(sampling.samples, new_id, dst_fce, turning_radius))
        lb_dist[mod1(new_id-1, length(goals))] = 
            compute_distance_matrix(sampling.samples, new_id-1, dst_fce, turning_radius)
    end
    #@warn "TODO Can be optimized - we have parent"
    return GDIPSolver(
        solver.turning_radius, sampling, goals, INIT_RESOLUTION, 
        Integer[], Integer[], 
        0., Inf,
        lb_dist
    )
end

function GDIPSolver_return_length(solver, selected_samples, maneuver_function)
    sampling = solver.sampling
    n = length(solver.sampling.samples)
    len = 0.
    configurations = Sample[]
    for a in 1:n
        b = (a%n)+1
        g1 = sampling.samples[a][selected_samples[a]]
        g2 = sampling.samples[b][selected_samples[b]]
        len += maneuver_function(g1, g2, solver.turning_radius)
        push!(configurations, g1)
    end
    return len, configurations
end

"""Plot the actual sampling, lower and upper bound path

    Returns:
        (double, double) - lower bound, upper bound
"""
function GDIPSolver_return_bounds(solver; lb_dist = nothing)
    solver.lower_path = GDIPSolver_find_shortest_tour(solver, lowerPathGDIP; dist = lb_dist, optimal = true)[1]
    solver.upper_path = GDIPSolver_find_shortest_tour(solver, upperPathGDIP, optimal = true)[1]

    lb, lb_samples = GDIPSolver_return_length(solver, solver.lower_path, lowerPathGDIP)
    ub, ub_samples = GDIPSolver_return_length(solver, solver.upper_path, upperPathGDIP)
    ub_states = getFeasibleState.(ub_samples)

    return (lb, ub, lb_samples, ub_samples, ub_states)
end

"""Select the samples which represent the shortest lower bound tour

    Returns:
        indexes of the samples (vector 1 x n)
"""
function GDIPSolver_find_shortest_tour(solver, fce; dist = nothing, selected = nothing, optimal = false)
    
    if dist === nothing || selected === nothing
        @timeit to "Compute dst" distances = compute_distances(solver.sampling.samples, fce, turning_radius = solver.turning_radius)
    else
        @timeit to "Update" distances = update_distances(solver.sampling.samples, dist, selected, fce, turning_radius = solver.turning_radius)
    end
    stable = false
    if OPTIMIZED_DTRP && !optimal
        for len in [1,2,4,8]
            if !stable
                @timeit to "Find-$len" selected_samples, stable = find_shortest_tour_simplified(distances, len)
            end
        end
    else
        @timeit to "Find-X" selected_samples = find_shortest_tour(distances)
    end
    return selected_samples, distances
end

##################################################
# Main loop over selected scenarios
##################################################
function refine_dtrp(solver::GDIPSolver, max_resolution, in_time_limit)
    while true
        refined = true
        selected_samples = nothing
        @timeit to "Iteration" begin
            while refined
                !isnothing(in_time_limit) && !in_time_limit() && return nothing
                @timeit to "Tour" selected_samples, solver.lb_dist = GDIPSolver_find_shortest_tour(
                    solver, lowerPathGDIP, dist = solver.lb_dist, selected = selected_samples)
                !isnothing(in_time_limit) && !in_time_limit() && return nothing
                @timeit to "Refine" refined = refine_samples(solver.sampling, selected_samples, solver.resolution)
            end
        end        
        solver.resolution *= 2
        if solver.resolution > max_resolution
            break
        end
    end
    !isnothing(in_time_limit) && !in_time_limit() && return nothing
    @timeit to "Final" (lb, ub, lb_samples, ub_samples, ub_states) = GDIPSolver_return_bounds(solver; lb_dist = solver.lb_dist)
    return (lb, ub, lb_samples, ub_samples, ub_states, solver.sampling)
end
