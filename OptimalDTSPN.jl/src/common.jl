##################################################
# Imports and includes, constants and globals
##################################################

import Base.+
import Base.-
import Base.*
import Base./
import Base./

const RES = 8
const MAX_RES = 2^20

##################################################
# Structures
##################################################
struct Point
    coords::Vector{Float64}
    theta::Float64
end

Point() = Point([], 0.0)
Point(x::T, y::T; theta::T = 0.0) where T<:AbstractFloat = Point([x, y], theta)
Point(x::T, y::T, z::T; theta::T = 0.0) where T<:AbstractFloat = Point([x, y, z], theta)
Point(x::Array{T}; theta::T = 0.0) where T<:AbstractFloat = Point(copy(x), theta)
Point(p::Point) = Point(p.coords, p.theta)

PointFromDubins(x::Array{Float64}) = Point(x[1:2], x[3])
toArray(p::Point) = [p.coords[1], p.coords[2], p.theta]

@inline (+)(a::Point, b::Point) = Point(a.coords .+ b.coords)
@inline (-)(a::Point, b::Point) = Point(a.coords .- b.coords)
@inline (*)(a::Point, b::Real) = Point(a.coords .* b)
@inline (*)(a::Real, b::Point) = (*)(b, a)
@inline (/)(a::Point, b::Real) = Point(a.coords ./ b)

@inline (pointdot)(a::Point, b::Point) = dot(a.coords, b.coords)
@inline dist(x::Point) = sqrt(pointdot(x, x))
@inline dist(x::Array{Float64}) = sqrt(dot(x, x))
@inline dist(a::T, b::T) where T = dist(a-b) 

abstract type AbstractTarget end

##################################################
# Structures - DTSPN/ETSPN
##################################################

# Target region given by the location and its sensing radius
mutable struct TargetRegion <: AbstractTarget 
    center::Point 
    radius::Float64
end  
TargetRegion(x::T, y::T, z::T) where T<:AbstractFloat = TargetRegion(Point(x, y), z)


##################################################
# Help functions
##################################################

""" Load problem from file and return `Array{TargetRegion}` """
load(filename, sensing_radius) = map(readlines(open(filename))) do line
    # @show parse.(Float64, split(line))
    sl = split(line)
    if length(sl) == 2
        TargetRegion(parse.(Float64, split(line))..., sensing_radius)
    elseif length(sl) == 3
        TargetRegion(parse.(Float64, split(line))...)
    else
        @assert false
    end
end

function tour_length(points::Vector{Point}, turning_radius::Float64)
    n = length(points)
    sum(Dubins(points[a], points[mod1(a+1, n)], turning_radius) for a=1:n)
end

@fastmath function tour_length(tour; closed=true)
    total_dist = 0.0
    for i in 1:length(tour) - 1
        total_dist += dist(tour[i], tour[i + 1])
    end
    if closed 
        total_dist += dist(tour[1], tour[length(tour)])
    end
    return total_dist
end

# Distance matrix
function distance_matrix(targets::Vector{TargetRegion})
    n = length(targets)
    return [i == j ? 0.0 : dist(targets[i].center, targets[j].center) for i=1:n, j=1:n] 
end

@inline function Dubins(p1::Point, p2::Point, turning_radius::Float64)
    gdip_init_dubins_maneuver(toArray(p1), toArray(p2), turning_radius)
    gdip_get_length()
end

include("dtrpsolver.jl")

function setTimer(timeroutput::TimerOutput)
    global to = timeroutput
end 
