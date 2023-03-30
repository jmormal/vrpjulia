#=
vrpobjects:
- Julia version: 
- Author: jmormal
- Date: 2023-03-30
=#


mutable struct Node
    ID::Int16
    x::Float64
    y::Float64
    demand::Float64
    isInterior::Bool
end


mutable struct Route
    ID::Int8
    cost::Float16
    nodes::Vector{Node}
    demand::Float64
end


mutable struct Solution
    ID::Int16
    routes::Vector{Route}
    cost::Float64
end