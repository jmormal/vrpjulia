#=
test1:
- Julia version: 
- Author: jj
- Date: 2023-04-02
=#
using Random
using Distributions
using ProfileView

v=Vector{Int64}(undef, 1000000)

println(v[1])
