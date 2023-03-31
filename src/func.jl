#=
func:
- Julia version: 
- Author: jmormal
- Date: 2023-03-30
=#
include("./vrpobjects.jl")
println("vrp_heuristic.jl")
using Distances
using BenchmarkTools
using Debugger

function CWS(solution::Solution, nodes::Vector{Node}, OrderedEdges, DictNodeToRoute, DCostNodes)
    vehCap= 100

    for k in 1:length(OrderedEdges)

    #     ijEdge = savingsList.pop(0) # select the next edge from the list
         iNode = nodes[OrderedEdges[k][1]]
         jNode = nodes[OrderedEdges[k][2]]

    #     # determine the routes associated with each node
        iRoute = solution.routes[DictNodeToRoute[OrderedEdges[k][1]]]
        jRoute = solution.routes[DictNodeToRoute[OrderedEdges[k][2]]]

    #     # check if merge is possible
        isMergeFeasible = checkMergingConditions(iNode, jNode, iRoute, jRoute,vehCap)
    #     # if all necessary conditions are satisfied, merge
        if isMergeFeasible == true

            # deterimine if iNode is at the beginning or end of iRoute
            iNodeBTrue = iRoute.nodes[1].ID == iNode.ID
            # determine if jNode is at the beginning or end of jRoute
            jNodeBTrue = jRoute.nodes[1].ID == jNode.ID








        end
end

end