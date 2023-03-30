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

            if iNodeBTrue==jNodeBTrue
                if iNodeBTrue==true
                    reverse!(iRoute.nodes)
                else
                    reverse!(jRoute.nodes)
                end
            else
                if iNodeBTrue==true
                    reverse!(iRoute.nodes)
                    reverse!(jRoute.nodes)
                end
            end
            if length(iRoute.nodes) > 1
               iNode.isInterior = true
            end
            if length(jRoute.nodes) > 1
               jNode.isInterior = true
            end
            iRoute.nodes = vcat(iRoute.nodes, jRoute.nodes)
            jRoute.nodes= []
            iRoute.demand += jRoute.demand
            iRoute.cost += jRoute.cost - DCostNodes[iNode.ID] - DCostNodes[jNode.ID]

            #         # iRoute will contain either edge (depot, i) or edge (i, depot)
            #         iEdge = getDepotEdge(iRoute, iNode) # iEdge is either (0,i) or (i,0)
            #         # remove iEdge from iRoute and update iRoute cost
            #         iRoute.edges.remove(iEdge)
            #         iRoute.cost -= iEdge.cost
            #         # if there are multiple edges in iRoute, then i will be interior
            #         if len(iRoute.edges) > 1: iNode.isInterior = True
            #         # if new iRoute does not start at 0 it must be reversed


            # Change the route of the nodes in jRoute
            for i in 1:length(jRoute.nodes)
                println("i: ", i)
                println("jRoute.nodes[i].ID: ", jRoute.nodes[i].ID)
                println("iRoute.ID: ", iRoute.ID)
                println("DictNodeToRoute[jRoute.nodes[i].ID]: ", DictNodeToRoute[jRoute.nodes[i].ID])
                DictNodeToRoute[jRoute.nodes[i].ID] = iRoute.ID
            end





    end
end

end