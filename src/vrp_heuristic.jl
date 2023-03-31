#=
vrp_heuristic:
- Julia version: 
- Author: jmormal
- Date: 2023-03-30
=#
include("./vrpobjects.jl")
println("vrp_heuristic.jl")
using Distances
using BenchmarkTools
using Debugger
function read_data(filename)

# with open(fileName) as instance:
#     i = 0
#     nodes = []
#     for line in instance:
#         # array data with node data: x, y, demand
#         data = [float(x) for x in line.split()]
#         aNode = Node(i, data[0], data[1], data[2])
#         nodes.append(aNode)
#         i += 1
# in julia this is
    nodes = Node[]
    i = 1
    for line in eachline(filename)
        data = [parse(Float32, x) for x in split(line)]
        aNode = Node(i, data[1], data[2], data[3], false)
        push!(nodes, aNode)
        i += 1
    end

    # Now we are going to define the array X, where there are length(nodes) columns. The first row
    # is the x postion and the second row is the y position. This is used to calculate the distance
    # between nodes.

    X= zeros(2, length(nodes))
    for i in 1:length(nodes)
        X[1,i] = nodes[i].x
        X[2,i] = nodes[i].y
    end
    D= pairwise(Euclidean(), X, dims=2)

#     DCostNodes= Dict{Int64,Float64}( nodes[i].ID => D[1,i] for i in 2:length(nodes))
    DCostNodes= Dict{Int64,Float64}()
    for i in 2:length(nodes)
        DCostNodes[nodes[i].ID] = D[1,i]
    end
    # Now we are going to define the array S, the savings of each pair of nodes. The savings are
    # defined as the difference between the distance between the nodes and the distance between the
    # depot and the nodes. The depot is the first node in the list of nodes.
    # We are going to assign a vector of zeros to S. The length is (length(nodes)+1)*(length(nodes))/2

    S= zeros(floor(Int,((length(nodes)-2)*(length(nodes)-1)/2)))
    # We are going to define D1.
    # It is a dictionary from int => (int, int)

    D1 = Dict{Int64, Tuple{Int64, Int64}}()
    k=0
    for i in 2:length(nodes)
        for j in i+1:length(nodes)
            if j <= length(nodes)
                k+=1
                S[k] = D[1,i] + D[1,j] - D[i,j]
                D1[k] = (i,j)
            end
        end
    end
    order= sortperm(S, rev=true)
    OrderedEdges= Vector{Tuple{Int64, Int64}}()

    for i in 1:length(order)
        push!(OrderedEdges, D1[order[i]])
    end
    sort!(S)
    routes= Route[]
    cost= 0
    DictNodeToRoute= Dict{Int64, Int64}()
    for k in 2:length(nodes)
        cost += D[1,k]*2
        push!(routes, Route(k, D[1,k]*2, [nodes[k]], nodes[k].demand))
        DictNodeToRoute[nodes[k].ID] = k-1

    end
    solution= Solution(0, routes, cost )
    return solution, nodes, OrderedEdges, DictNodeToRoute, DCostNodes,S, D ,routes, cost

end


# def checkMergingConditions(iNode, jNode, iRoute, jRoute):
#     # condition 1: iRoute and jRoure are not the same route object
#     if iRoute == jRoute: return False
#     # condition 2: both nodes are exterior nodes in their respective routes
#     if iNode.isInterior == True or jNode.isInterior == True: return False
#     # condition 3: demand after merging can be covered by a single vehicle
#     if vehCap < iRoute.demand + jRoute.demand: return False
#     # else, merging is feasible
#     return True

function checkMergingConditions(iNode::Node, jNode::Node, iRoute::Route, jRoute::Route, vehCap::Int)
    # condition 1: iRoute and jRoure are not the same route object
#     if length(iRoute.nodes) == 0 || length(jRoute.nodes) == 0
#         return false
#     end
    if iRoute == jRoute
        return false
    end
    # condition 2: both nodes are exterior nodes in their respective routes
    if iNode.isInterior == true || jNode.isInterior == true
        return false
    end
    # condition 3: demand after merging can be covered by a single vehicle
    if vehCap < iRoute.demand + jRoute.demand
        return false
    end
    # else, merging is feasible

    return true
end



function CWS(solution::Solution, nodes::Vector{Node}, OrderedEdges, DictNodeToRoute, DCostNodes,S, D)

    vehCap= 100

    for k in 1:length(OrderedEdges)
    #     # select the next edge from the list
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


            for i in 1:length(jRoute.nodes)

                DictNodeToRoute[jRoute.nodes[i].ID] = DictNodeToRoute[iRoute.nodes[1].ID]

            end

            iRoute.nodes = vcat(iRoute.nodes, jRoute.nodes)
            jRoute.nodes= []
            iRoute.demand += jRoute.demand
            iRoute.cost += jRoute.cost - DCostNodes[iNode.ID ] - DCostNodes[jNode.ID] + D[iNode.ID, jNode.ID]





        end
    end
    solution.cost= 0
    for i in 1:length(solution.routes)
        if length(solution.routes[i].nodes) > 0
            solution.cost += solution.routes[i].cost
        end
    end

end



function main()

    instanceName = "A-n80-k10" # name of the instance
#     instanceName = "A-n32-k5" # name of the instance

    fileName = "data/" * instanceName * "_input_nodes.txt"
    solution, nodes, OrderedEdges, DictNodeToRoute, DCostNodes,S ,D ,routes, cost = read_data(fileName)
    for _ in 1:50000
        # Make a copy of the solution
        solution=Solution(0,  routes, cost )

        CWS(solution, nodes, OrderedEdges, DictNodeToRoute, DCostNodes,S,D)
   end

end


@time main()
# b=@benchmark main() seconds=1 time_tolerance=0.01 gctrial=true samples=1000 evals=1
#
# println("Benchmarking vrp_heuristic.jl")
# println("BenchmarkTools.Trial: ")
# show(b)
# # show the minimum time
# println("BenchmarkTools.Trial minimum time: ")
# println(minimum(b.times))
# # show the median time
# println("BenchmarkTools.Trial median time: ")
# println(b.times[round(Int, length(b.times)/2)])
# # show the maximum time
# println("BenchmarkTools.Trial maximum time: ")
# println(maximum(b.times), " arg ", argmax(b.times)/length(b.times))
# # show the mean time
# println("BenchmarkTools.Trial mean time: ")
# println(mean(b.times))
# # show the standard deviation
# #
