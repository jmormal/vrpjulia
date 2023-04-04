#=
vrp_heuristic:
- Julia version: 
- Author: jmormal
- Date: 2023-03-30
=#
include("./vrpobjects.jl")
using Distances
using BenchmarkTools
using Debugger
using Distributions

using ProfileView
using Base.Threads
using Random
using LoopVectorization
function shuffled_list(n::Int64, index1)
    beta=10
    c_list=collect(1:n)
    o_list=Vector{Int64}(undef,n)
    k=rand(1:length(index1))
    k1=0
    for j in 1:length(c_list)
        k+=1
        k1+=1
        index=index1[k%length(index1)+1] % n+1
        o_list[k1]=c_list[index]
        deleteat!(c_list,index)
    end
    return o_list
end
function shuffled_list_knuth(n::Int64, index1)
    beta=10
    c_list=collect(1:n)
    k=0
    k1=rand(1:length(index1))

    for j in 1:n

        k+=1
        k1+=1
        @inbounds index=(index1[k1%length(index1)+1]) % (n-k+1)+k
        @inbounds c_list[k],c_list[index]=c_list[index],c_list[k]
    end
    return c_list
end

#= function knuth_shuffle_geometric(n::Int)
    a = collect(1:n)
    for i in 1:n
        j = i + rand(Geometric(0.2))
        if j > n
            j = n
        end
        a[i], a[j] = a[j], a[i]
    end
    return a
end
 =#
function read_data(filename)

# with open(fileName) as instance:g
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
        aNode = Node(i, data[1], data[2], data[3], false, 0)
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
#     DCostNodes= Dict{Int64,Float64}()
    DCostNodes= Dict{Int64,Float64}(nodes[i].ID => D[1,i] for i = 2:length(nodes))
#     for i in 2:length(nodes)
#
#     end
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
        nodes[k].route=k-1

    end
    solution= Solution(0, routes, cost ,nodes)
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





function CWS(solution::Solution, nodes::Vector{Node}, OrderedEdges, DictNodeToRoute, DCostNodes,S, D, best_sol,index1)
    vehCap= 100

    for k = shuffled_list_knuth(length(OrderedEdges),index1)
#         k+=rand((-5:5))
#         if k<0:
#             k=1
#         end
#         if k>length(OrderedEdges)
#             k=length(OrderedEdges)
#         end
    #     # select the next edge from the list
    #     ijEdge = savingsList.pop(0) # select the next edge from the list
         @inbounds i1=OrderedEdges[k][1]
         @inbounds i2=OrderedEdges[k][2]
         @inbounds iNode = nodes[i1]
         @inbounds jNode = nodes[i2]

        iNR1=iNode.route
        iNR2=jNode.route
    #     # determine the routes associated with each node
        @inbounds iRoute = solution.routes[iNR1]
        @inbounds jRoute = solution.routes[iNR2]
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


            for n in jRoute.nodes
                n.route= iNode.route
#                 DictNodeToRoute[jRoute.nodes[i].ID] = DictNodeToRoute[iRoute.nodes[1].ID]

            end
            for node in jRoute.nodes
                push!(iRoute.nodes, node)    
                # deleteat!(jRoute.nodes,node)
            end

            # push!(iRoute.nodes, [node for node in jRoute.nodes])
            # iRoute.nodes = vcat(iRoute.nodes, jRoute.nodes)
            empty!(jRoute.nodes)
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
    if solution.cost<best_sol
        best_sol=solution.cost

    end
    return best_sol

end

function CWS1(solution::Solution, nodes::Vector{Node}, OrderedEdges, DictNodeToRoute, DCostNodes,S, D)
    vehCap= 100
    k=1
    while length(OrderedEdges)!=0
    #     # select the next edge from the list
    #     ijEdge = savingsList.pop(0) # select the next edge from the list
         @inbounds i1=OrderedEdges[k][1]
         @inbounds i2=OrderedEdges[k][2]
         @inbounds iNode = nodes[i1]
         @inbounds jNode = nodes[i2]

        iNR1=iNode.route
        iNR2=jNode.route
    #     # determine the routes associated with each node
        @inbounds iRoute = solution.routes[iNR1]
        @inbounds jRoute = solution.routes[iNR2]
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


            for n in jRoute.nodes
                n.route= iNode.route
#                 DictNodeToRoute[jRoute.nodes[i].ID] = DictNodeToRoute[iRoute.nodes[1].ID]

            end

            iRoute.nodes = vcat(iRoute.nodes, jRoute.nodes)
            jRoute.nodes= []
            iRoute.demand += jRoute.demand
            @inbounds iRoute.cost += jRoute.cost - DCostNodes[iNode.ID ] - DCostNodes[jNode.ID] + D[iNode.ID, jNode.ID]




        end
        popfirst!(OrderedEdges)

    end
    solution.cost= 0
    for i in 1:length(solution.routes)
        @inbounds if length(solution.routes[i].nodes) > 0
            @inbounds solution.cost += solution.routes[i].cost
        end
    end

end



function main(instanceName, index1)

#     instanceName = "A-n80-k10" # name of the instance

    fileName = "data/" * instanceName
    solution, nodes, OrderedEdges, DictNodeToRoute, DCostNodes,S ,D ,routes, cost = read_data(fileName)
    num_routes=10000*3*4*10
    Edges=Vector{typeof(OrderedEdges)}(undef, 2000000)

    Solutions=Vector{Solution}(undef, num_routes)

    best_sol=deepcopy(Solution(0,  deepcopy(routes), cost,deepcopy(nodes) ))


#  for k in 1:num_routes
#         Solutions[k]=Solution(0,  deepcopy(routes), cost,deepcopy(nodes ))
#
#     end

#     best_sol=Solutions[1]
    best_sol_cost=100000000
    for k in 1:num_routes

        best_sol_cost=CWS(solution, solution.nodes, OrderedEdges,
         DictNodeToRoute, DCostNodes,S,D,best_sol_cost  , index1)

         for node in solution.nodes[2:end]
            node.isInterior=false
            node.route=node.ID-1
            empty!(solution.routes[node.route].nodes)
            push!(solution.routes[node.route].nodes,node)
            solution.routes[node.route].demand=node.demand
            solution.routes[node.route].cost=2*D[1,node.ID]
            DictNodeToRoute[node.ID] = node.ID-1
            solution.cost=0

        end
#
   end
   println(best_sol_cost)
#     best_sols=[Solutions[l] for l in 1:Threads.nthreads()]
#     thread_count=[0 for l in 1:Threads.nthreads()]
#     @threads for k in 1:num_routes
#
#     best_sols[Threads.threadid()]=CWS(Solutions[k], Solutions[k].nodes, OrderedEdges, DictNodeToRoute, DCostNodes,S,D,best_sols[Threads.threadid()]
#     )
#     thread_count[Threads.threadid()]+=1
#     end

#    println([b.cost for b in best_sols])
#    println(thread_count)
#
#   best_sol_costs=[100000000 for _ in 1:Threads.nthreads()]
#     @threads for k in 1:num_routes
#
#         best_sol_costs[Threads.threadid()]=CWS(solution, solution.nodes, OrderedEdges,
#          DictNodeToRoute, DCostNodes,S,D,best_sol_cost  )
#
#          for node in solution.nodes[2:end]
#             node.isInterior=false
#             node.route=node.ID-1
#             solution.routes[node.route].nodes=Node[node]
#             solution.routes[node.route].demand=node.demand
#             solution.routes[node.route].cost=2*D[1,node.ID]
#             DictNodeToRoute[node.ID] = node.ID-1
#             solution.cost=0
#

#         end
# #
#    end
#    println(best_sol_cost)
end
index1=Tuple(rand(Geometric(0.3),100000))
# println(index1)
for filename in filter(x -> occursin(r"\.txt$", x), readdir("data"))[1:1]
    

    println(filename)
       @time main(filename,index1)

   @profview main(filename,index1)
    readline()
#     #
#     @time main()
#      main()
#     readline()
    # b=@benchmark main($filename,$index1)  time_tolerance=0.01 gctrial=true samples=500 evals=1
    
    # println("Benchmarking vrp_heuristic.jl")
    # println("BenchmarkTools.Trial: ")
    # show(b)
    # # show the minimum time
    # println("BenchmarkTools.Trial minimum time: ")
    # println(minimum(b.times))
    # # show the maximum time
    # println("BenchmarkTools.Trial maximum time: ")
    # println(maximum(b.times), " arg ", argmax(b.times)/length(b.times))
    # # show the mean time
    # println("BenchmarkTools.Trial mean time: ")
    # println(mean(b.times))
    # # show the standard deviation
# #
end


