#=
test1:
- Julia version: 
- Author: jj
- Date: 2023-04-02
=#
using Random
using Distributions
using ProfileView

function shuffled_list(n::Int64)
    beta=10
    c_list=collect(1:n)
    o_list=Vector{Int64}(undef,n)
    k=0
    index1=rand(Geometric(0.25),length(c_list))

    for j in 1:length(c_list)
        k+=1
        index=index1[k] % length(c_list)+1
        o_list[k]=c_list[index]
        deleteat!(c_list,index)
    end
    return o_list
end
function creating_random()
    for _ in 1:10000
        shuffled_list(80*80)
    end
end
@time creating_random()
@profview creating_random()
readline()

