#Debugging

#Iterate vertices to get infection status
#neighbor_v = neighbors(space, 2)  #neighbors of vertex 1
#neighbor_infections = infected[neighbor_v] #Infection status of neighbors
#neighbor_v[findall(x -> x == 1, neighbor_infections)] #vertex ID of infected neighbors

#Social contagion
function social_contagion(social_graph, exposed, infected, θ)

    #Set up counters and variables 
    new_infections = sum(infected) #these are the number of new infections at current step (loop stops when this is 0)

    while new_infections > 0

       #Spread infection
       #1. compute number of neighbors that are infected
       I_neighbors = infected_neighbors(social_graph, infected) #infected neighbors of each vertex

       #2. Check which vertices are exposed and exceed threshold  
       exposed_and_exceed = exposed .&& (I_neighbors .>= θ) #vertices that are exposed (E) and have more than θ infected neighbors (I_n)
       
       #3. Check if there are new infections and update 
       new_infections = sum(exposed_and_exceed .| infected) - sum(infected)
       infected = exposed_and_exceed .| infected #update infection status

       #println("New infections: $new_infections")

   end

   return infected
end

#-------------------------
#Tests
#-------------------------
#Set up constants ---
#Graphs
s = 10
space = grid([s, s])
social = social_graph(s^2, 0.1)
p = 0.6
θ = 1

#Initial conditions
num_initial_infections = 10
initial_n = sample(1:nv(space), num_initial_infections, replace = false)
E = percolate_from_locations(space, initial_n, p)
infected = zeros(nv(space)) |> BitVector #infection status should be binary 
infected[initial_n] .= 1


#Vary parameter values ----
using Plots
#Vary θ
θ = collect(0:5)
social_infection = zeros(length(θ))

#Social spread
[social_infection[i] = social_contagion2(social, E, infected[:], θ[i]) |> sum for i in 1:length(θ)]

plot(θ,social_infection)


#Social + Spatial spread
result = []
reps = 10
for th in θ
    #[spatial_social_contagion(space, social, initial_n, p, th) for r in 1:reps] |> mean
    result = append!(result, [spatial_social_contagion(space, social, initial_n, p, th)  for r in 1:reps] |> mean)
end

plot(θ,result)

#Vary social ties
num_nodes_removed = 2
rand_node = sample(1:nv(social), num_nodes_removed, replace = false)
social_removed = remove_social_ties(social, rand_node)

neighbors(social, 47)
neighbors(social_removed, 47)




#DEBUG: Social tie removal
social_test = copy(social)
for v in rand_node
    neighbors_to_remove = neighbors(social_test, v)
    num_edges = length(neighbors_to_remove)
    #Number of neighbors before removal
    println("BEFORE removal: vertex $v has $num_edges num_edges")
    
    for n in neighbors_to_remove
        rem_edge!(social_test, v, n)
    end
    num_edges = neighbors(social_test, v) |> length
    println("AFTER removal: vertex $v has $num_edges num_edges")

end

social_test = copy(social)
neighbor_vertices = neighbors(social_test, 47)
#rem_edge!(social_test, 47, neighbor_vertices)

for n in 1:length(neighbor_vertices)
    node = neighbor_vertices[n]
    rem_edge!(social_test, 47, node)
    println("Removing node: $node")
end

[rem_edge!(social_test, 47, n) for n in neighbor_vertices]
neighbors(social_test, 47)
rem_edge!(social_test, 47, 56) 

#NOtE: Edge removal works  manually but has issues inside a loop. 