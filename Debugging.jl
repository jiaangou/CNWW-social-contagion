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
[social_infection[i] = social_contagion(social, E, infected[:], θ[i]) |> sum for i in 1:length(θ)]

social_contagion(social, E, infected[:], θ[1])

plot(θ,social_infection, xlabel = "θ", ylabel = "Number of infected", legend = false)


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

#NOTE: Edge removal works  manually but has issues inside a loop. 

function remove_all_edges_of_vertices(graph, vertices)
    # Create a copy of the graph
    graph_copy = copy(graph)

    # Remove edges connected to the specified vertices
    for vertex in vertices
        for edge in collect(edges(graph_copy))
            if src(edge) == vertex || dst(edge) == vertex
                rem_edge!(graph_copy, edge)
            end
        end
    end

    return graph_copy
end

initial_n
social
new_g = remove_all_edges_of_vertices(social, initial_n)
[neighbors(new_g, i) for i in initial_n]
[neighbors(social, i) for i in initial_n]


####################################
####################################
#Debug treatments: social ties & source node
#start with N inifected + source nodes
#remove source or social ties of r nodes (either randomly or all possible removals)

#Graphs
s = 10
space = grid([s, s])
social = social_graph(s^2, 0.1)
p = 0.6
θ = 1

#Initial conditions
num_initial_infections = 20
initial_n = sample(1:nv(space), num_initial_infections, replace = false)
E = percolate_from_locations(space, initial_n, p)
infected = zeros(nv(space)) |> BitVector #infection status should be binary 
infected[initial_n] .= 1


#Treatments 
r_removed = 1
selected_node = sample(initial_n, r_removed, replace = false)
source_removed_nodes = initial_n[initial_n .!= selected_node] #remove selected node/same
ties_removed_graph = remove_all_edges_of_vertices(social, selected_node)

spatial_social_contagion(space, social, initial_n, p, θ) #control
spatial_social_contagion(space, social, source_removed_nodes, p, θ) #source removed
spatial_social_contagion(space, ties_removed_graph, initial_n, p, θ) #social ties removed

#Iterator
max_num_nodes_removed = 10

#Output variables
control = zeros(max_num_nodes_removed + 1)
source_removal = zeros(max_num_nodes_removed + 1)
social_removal = zeros(max_num_nodes_removed + 1)

#Number of replicates
reps = 10

#selected_nodes = sample(initial_n, n, replace = false)
#source_removed_nodes = filter(x -> !(x in selected_nodes), initial_n)#remove selected node/same
#ties_removed_graph = remove_all_edges_of_vertices(social, selected_node)
#control[1] = [spatial_social_contagion(space, social, initial_n, p, θ) for r in 1:reps] |> mean #control
#source_removal[1] = [spatial_social_contagion(space, social, source_removed_nodes, p, θ) for r in 1:reps] |> mean #source removed
#social_removal[1] = [spatial_social_contagion(space, ties_removed_graph, initial_n, p, θ) for r in 1:reps] |> mean #social ties removed


for n in collect(0:1:max_num_nodes_removed)

      #treatments
      selected_nodes = sample(initial_n, n, replace = false)
      source_removed_nodes = filter(x -> !(x in selected_nodes), initial_n)#remove selected node/same
      ties_removed_graph = remove_all_edges_of_vertices(social, selected_node)

      #compute results
      control[n+1] = [spatial_social_contagion(space, social, initial_n, p, θ) for r in 1:reps] |> mean #control
      source_removal[n+1] = [spatial_social_contagion(space, social, source_removed_nodes, p, θ) for r in 1:reps] |> mean #source removed
      social_removal[n+1] = [spatial_social_contagion(space, ties_removed_graph, initial_n, p, θ) for r in 1:reps] |> mean #social ties removed

end


plot(collect(0:1:max_num_nodes_removed), control)
plot!(collect(0:1:max_num_nodes_removed), source_removal)
plot!(collect(0:1:max_num_nodes_removed), social_removal)



#Simulate parameters -----------
#Number of nodes
source_n = 100
result = zeros(source_n)
reps = 10

#Setup response variables
control = zeros(source_n)
source_removal = zeros(source_n)
social_removal = zeros(source_n)

for n in collect(1:source_n)
    #draw random nodes
    initial_nodes = sample(1:nv(G_g), n, replace = false) #sample nodes to expose 

    #treatments
    rand_node = sample(initial_nodes) #randomly select a node to be manipulated
    source_removed = initial_nodes[initial_nodes .!= rand_node] #remove the random node as a source
    G_s_removed = remove_social_ties(G_s, rand_node) #remove the social ties of that random node

    #run multiple replicates and average
    control[n] = [spatial_social_contagion(G_g, G_s, initial_nodes, p, θ) for r in 1:reps] |> mean
    source_removal[n] = [spatial_social_contagion(G_g, G_s, source_removed, p, θ) for r in 1:reps] |> mean
    social_removal[n] = [spatial_social_contagion(G_g, G_s_removed, initial_nodes, p, θ) for r in 1:reps] |> mean

end
