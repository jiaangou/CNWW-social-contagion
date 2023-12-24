using SparseArrays
using StatsBase
using Graphs
using IterTools
using GraphPlot, Compose


#Functions --------------------------------

#1. Social graph generation
function social_graph(N_v, p::Float64)
    g = Graph(N_v, 0)

    # Ensure p is between 0 and 1
    @assert 0 <= p <= 1 "Probability p must be between 0 and 1"

    n = nv(g)  # Number of vertices in the graph

    # Iterate over all pairs of vertices
    for i in 1:n
        for j in (i+1):n  # Avoid self-loops and double counting edges
            if rand() < p  # Connect with probability p
                add_edge!(g, i, j)
            end
        end
    end
    return g
end


#2. Supply percolation 
function percolate_from_locations(lattice, initial_nodes, p)
    # Input: lattice, starting locations, and probability
    # Output: vertices that are activated by percolation process

    # Initialize all nodes as unoccupied
    occupied = falses(nv(lattice))

    # Set initial nodes as occupied
    for node in initial_nodes
        occupied[node] = true
    end

    # Perform the percolation process
    for node in vertices(lattice)
        if occupied[node]
            for neighbor in neighbors(lattice, node)
                if rand() < p && !occupied[neighbor]
                    occupied[neighbor] = true
                end
            end
        end
    end

    return occupied
end

#3.Infection status of social neighbors
function infected_neighbors(SocialGraph, infected)

    v_n = nv(SocialGraph)
    
    #Set up variable for infected neighbors of each node
    infected_neighbors = zeros(v_n)

    #Iterate over vertices and count its the number of infected neighbors
    for v in 1:v_n
        ns = neighbors(SocialGraph, v) #get neighbors
        infected_neighbors[v] = sum(infected[ns]) #count number of infected neighbors

    end

    return infected_neighbors
end


#4. Social contagion
function social_contagion(SocialGraph, exposed, infected, θ)

    #Set up counters and variables 
    new_infections = sum(infected) #these are the number of new infections at current step (loop stops when this is 0)

    while new_infections > 0

       #Spread infection
       #1. compute number of neighbors that are infected
       I_neighbors = infected_neighbors(SocialGraph, infected) #infected neighbors of each vertex

       #2. Check which vertices are exposed and exceed threshold  
       exposed_and_exceed = exposed .&& (I_neighbors .>= θ) #vertices that are exposed (E) and have more than θ infected neighbors (I_n)
       
       #3. Check if there are new infections and update 
       new_infections = sum(exposed_and_exceed .| infected) - sum(infected)
       infected = exposed_and_exceed .| infected #update infection status

       #println("New infections: $new_infections")

   end

   return infected
end


#5. Simulator (spatial + social contagion)
function spatial_social_contagion(SpaceGraph, SocialGraph, n0, p, θ)
    
    #n0 - initial nodes that are both infected and exposed
    #p - percolation probability
    #θ - activation threshold
    
    # Ensure SpaceGraph and SocialGraph have the same number of vertices
    @assert  nv(SpaceGraph) == nv(SocialGraph) "Number of vertices should be equal for both graphs"

    #number of vertices
    n_v = nv(SpaceGraph) 

    #Setup initial infection status
    I0 = zeros(n_v) |> BitVector#infection vector
    I0[n0] .= 1   #initial infected nodes
    #I_neighbors = infected_neighbors(SocialGraph, I0) #number of infected neighbors
    
    #Percolate supply in space graph
    E = percolate_from_locations(SpaceGraph, n0, p) #Exposure state of vertices
    
    #Percolate infection in social graph
    I = social_contagion(SocialGraph, E, I0[:], θ)

    #Summarise results
    prop_infected = sum(I)/length(I)

    return prop_infected
end


# Experimental treatments ------
# 1. Removing social ties
function remove_social_ties(graph, vertex)
    g = copy(graph)

    #Remove neighbors of each vertex
    for v in vertex
        neighbors_to_clear = neighbors(g, v)
        
        for neighbor in neighbors_to_clear
            rem_edge!(g, v, neighbor)
        end
        
    end
    
    return g
end

########################################################################
#Set up 

#Graphs
N = 100
G_g = grid([N, N]) 
G_s = social_graph(N^2, 0.1)

#Initial infections
N_initial = 100 #number of initial exposed
initial_nodes = sample(1:nv(G_g), N_initial, replace = false) #sample nodes to expose 

#Percolation parameters
p = 0.6 #percolation probability 
θ = 2 #activation threshold

#Test
spatial_social_contagion(G_g, G_s, initial_nodes, p, θ)

#Logic -------
#Select random node
rand_node = sample(initial_nodes)
source_removed = initial_nodes[initial_nodes .!= rand_node]
G_s_removed = remove_social_ties(G_s, rand_node)

spatial_social_contagion(G_g, G_s, initial_nodes, p, θ) #control
spatial_social_contagion(G_g, G_s, source_removed, p, θ) #source removed
spatial_social_contagion(G_g, G_s_removed, initial_nodes, p, θ) #social ties removed

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

#-----------------------------------------------------
N = 100
G_g = grid([N, N]) 
G_s = social_graph(N^2, 0.1)

#Initial infections
N_initial = 500 #number of initial exposed
initial_nodes = sample(1:nv(G_g), N_initial, replace = false) #sample nodes to expose 


#Removing an increasing number of source nodes OR social ties
n_removed = 200
control_state = zeros(n_removed)
source_removal_state = zeros(n_removed)
social_removal_state = zeros(n_removed)

reps = 10

#parameters: adjust percolation probability, edge density (to account for the size of network)
for n in 1:n_removed

    removed_nodes = sample(initial_nodes, n, replace = false) #random select nodes to remove
    remaining_nodes = filter(x -> !(x in removed_nodes), initial_nodes)
    social_g = remove_social_ties(G_s, removed_nodes) #remove social ties of the nodes
    
    #compute
    control_state[n] = [spatial_social_contagion(G_g, G_s, initial_nodes, p, θ) for r in 1:reps] |> mean
    source_removal_state[n] = [spatial_social_contagion(G_g, G_s, remaining_nodes, p, θ) for r in 1:reps] |> mean
    social_removal_state[n] = [spatial_social_contagion(G_g, social_g, initial_nodes, p, θ) for r in 1:reps] |> mean

end


#--------------------------
#Social ties debugging
edge_p = collect(0.0001:0.0005:0.1)
threshold = collect(2:5)
#results = zeros(llength(edge_p) * length(threshold))
results = []

for i in 1:length(edge_p)
    
    social_g  = social_graph(N^2, edge_p[i])

    for t in threshold
        θ = t
        results = append!(results, [spatial_social_contagion(G_g, social_g, initial_nodes, p, θ) for r in 1:reps] |> mean)
    end

end

#plot(edge_p, results, xlabel = "Edge probability", ylabel = "Proportion of infected", label = "")


# Prepare a plot
#plot(xlabel="Edge probability", ylabel="Infection %")
#findall(x -> x == 2, threshold)
plot(edge_p, results[1:200], label = "Threshold: 2")
plot!(edge_p, results[201:400], label = "Threshold: 3")
plot!(edge_p, results[401:600], label = "Threshold: 4")
plot!(edge_p, results[601:800], label = "Threshold: 5")


#Note: higher threshold is leading to more infections (i.e. IT SHOULD DECREASE. BUG!!!)

#Visualize graph
using Plots
plot(control_state, label = "Control")
plot!(source_removal_state, label = "Source removed")
plot!(social_removal_state, label = "Social ties removed")

#Increase number of initial infected should shift where the transintion point is 



