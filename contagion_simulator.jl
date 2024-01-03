using SparseArrays
using StatsBase
using Graphs
using IterTools
using GraphPlot, Compose
using Plots

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


########################################################################
#Set up 
#Graphs
N = 50
social_edge_p = 0.1
G_g = Graphs.grid([N, N])
G_s = social_graph(N^2, social_edge_p)
#[neighbors(G_s, n)|>length for n in 1:nv(G_s)] |> mean #average number of neighbors

#Initial infections
N_initial =  100 #number of initial exposed
initial_nodes = sample(1:nv(G_g), N_initial, replace = false) #sample nodes to expose 
#Percolation parameters
p = 0.6 #percolation probability 
θ = 10 #activation threshold

#Test 
spatial_social_contagion(G_g, G_s, initial_nodes, p, θ)
#out = zeros(1000)
#[out[i] = spatial_social_contagion(G_g, G_s, sample(1:nv(G_g), N_initial, replace = false), p, θ) for i in 1:1000]
#histogram(out, label = false)
#annotate!(0.2, 150, Plots.text("Distribution of infection rates over random initial conditions", :black, :left, 8))
#png("infection-distribution.png")


#Simulate parameters ----------------
#Removing an increasing number of source nodes OR social ties
n_removed = 90
control_state = zeros(n_removed + 1)
source_removal_state = zeros(n_removed + 1)
social_removal_state = zeros(n_removed +1)

#Replicates
reps = 50
control = zeros(reps)
source = zeros(reps)
social = zeros(reps)

#parameters: adjust percolation probability, edge density (to account for the size of network)
for n in 0:n_removed

    #replicates
    for r in 1:reps

        #Draw random nodes to infect (and receive treatment)
        initial_nodes = sample(1:nv(G_g), N_initial, replace = false) #sample nodes to expose 
        removed_nodes = sample(initial_nodes, n, replace = false) #random select nodes to remove
        remaining_nodes = filter(x -> !(x in removed_nodes), initial_nodes)
        social_g = remove_all_edges_of_vertices(G_s, removed_nodes) #remove social ties of the nodes

        #Compute infection rates
        control[r] = spatial_social_contagion(G_g, G_s, initial_nodes, p, θ)
        source[r] = spatial_social_contagion(G_g, G_s, remaining_nodes, p, θ)
        social[r] =  spatial_social_contagion(G_g, social_g, initial_nodes, p, θ)

    end

    #Summarise replicates    
    control_state[n+1] = control |> mean
    source_removal_state[n+1] = source |> mean
    social_removal_state[n+1] = social |> mean

end

#Plot results ----
xs = collect(0:1:n_removed)
total_nodes = N^2
plot(xs, control_state, label = "Control", xlabel = "Number of nodes removed", ylabel = "Proportion of infected")
plot!(source_removal_state, label = "Source removed")
plot!(social_removal_state, label = "Social ties removed")
annotate!(1, 0.30, Plots.text("Initial infection: $N_initial", :black, :left, 8))
annotate!(1, 0.27, Plots.text("Number of nodes: $total_nodes", :black, :left, 8))
annotate!(1, 0.24, Plots.text("θ: $θ", :black, :left, 8))
annotate!(1, 0.21, Plots.text("Percolation p: $p", :black, :left, 8))
annotate!(1, 0.18, Plots.text("Social edge p: $social_edge_p", :black, :left, 8))
#png("social-contagion.png")

