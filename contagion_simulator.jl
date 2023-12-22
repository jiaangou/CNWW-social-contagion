#Generate NxN grid
#using LightGraphs
using SparseArrays
using StatsBase
using Graphs
using IterTools
using GraphPlot, Compose
#using Plots


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
function infected_neighbors(social_graph, infected)

    v_n = nv(social_graph)
    
    #Set up variable for infected neighbors of each node
    infected_neighbors = zeros(v_n)

    #Iterate over vertices and count its the number of infected neighbors
    for v in 1:v_n
        ns = neighbors(social_graph, v) #get neighbors
        infected_neighbors[v] = sum(infected[ns]) #count number of infected neighbors

    end

    return infected_neighbors
end


#4. Social contagion
function social_contagion(social_graph, exposed, infected, θ)

    #Set up counters and variables 
    #I = I0  #Assign initial condition to infection variable
    new_infections = sum(infected) #these are the number of new infections at current step (loop stops when this is 0)

    while new_infections > 0

       #Spread infection
       #1. compute number of neighbors that are infected
       I_neighbors = infected_neighbors(social_graph, infected) #infected neighbors of each vertex

       #2. Check which vertices are exposed and exceed threshold  
       exposed_and_exceed = exposed .&& (I_neighbors .>= θ) #vertices that are exposed (E) and have more than θ infected neighbors (I_n)

       #3. Check if these vertices are already infected
       already_infected = sum(exposed_and_exceed .&& infected)  #of those exposed and exceeded, how many are already infected
       new_infections = sum(exposed_and_exceed) - already_infected

       #Update infections
       if new_infections > 0
           infected[exposed_and_exceed] .= 1 #updated infections
       end

       #print to check
       #println("New infection: $new_infections")
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
    I0 = zeros(n_v) #infection vector
    I0[n0] .= 1   #initial infected nodes
    I_neighbors = infected_neighbors(SocialGraph, I0) #number of infected neighbors
    
    #Percolate supply in space graph
    E = percolate_from_locations(SpaceGraph, n0, p) #Exposure state of vertices
    
    #Percolate infection in social graph
    I = social_contagion(SocialGraph, E, I0[:], θ)

    #Summarise results
    prop_infected = sum(I)/length(I)

    return prop_infected
end




########################################################################
#Set up 

#Graphs
N = 10
G_g = grid([N, N]) 
G_s = social_graph(N^2, 0.1)

#Initial infections
N_initial = 10 #number of initial exposed
initial_nodes = sample(1:nv(G_g), N_initial, replace = false) #sample nodes to expose 

#Percolation parameters
p = 0.6 #percolation probability 
θ = 3 #activation threshold

#Test
spatial_social_contagion(G_g, G_s, initial_nodes, p, θ)

