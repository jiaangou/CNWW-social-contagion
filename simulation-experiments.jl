include("contagion_simulator.jl")

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
θ = 5 #activation threshold

#Test 
spatial_social_contagion(G_g, G_s, initial_nodes, p, θ)
out = zeros(1000)
[out[i] = spatial_social_contagion(G_g, G_s, sample(1:nv(G_g), N_initial, replace = false), p, θ) for i in 1:1000]
histogram(out, label = false)
#annotate!(0.2, 150, Plots.text("Distribution of infection rates over random initial conditions", :black, :left, 8))
#png("infection-distribution.png")


#######################################
#######################################
#Activation threshold -----------------
thetas = 1:30
t_out = zeros(length(thetas))

for th in thetas
    #t_out[th] = [spatial_social_contagion(G_g, G_s, sample(1:nv(G_g), N_initial, replace = false), p, th) for r in 1:50] |> mean
    t_out[th] = [spatial_social_contagion(G_g, G_s, sample(1:nv(G_g), 50, replace = false), p, th) for r in 1:50] |> mean

end

plot(collect(thetas), t_out, label = false, ylabel = "Proportion of infected", xlabel = "Activation threshold (θ)")
#######################################
#######################################


#######################################
#######################################
#Percolation probability  -----------------
perc_p = 0.1:0.05:0.8
p_out = zeros(length(perc_p))
#p_out2 = zeros(length(perc_p))

for p in 1:length(perc_p)
    p_out[p] = [spatial_social_contagion(G_g, G_s, sample(1:nv(G_g), N_initial, replace = false), perc_p[p], θ) for r in 1:50] |> mean
end

#for p in 1:length(perc_p)
#    p_out2[p] = [spatial_social_contagion(G_g, G_s, sample(1:nv(G_g), 50, replace = false), perc_p[p], θ) for r in 1:50] |> mean
#end

plot(collect(perc_p), p_out, label = false, ylabel = "Proportion of infected", xlabel = "Percolation probability (p)", labels ="Initial nodes: 100" )
plot!(collect(perc_p), p_out2, label = false, color = "red", labels = "Initial nodes: 50")
#png("perco-p.png")
#######################################
#######################################


#######################################
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


    #Plot


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

#######################################
#######################################
#Treatments: theta, initial infections, percolation p, edge density 
initial = 50:5:100 |> collect
percolation_p = 0.2:0.1:0.8 |> collect
edge_p = 0.01:0.02:0.2 |> collect
thetas = 5:3:25 |> collect

treatment = (Iterators.product(initial, percolation_p, edge_p, thetas)|> collect)[:]


theta = 5:1:25
treatment = collect(Iterators.product(1:10, 1:10, 1:10))
treatment[:] |> length
treatment[2][1]


x_vals = [t[1] for t in treatment]
y_vals = [t[2] for t in treatment]
scatter(x_vals, y_vals, label="Grid points", color = "black", labels = false)



