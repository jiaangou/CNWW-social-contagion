#Degree distribution
function degree_distribution(G)
    num_vertices = nv(G)
    
    neighbors = zeros(num_vertices)
    for n in 1:num_vertices
        neighbors[n] = length(all_neighbors(G, n))
    end

    return neighbors

end


