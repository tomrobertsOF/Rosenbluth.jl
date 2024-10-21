
struct SAW2D <: RosenbluthSampleable
    history::Vector{Tuple{Int,Int}}
end
SAW2D() = SAW2D([])
function Rosenbluth.atmosphere(model::SAW2D)::Int
    if isempty(model.history)
        return 1
    end

    directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    endpoint = model.history[end]
    neighbors = [(endpoint[1] + direction[1], endpoint[2] + direction[2]) for direction in directions]
    return sum([neighbor ∉ model.history for neighbor in neighbors])
end
function Rosenbluth.grow!(model::SAW2D)
    if isempty(model.history)
        push!(model.history, (0, 0))
        return
    end

    directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    endpoint = model.history[end]
    neighbors = [(endpoint[1] + direction[1], endpoint[2] + direction[2]) for direction in directions]
    neighbors = [neighbor for neighbor in neighbors if neighbor ∉ model.history]
    if length(neighbors) == 0
        throw(ArgumentError("No neighbors to grow to"))
    end
    neighbor = neighbors[rand(1:length(neighbors))]
    push!(model.history, neighbor)
    return
end
function Rosenbluth.size(model::SAW2D)
    return length(model.history)
end
