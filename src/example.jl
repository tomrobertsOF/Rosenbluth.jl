using Revise
using Rosenbluth

import Base: iterate, length, +, ==

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


mutable struct BinaryTree <: GARMSampleable
    n::Int
    k::Int
end
BinaryTree() = BinaryTree(0, 0)
function Rosenbluth.positive_atmosphere(model::BinaryTree)
    return model.n + 1
end
function Rosenbluth.negative_atmosphere(model::BinaryTree)
    return model.k
end
function Rosenbluth.size(model::BinaryTree)
    return model.n
end
function Rosenbluth.grow!(model::BinaryTree)
    if rand() > 2 * model.k // (model.n + 1)
        model.k += 1
    end
    model.n += 1
    return
end

Point2D = Tuple{Int, Int}
function ==(a::Point2D, b::Point2D)
    return all(a .== b)
end
+(a::Point2D, b::Point2D) = a .+ b

neighbors(p::Point2D) = [p + (1, 0), p + (0, 1), p + (-1, 0), p + (0, -1)]

function getoccupiedneighbors(p::Point2D, occupied::Set{Point2D})
    return Set{Point2D}(filter(neighbor -> neighbor in occupied, neighbors(p)))
end

struct SiteTree <: GARMSampleable
    occupied::Set{Point2D}
    growth_candidates::Set{Point2D}
    shrink_candidates::Set{Point2D}
end
function Base.iterate(model::SiteTree, state=1)
    if state > 3 || state < 1
        return nothing
    end
    if state == 1
        return model.occupied, 2
    elseif state == 2
        return model.growth_candidates, 3
    elseif state == 3
        return model.shrink_candidates, nothing
    end
end
SiteTree() = SiteTree(Set{Point2D}(), Set{Point2D}([(0,0)]), Set{Point2D}())
function Rosenbluth.positive_atmosphere(model::SiteTree)
    return length(model.growth_candidates)
end
function Rosenbluth.negative_atmosphere(model::SiteTree)
    return length(model.shrink_candidates)
end
function Rosenbluth.size(model::SiteTree)
    return length(model.occupied)
end
function Rosenbluth.grow!(model::SiteTree)
    new_site = rand(model.growth_candidates)

    occupied, growth_candidates, shrink_candidates = model
    push!(occupied, new_site)
    # Remove from a_plus
    delete!(growth_candidates, new_site)
    # Remove newly invalid sites from a_plus. I.e. unoccupied neighbors of new_node
    new_neighbors = neighbors(new_site)
    for neighbor in new_neighbors
        delete!(growth_candidates, neighbor)
    end
    # Remove neighbor from a_minus. There will only be one by construction, except for the first step where there will be none
    occupied_neighbors = intersect(occupied, new_neighbors)
    if !isempty(occupied_neighbors)
        neighbor = first(occupied_neighbors)
        if neighbor in shrink_candidates && (neighbor != (0, 0) || length(getoccupiedneighbors(neighbor, occupied)) > 1)
            delete!(shrink_candidates, neighbor)
        end
    end
    # Add new site to a_minus
    push!(shrink_candidates, new_site)
    # Add new valid neighbors to a_plus
    for neighbor in new_neighbors
        # Valid growth sites are unoccupied and have exactly one occupied neighbor
        if !in(neighbor, occupied) && length(getoccupiedneighbors(neighbor, occupied)) == 1
            push!(growth_candidates, neighbor)
        end
    end
end


w, s = sample(SiteTree, 50, 10000, prune_enrich=true)