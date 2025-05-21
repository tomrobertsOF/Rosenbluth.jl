import Base: iterate, length, +, ==

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
    history::Vector{Point2D}
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
SiteTree() = SiteTree(Vector{Point2D}(), Set{Point2D}(), Set{Point2D}([(0,0)]), Set{Point2D}())
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

    push!(model.history, new_site)

    push!(model.occupied, new_site)
    # Remove from a_plus
    delete!(model.growth_candidates, new_site)
    # Remove newly invalid sites from a_plus. I.e. unoccupied neighbors of new_node
    new_neighbors = neighbors(new_site)
    for neighbor in new_neighbors
        delete!(model.growth_candidates, neighbor)
    end
    # Remove neighbor from a_minus. There will only be one by construction, except for the first step where there will be none
    occupied_neighbors = intersect(model.occupied, new_neighbors)
    if !isempty(occupied_neighbors)
        neighbor = first(occupied_neighbors)
        if neighbor in model.shrink_candidates && (neighbor != (0, 0) || length(getoccupiedneighbors(neighbor, model.occupied)) > 1)
            delete!(model.shrink_candidates, neighbor)
        end
    end
    # Add new site to a_minus
    push!(model.shrink_candidates, new_site)
    # Add new valid neighbors to a_plus
    for neighbor in new_neighbors
        # Valid growth sites are unoccupied and have exactly one occupied neighbor
        if !in(neighbor, model.occupied) && length(getoccupiedneighbors(neighbor, model.occupied)) == 1
            push!(model.growth_candidates, neighbor)
        end
    end
end

function Rosenbluth.shrink!(model::SiteTree)
    removed_site = pop!(model.history)

    delete!(model.occupied, removed_site)
    # Remove from a_minus
    delete!(model.shrink_candidates, removed_site)
    # Add removed site to a_plus
    push!(model.growth_candidates, removed_site)
    
    for neighbor in neighbors(removed_site)
        if neighbor in model.growth_candidates
            delete!(model.growth_candidates, neighbor)
        end

        if length(getoccupiedneighbors(neighbor, model.occupied)) == 1
            if neighbor in model.occupied
                push!(model.shrink_candidates, neighbor)
            else
                push!(model.growth_candidates, neighbor)
            end
        end
    end
end

function Rosenbluth.max_aplus(::Type{SiteTree}, max_size::Int)
    2 * (max_size + 1) + 4
end
function Rosenbluth.max_aminus(::Type{SiteTree}, max_size::Int)
    (max_size + 2) ÷ 4 + (max_size + 1) ÷ 4 + 2
end
