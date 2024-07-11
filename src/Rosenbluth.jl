module Rosenbluth

# Write your package code here.
export RosenbluthSampleable, GARMSampleable, rosenbluth, garm, perm, pegarm, sample
abstract type GARMSampleable end
abstract type RosenbluthSampleable <: GARMSampleable end

function atmosphere(model::RosenbluthSampleable)::Int
    throw(ArgumentError("atmosphere not implemented for $(typeof(model))"))
end

function positive_atmosphere(model::RosenbluthSampleable)::Int
    return atmosphere(model)
end
function negative_atmosphere(model::RosenbluthSampleable)::Int
    return 1
end


function positive_atmosphere(model::GARMSampleable)::Int
    throw(ArgumentError("positive_atmosphere not implemented for $(typeof(model))"))
end
function negative_atmosphere(model::GARMSampleable)::Int
    throw(ArgumentError("negative_atmosphere not implemented for $(typeof(model))"))
end

function grow!(model::RosenbluthSampleable)
    throw(ArgumentError("grow! not implemented for $(typeof(model))"))
end

function size(model::GARMSampleable)
    throw(ArgumentError("size not implemented for $(typeof(model))"))
end

function sample(::Type{T}, max_size::Int, num_samples::Int; prune_enrich=false) where {T<:GARMSampleable}
    if T <: RosenbluthSampleable
        if prune_enrich
            return perm(T, max_size, num_samples)
        else
            return rosenbluth(T, max_size, num_samples)
        end
    else
        if prune_enrich
            return pegarm(T, max_size, num_samples)
        else
            return garm(T, max_size, num_samples)
        end
    end
end

function rosenbluth(::Type{T}, max_size::Int, num_samples::Int) where {T<:RosenbluthSampleable}

    weights = zeros(Float64, max_size)

    samples = zeros(Int, max_size)

    for _ in 1:num_samples
        model = T()
        weight = 1.0
        n = size(model)
        while n < max_size
            weight *= atmosphere(model)
            if weight == 0
                break
            end
            grow!(model::T)
            n = size(model)

            weights[n] += weight
            samples[n] += 1
        end
    end

    return weights ./ num_samples, samples
end

function garm(::Type{T}, max_size::Int, num_samples::Int) where {T<:GARMSampleable}

    weights = zeros(Float64, max_size)

    samples = zeros(Int, max_size)

    for _ in 1:num_samples
        model = T()
        weight = 1.0

        n = size(model)
        while n < max_size

            weight *= positive_atmosphere(model)
            if weight == 0
                break
            end
            grow!(model::T)

            weight /= negative_atmosphere(model)

            n = size(model)

            weights[n] += weight
            samples[n] += 1

        end
    end

    return weights ./ num_samples, samples
end

function perm(::Type{T}, max_size::Int, num_tours::Int) where {T<:RosenbluthSampleable}
    weights = zeros(Float64, max_size)
    samples = zeros(Int, max_size)

    for _ in 1:num_tours
        enrichment_stack = Vector{Tuple{T,Float64}}()
        push!(enrichment_stack, (T(), 1.0))

        while !isempty(enrichment_stack)
            model, weight = pop!(enrichment_stack)

            weight *= atmosphere(model)
            if weight == 0
                continue
            end
            grow!(model)

            n = size(model)
            weights[n] += weight
            samples[n] += 1

            if n >= max_size
                continue
            end

            ratio = weight * samples[1] / weights[n]
            p = ratio % 1
            copies = floor(Int, ratio)
            if rand() < p
                copies += 1
            end

            for _ in 1:copies
                push!(enrichment_stack, (deepcopy(model), weight / ratio))
            end
        end
    end

    return weights ./ num_tours, samples
end

function pegarm(::Type{T}, max_size::Int, num_tours::Int) where {T<:GARMSampleable}
    weights = zeros(Float64, max_size)
    samples = zeros(Int, max_size)

    for _ in 1:num_tours
        enrichment_stack = Vector{Tuple{T,Float64}}()
        push!(enrichment_stack, (T(), 1.0))

        while !isempty(enrichment_stack)
            model, weight = pop!(enrichment_stack)

            weight *= positive_atmosphere(model)
            if weight == 0
                continue
            end
            grow!(model)
            weight /= negative_atmosphere(model)

            n = size(model)
            weights[n] += weight
            samples[n] += 1

            if n >= max_size
                continue
            end

            ratio = weight * samples[1] / weights[n]
            p = ratio % 1
            copies = floor(Int, ratio)
            if rand() < p
                copies += 1
            end

            for _ in 1:copies
                push!(enrichment_stack, (deepcopy(model), weight / ratio))
            end
        end
    end

    return weights ./ num_tours, samples
end

end