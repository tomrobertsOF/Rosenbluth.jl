module Rosenbluth

# Write your package code here.
export RosenbluthSampleable, GARMSampleable, rosenbluth, garm, perm, pegarm, sample, growshrinkgarm, flatgarm
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

function grow!(model::GARMSampleable)
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

function max_aplus(::Type{T}, max_size::Int) where {T<:GARMSampleable}
    throw(ArgumentError("max_aplus not implemented for $(typeof(T))"))
end
function max_aminus(::Type{T}, max_size::Int) where {T<:GARMSampleable}
    throw(ArgumentError("max_aminus not implemented for $(typeof(T))"))
end

function bin_dimensions(::Type{T}, max_size::Int) where {T<:GARMSampleable}
    throw(ArgumentError("bin_dimensions not implemented for $(typeof(T))"))
end

function flattour!(::Type{T}, max_size::Int, weights, samples, bin_function::Function) where {T<:GARMSampleable}
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
        bin_index = bin_function(model)
        weights[bin_index...] += weight
        samples[bin_index...] += 1

        if (n == 1)
            started_tours += 1
        end

        if n >= max_size
            continue
        end

        ratio = weight * started_tours / weights[bin_index...]
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

function flatgarm(::Type{T}, max_size::Int, num_tours::Int, results_dimensions::Tuple, bin_function::Function) where {T<:GARMSampleable}
    weights = zeros(Float64, results_dimensions)
    samples = zeros(Int, results_dimensions)
    started_tours = 0

    for t in 1:num_tours
        println("Tour: ", t)
        flattour!(T, max_size, weights, samples, bin_function)
    end

    return weights ./ num_tours, samples
end

function atmosphericflattening(::Type{T}, max_size::Int, num_tours::Int) where {T<:GARMSampleable}
    results_dimensions = (max_size, max_aplus(T, max_size), max_aminus(T, max_size))

    return flatgarm(T, max_size, num_tours, results_dimensions, (model::T) -> (size(model), positive_atmosphere(model), negative_atmosphere(model)))
end

function shrink!(model::GARMSampleable)
    throw(ArgumentError("shrink! not implemented for $(typeof(model))"))
end

function growshrinktour!(::Type{T}, max_size::Int, weights, samples) where {T<:GARMSampleable}
    model = T()
    weight = zeros(Float64, max_size)
    copies = zeros(Int, max_size)

    enrichment_stack = Vector{Bool}()
    push!(enrichment_stack, true)

    while (size(model) != 1 || !isempty(enrichment_stack))
        
        n = size(model)
        prev_aplus = positive_atmosphere(model)

        if prev_aplus > 0
            if n > 0
                copies[n] -= 1
            end
            pop!(enrichment_stack)
            grow!(model)
            n = size(model)

            weight[n] = (n > 1 ? weight[n-1] : 1)
            weight[n] *= prev_aplus
            weight[n] /= negative_atmosphere(model)

            weights[n] += weight[n]
            samples[n] += 1
        end

        if (n >= max_size || prev_aplus == 0)
            copies[n] = 0
        else
            ratio = weight[n] * samples[1] / weights[n]
            p = ratio % 1
            copies[n] = floor(Int, ratio)
            if rand() < p
                copies[n] += 1
            end

            for _ in 1:copies[n]
                push!(enrichment_stack, true)
            end

            weight[n] /= ratio
        end
        if (copies[n] == 0)
            while (n > 1 && copies[n] == 0)
                shrink!(model)
                n = size(model)
            end
        end
    end
end

function growshrinkgarm(::Type{T}, max_size::Int, num_tours::Int) where {T<:GARMSampleable}
    weights = zeros(Float64, max_size)
    samples = zeros(Int, max_size)

    for _ in 1:num_tours
        growshrinktour!(T, max_size, weights, samples)
    end

    return weights ./ num_tours, samples
end

include("Models.jl")

end