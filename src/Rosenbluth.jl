module Rosenbluth

# Write your package code here.
export RosenbluthSampleable, GARMSampleable, rosenbluth, garm, perm, pegarm, sample, growshrinkgarm, flatgarm, garmsample, PruneEnrichMethod
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

function shrink!(model::GARMSampleable)
    throw(ArgumentError("shrink! not implemented for $(typeof(model))"))
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


function sample(::Type{T}, max_size::Int, num_samples::Int; prune_enrich=false) where {T<:GARMSampleable}
    if prune_enrich
        return pegarm(T, max_size, num_samples)
    else
        return garm(T, max_size, num_samples)
    end
end

function garm(::Type{T}, max_size::Int, num_samples::Int; logging=true) where {T<:GARMSampleable}
    @debug "garm called"
    weights = zeros(Float64, max_size)
    samples = zeros(Int, max_size)

    for t in 1:num_samples
        if logging && t % (num_samples ÷ 20) == 0
            println("Tour: ", t)
        end

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

function pegarm(::Type{T}, max_size::Int, num_tours::Int; logging=true) where {T<:GARMSampleable}
    @debug "pegarm called"
    return flatgarm(T, max_size, num_tours, (max_size,), @inline (model::T) -> (size(model),); logging=logging)
end

function flattour!(::Type{T}, max_size::Int, weights, samples, started_tours, bin_function::Function) where {T<:GARMSampleable}
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


function flatgarm(::Type{T}, max_size::Int, num_tours::Int, results_dimensions::Tuple, bin_function::Function; logging=true) where {T<:GARMSampleable}
    weights = zeros(Float64, results_dimensions)
    samples = zeros(Int, results_dimensions)

    for t in 1:num_tours
        if logging && t % (num_tours ÷ 20) == 0
            println("Tour: ", t)
        end
        flattour!(T, max_size, weights, samples, t, bin_function)
    end

    return weights ./ num_tours, samples
end

function atmosphericflattening(::Type{T}, max_size::Int, num_tours::Int; logging=true) where {T<:GARMSampleable}
    results_dimensions = (max_size, max_aplus(T, max_size), max_aminus(T, max_size))

    return flatgarm(T, max_size, num_tours, results_dimensions, @inline (model::T) -> (size(model), positive_atmosphere(model), negative_atmosphere(model)); logging=logging)
end

function flatgrowshrinktour!(::Type{T}, max_size::Int, weights, samples, normalize_function::Function, bin_function::Function) where {T<:GARMSampleable}
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
            bin_index = bin_function(model)

            weight[n] = (n > 1 ? weight[n-1] : 1)
            weight[n] *= prev_aplus
            weight[n] /= negative_atmosphere(model)

            weights[bin_index...] += weight[n]
            samples[bin_index...] += 1
        end

        if (n >= max_size || prev_aplus == 0)
            copies[n] = 0
        else
            ratio = weight[n] * normalize_function(n) / weights[bin_index...]
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

function flatgrowshrinktour!(::Type{T}, max_size::Int, weights, samples, started_tours::Int, bin_function::Function) where {T<:GARMSampleable}
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
            bin_index = bin_function(model)

            weight[n] = (n > 1 ? weight[n-1] : 1)
            weight[n] *= prev_aplus
            weight[n] /= negative_atmosphere(model)

            weights[bin_index...] += weight[n]
            samples[bin_index...] += 1
        end

        if (n >= max_size || prev_aplus == 0)
            copies[n] = 0
        else
            ratio = weight[n] * started_tours / weights[bin_index...]
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

function growshrinkstacktour!(::Type{T}, max_size::Int, weights, samples, started_tours, bin_function) where {T}
    model = T()
    enrichment_stack = Vector{Tuple{Int, Float64}}()
    push!(enrichment_stack, (size(model), 1.0))

    bin_index = bin_function(model)

    while !isempty(enrichment_stack)
        n, current_weight = pop!(enrichment_stack)
        while size(model) > n
            shrink!(model)
        end

        prev_aplus = positive_atmosphere(model)

        if prev_aplus > 0
            grow!(model)
            n = size(model)
            bin_index = bin_function(model)

            current_weight *= prev_aplus
            current_weight /= negative_atmosphere(model)

            weights[bin_index...] += current_weight
            samples[bin_index...] += 1
        end

        if n >= max_size || prev_aplus == 0
            continue
        else
            ratio = current_weight * started_tours / weights[bin_index...]
            copies = floor(Int, ratio) + Int(rand() < (ratio % 1))

            append!(enrichment_stack, [(n, current_weight / ratio) for _ in 1:copies])
        end
    end
end

function growshrinkgarm(::Type{T}, max_size::Int, num_tours::Int; logging=true) where {T<:GARMSampleable}
    @debug "growshrinkgarm called"
    return growshrinkflatgarm(T, max_size, num_tours, (max_size,), @inline (model::T) -> (size(model),); logging=logging)
end

function growshrinkflatgarm(::Type{T}, max_size::Int, num_tours::Int, results_dimensions::Tuple, bin_function::Function; logging=true) where {T<:GARMSampleable}
    weights = zeros(Float64, results_dimensions)
    samples = zeros(Int, results_dimensions)

    for t in 1:num_tours
        if logging && t % (num_tours ÷ 20) == 0
            println("Tour: ", t)
        end
        flatgrowshrinktour!(T, max_size, weights, samples, t, bin_function)
    end

    return weights ./ num_tours, samples
end

function growshrinkatmosphericflattening(::Type{T}, max_size::Int, num_tours::Int; logging=true) where {T<:GARMSampleable}
    results_dimensions = (max_size, max_aplus(T, max_size), max_aminus(T, max_size))

    return growshrinkflatgarm(T, max_size, num_tours, results_dimensions, @inline (model::T) -> (size(model), positive_atmosphere(model), negative_atmosphere(model)); logging=logging)
end

function isspecialized(f::Function, T::Type)
    return any([method.sig <: Tuple{Any,T,Vararg} for method in methods(f)])
end

function isshrinkable(T::Type)
    return isspecialized(shrink!, T)
end

const PruneEnrichMethod = (NONE=:none, STANDARD=:standard, ATMOSPHERIC_FLATTENING=:atm_flat)

function get_sampler(T::Type, prune_enrich::Symbol)
    if prune_enrich == PruneEnrichMethod.NONE
        return garm
    elseif prune_enrich == PruneEnrichMethod.STANDARD
        return isshrinkable(T) ? growshrinkgarm : pegarm
    elseif prune_enrich == PruneEnrichMethod.ATMOSPHERIC_FLATTENING
        return isshrinkable(T) ? growshrinkatmosphericflattening : atmosphericflattening
    end
end

function garmsample(::Type{T}, max_size::Int, num_tours::Int; prune_enrich_method=:none, logging=true) where {T<:GARMSampleable}
    return get_sampler(T, prune_enrich_method)(T, max_size, num_tours; logging=logging)
end

include("Models.jl")

end
# Ideas:
# - disentangle results data from pruning/enrichment target data
#   - this could be done with a function that takes some required input and returns the target weight
#   - the problem may be determining what the required input should be to cover all cases
#   - the cases we may need to cover so far are:
#     - atmospheric flattening
#
# Use Cases:
# - Population control
#   - ratio from binned targets, then scale by average ratio across all bins so far (i.e. typical configurations make 1 copy)
# - Atmospheric flattening