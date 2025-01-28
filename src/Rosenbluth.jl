module Rosenbluth

import StatsBase

export RosenbluthSampleable, GARMSampleable, sample, PruneEnrichMethod

include("Sampleable.jl")

# Constants
"""
    PruneEnrichMethod

A constant defining the available pruning and enrichment methods:
- `NONE`: No pruning or enrichment.
- `STANDARD`: Standard pruning and enrichment.
- `ATMOSPHERIC_FLATTENING`: Atmospheric flattening.
"""
const PruneEnrichMethod = (NONE=:none, STANDARD=:standard, ATMOSPHERIC_FLATTENING=:atm_flat)

# Public API
"""
    sample(::Type{T}, max_size::Int, num_tours::Int; prune_enrich_method=:none, logging=true) where {T<:GARMSampleable}

Sample a model of type `T` using the GARM algorithm with optional pruning and enrichment methods.

# Arguments
- `T`: The type of the model to be sampled.
- `max_size`: The maximum size of the model.
- `num_tours`: The number of tours to perform.
- `prune_enrich_method`: The pruning and enrichment method to use (`:none`, `:standard`, or `:atm_flat`).
- `logging`: A boolean indicating whether to log progress.

# Returns
A tuple containing the weights and samples.
"""
function sample(::Type{T}, max_size::Int, num_tours::Int; prune_enrich_method=:none, logging=true) where {T<:GARMSampleable}
    isSampleable(T) || throw(ArgumentError("Type $(T) is not sampleable. It must implement the GARMSampleable or RosenbluthSampleable interface."))

    return get_sampler(T, prune_enrich_method)(T, max_size, num_tours; logging=logging)
end

# Internal Functions
function get_sampler(T::Type, prune_enrich::Symbol)
    if prune_enrich == PruneEnrichMethod.NONE
        return garm
    elseif prune_enrich == PruneEnrichMethod.STANDARD
        return isshrinkable(T) ? growshrinkgarm : pegarm
    elseif prune_enrich == PruneEnrichMethod.ATMOSPHERIC_FLATTENING
        return isshrinkable(T) ? growshrinkatmosphericflattening : atmosphericflattening
    end
end

function garm(::Type{T}, max_size::Int, num_samples::Int; logging=true) where {T<:GARMSampleable}
    @debug "garm called"
    weights = zeros(Float64, max_size)
    samples = zeros(Int, max_size)

    for t in 1:num_samples
        if logging && num_samples > 20 && t % (num_samples ÷ 20) == 0
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

mutable struct Particle{T<:GARMSampleable}
    model::T
    weight::Float64
end
function Particle{T}() where {T<:GARMSampleable}
    Particle{T}(T(), 1.0)
end

function breadthfirstGARM(::Type{T}, max_size::Int, num_samples::Int; logging::Bool=true) where {T<:GARMSampleable}
    weights = zeros(Float64, max_size)
    samples = zeros(Int, max_size)

    particles::Vector{Particle{T}} = [Particle{T}() for _ in 1:num_samples]

    for n in 1:max_size
        Threads.@threads for particle in particles
            step!(particle)
        end

        weights[n] = sum(p.weight for p in particles)
        samples[n] = count(p -> p.weight != 0, particles)
    end

    return weights ./ num_samples, samples
end

function step!(particle::Particle{T}) where {T}
    particle.weight *= positive_atmosphere(particle.model)
    if particle.weight == 0
        return
    end
    grow!(particle.model)
    particle.weight /= negative_atmosphere(particle.model)
    return
end

function breadthfirstPEGARM(::Type{T}, max_size::Int, num_samples::Int; logging::Bool=true) where {T<:GARMSampleable}
    weights = zeros(Float64, max_size)
    samples = zeros(Int, max_size)

    particles::Vector{Particle{T}} = [Particle{T}() for _ in 1:num_samples]

    locks = [Threads.SpinLock() for _ in 1:num_samples]

    for n in 1:max_size
        # Step
        Threads.@threads for particle in particles
            step!(particle)
        end

        # Record
        weights[n] = sum(p.weight for p in particles)
        samples[n] = count(p -> p.weight != 0, particles)

        # Resample
        indices = StatsBase.sample(1:length(particles), StatsBase.Weights([p.weight for p in particles]), length(particles))
        should_copy = zeros(Bool, length(particles))
        new_particles = Vector{Particle{T}}(undef, length(particles))

        Threads.@threads for i in 1:length(particles)

            # Thread-safe update array returning the previous value
            old_copy = @lock locks[indices[i]] begin
                prev, should_copy[indices[i]] = should_copy[indices[i]], true
                prev
            end

            # Optimization: only need to copy if we are resampling the particle more than once,
            # otherwise we can just assign the particle
            new_particles[i] = if old_copy
                deepcopy(particles[indices[i]])
            else
                particles[indices[i]]
            end
            new_particles[i].weight = weights[n] / num_samples
        end

        particles = new_particles
    end

    return weights ./ num_samples, samples
end

function pegarm(::Type{T}, max_size::Int, num_tours::Int; logging=true) where {T<:GARMSampleable}
    @debug "pegarm called"
    return flatgarm(T, max_size, num_tours, (max_size,), @inline (model::T) -> (size(model),); logging=logging)
end

function flatgarm(::Type{T}, max_size::Int, num_tours::Int, results_dimensions::Tuple, bin_function::Function; logging=true) where {T<:GARMSampleable}
    weights = zeros(Float64, results_dimensions)
    samples = zeros(Int, results_dimensions)

    for t in 1:num_tours
        if logging && num_tours > 20 && t % (num_tours ÷ 20) == 0
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

# Tour Functions
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

function flatgrowshrinktour!(::Type{T}, max_size::Int, weights, samples, normalize_function::Function, bin_function::Function) where {T<:GARMSampleable}
    model = T()
    weight = zeros(Float64, max_size)
    copies = zeros(Int, max_size)

    outstanding_copies = 1

    while (size(model) != 1 || outstanding_copies > 0)

        n = size(model)
        prev_aplus = positive_atmosphere(model)

        if prev_aplus > 0
            outstanding_copies -= 1
            if n > 0
                copies[n] -= 1
            end

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

            outstanding_copies += copies[n]

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
    flatgrowshrinktour!(T, max_size, weights, samples, (::Int) -> started_tours, bin_function)
end

function growshrinkstacktour!(::Type{T}, max_size::Int, weights, samples, started_tours, bin_function) where {T}
    model = T()
    enrichment_stack = Vector{Tuple{Int,Float64}}()
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