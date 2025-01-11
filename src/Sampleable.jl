# This file contains the default interface for the GARMSampleable and RosenbluthSampleable types.
# The user can implment these methods for their own models to use the GARMSampleable interface.


"""
    GARMSampleable

An abstract type that serves as a base for types that can be sampled using the GARM (Generalized Atmospheric Rosenbluth Method) algorithm. Types that inherit from `GARMSampleable` should implement the necessary methods to facilitate sampling.

At a minimum, the following methods should be implemented:

- `positive_atmosphere(model::GARMSampleable)::Int`
- `negative_atmosphere(model::GARMSampleable)::Int`
- `grow!(model::GARMSampleable)`
- `size(model::GARMSampleable)`

These methods define the basic interface required for the GARM algorithm to interact with the model.

Additionally, the following methods can be implemented to provide additional functionality:

If shrinking is supported, the following method should be implemented, allowing for more effieicnt sampling by avoiding unnecessary copies:
- `shrink!(model::GARMSampleable)`

These two methods are required to use the atmospheric flattening enrichment strategy:
- `max_aplus(::Type{T}, max_size::Int) where {T<:GARMSampleable}` 
- `max_aminus(::Type{T}, max_size::Int) where {T<:GARMSampleable}`
"""
abstract type GARMSampleable end


"""
    RosenbluthSampleable

An abstract type that inherits from `GARMSampleable` and serves as a base for types that can be sampled using the Rosenbluth method. Types that inherit from `RosenbluthSampleable` should implement the necessary methods to facilitate sampling specific to the Rosenbluth algorithm.

At a minimum, the following methods should be implemented:

- `atmosphere(model::RosenbluthSampleable)::Int`
- `grow!(model::GARMSampleable)`
- `size(model::GARMSampleable)`

These methods define the basic interface required for the Rosenbluth algorithm to interact with the model.
"""
abstract type RosenbluthSampleable <: GARMSampleable end


function atmosphere(model::RosenbluthSampleable)::Int
    throw(ArgumentError("atmosphere not implemented for $(typeof(model))"))
end

function positive_atmosphere(model::RosenbluthSampleable)::Int
    return atmosphere(model)
end
function negative_atmosphere(::RosenbluthSampleable)::Int
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

function max_aplus(::Type{T}) where {T<:GARMSampleable}
    return n -> max_aplus(T, n)
end
function max_aplus(::Type{T}, max_size::Int) where {T<:GARMSampleable}
    throw(ArgumentError("max_aplus not implemented for $(typeof(T))"))
end

function max_aminus(::Type{T}) where {T<:GARMSampleable}
    return n -> max_aminus(T, n)
end
function max_aminus(::Type{T}, max_size::Int) where {T<:GARMSampleable}
    throw(ArgumentError("max_aminus not implemented for $(typeof(T))"))
end