
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

function Rosenbluth.max_aplus(::Type{BinaryTree}, max_size::Int)
    max_size + 1
end
function Rosenbluth.max_aminus(::Type{BinaryTree}, max_size::Int)
    ceil(Int, max_size / 2)
end
