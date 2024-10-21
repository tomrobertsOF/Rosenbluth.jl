module Models

using ..Rosenbluth

import Base: iterate, length, +, ==

include("models/BinaryTree.jl")
include("models/SiteTree.jl")
include("models/SAW2D.jl")

end