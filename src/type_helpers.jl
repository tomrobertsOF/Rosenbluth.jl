function isspecialized(f::Function, T::Type)
    return any([method.sig <: Tuple{Any,T,Vararg} for method in methods(f)])
end

function isshrinkable(T::Type)
    return isspecialized(shrink!, T)
end
