function same_size(A::AbstractArray,B::AbstractArray)
    return size(A) == size(B)
end

function same_size(A::AbstractArray,B::AbstractArray,C::AbstractArray)
    return size(A) == size(B) == size(C)
end

function same_size(A::AbstractArray,B::AbstractArray,dim)
    return size(A,dim) == size(B,dim)
end

function same_size(A::AbstractArray,B::AbstractArray,C::AbstractArray,dim)
    return size(A,dim) == size(B,dim) == size(C,dim)
end

function eight_components(A::AbstractArray)
    return size(A,1) == 8
end

function eight_components(A::AbstractArray, B::AbstractArray)
    return size(A,1) == 8 && size(B,1) == 8
end

function not_simultaneously_nothing(a, b)
    return a!==nothing || b!==nothing
end