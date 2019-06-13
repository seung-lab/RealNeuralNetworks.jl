module Mathes 

import Statistics: mean 
import LinearAlgebra: svd, diagm, diag

"""
    pca(data)
perform PCA using SVD
inputs:
    - data: M x N matrix of input data. (M dimensions, N trials)
outputs:
    - PC: each column is a principle component
    - V: M x 1 matrix of variances
"""
function pca(data::Array{T,2}) where T
    X = data .- mean(data, dims=2)
    Y = X' ./ sqrt(T(size(X,2)-1))
    U,S,PC = svd(Y)
    S = diagm(0=>S)
    V = S .* S
    
    # find the least variance vector
    indexList = sortperm(diag(V); rev=true)

    PCs = map(x->PC[:,x], indexList)
    return PCs, diag(V)[indexList]
end 

end # end of module
