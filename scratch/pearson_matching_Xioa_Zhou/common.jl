import Statistics: mean, std, quantile
import FastGaussQuadrature: gausshermite
import Polynomials: Polynomial
import IntervalRootFinding: roots
import Memoize: @memoize

using Distributions
using IntervalArithmetic


margins = [
    (Beta(2, 3)),
    (Normal(12, 2)),
    (Exponential(3.14))
]


"""
The Probabilists' version of the Hermite polynomials
"""
@memoize function He(x, n::Int)
    if n == 0
        return length(x) > 1 ? ones(length(x)) : 1
    elseif n == 1
        return x
    else
        return x .* He(x, n-1) .- (n-1) .* He(x, n-2)
    end
end


"""
The Physicists' version of the Hermite polynomials
"""
function H(x, n::Int)
    return 2^(n/2) * He(x*√2, n)
end


"""
	z2x(d::Distribution, x::AbstractArray)
Convert samples from a standard normal distribution to a given marginal
distribution
"""
function z2x(d::Distribution, x::AbstractArray)
    quantile.(d, cdf.(Normal(0,1), x))
end

