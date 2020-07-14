import Statistics: mean, std, quantile, clampcor
import FastGaussQuadrature: gausshermite
import Polynomials: Polynomial
import IntervalRootFinding: roots
import Memoize: @memoize

using Distributions
using IntervalArithmetic


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


"""
    get_coefs()

Get the coefficients for the Hermite Polynomial expansion for F⁻¹[Φ(⋅)]

# Notes
The paper describes using Guass-Hermite quadrature using the Probabilists'
version of the Hermite polynomials, while the package `FastGaussQuadrature.jl`
uses the Physicists' version. Because of this, we need to do a rescaling of the
input and the output:

∑(w⋅Hₑ(t, k) / k!) ⇒ (1/√π)⋅∑(w⋅Hₑ(t√2, k) / k!)
"""
function get_coefs(margin, n)
    c = Array{Float64, 1}(undef, n+1)
    m = n+4
    t, w = gausshermite(m)
    for k = 0:1:n
        # need to do a change of variable
        X = z2x(margin, t * √2)
        c[k+1] = (1 / √π) * sum(w .* He(t * √2, k) .* X) / factorial(k)
    end
    c
end


function ρz_bounds(dA::UnivariateDistribution, dB::UnivariateDistribution; n::Int=7)
    μA = mean(dA)
    σA = std(dA)
    μB = mean(dB)
    σB = std(dB)

    # Eq (25)
    k = 0:1:n
    a = get_coefs(dA, n)
    b = get_coefs(dB, n)

    # Eq (22)
    c1 = -μA * μB
    c2 = 1 / (σA * σB)
    kab = factorial.(k) .* a .* b
    ρx_l = c1 * c2 + c2 * sum((-1) .^ k .* kab)
    ρx_u = c1 * c2 + c2 * sum(kab)

    (clampcor(ρx_l), clampcor(ρx_u))
end
