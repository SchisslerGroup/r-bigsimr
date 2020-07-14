include("common.jl")

using Test

function Hϕ(x::T, n::Int) where T<:Real
    if isinf(x)
        zero(T)
    else
        He(x, n) * pdf(Normal(), x)
    end
end


function Gn0m(n::Int, A, α, dB, σAσB_inv)
    M = length(A)
    accu = 0.0
    for r=1:M
        accu += A[r] * (Hϕ(α[r+1], n-1) - Hϕ(α[r], n-1))
    end
    m = n + 4
    t, w = gausshermite(m)
    X = z2x(dB, t * √2)
    S = (1 / √π) * sum(w .* He(t * √2, n) .* X)
    -σAσB_inv * accu * S
end


function solvePoly_pmOne(coef)
    n = length(coef) - 1
    P(x) = Polynomial(coef)(x)
    dP(x) = Polynomial((1:n) .* coef[2:end])(x)
    r = roots(P, dP, -1..1)
    if length(r) == 0
        NaN
    else
        mid(r[1].interval)
    end
end


function rho_z_mixed(ρx, dA::DiscreteUnivariateDistribution, dB::ContinuousUnivariateDistribution, n::Int=5)
    σA = std(dA)
    σB = std(dB)
    minA = minimum(dA)
    maxA = maximum(dA)
    TA = eltype(dA)
    maxA = isinf(maxA) ? TA(quantile(dA, 0.995)) : maxA
    A = range(minA, maxA, step=1.0)
    α = [-Inf; quantile.(Normal(), cdf.(dA, A))]

    c2 = 1 / (σA * σB)

    coef = [Gn0m(i, A, α, dB, c2) / factorial(i) for i=1:n]
    coef = [-ρx; coef]

    r = solvePoly_pmOne(coef)
    if isnan(r)
        ρx_l, ρx_u = ρz_bounds(dA, dB, n)
        clampcor(clamp(ρx, ρx_l, ρx_u))
    else
        r
    end
end

rho_z_mixed(ρx, dA::ContinuousUnivariateDistribution, dB::DiscreteUnivariateDistribution, n::Int=5) = rho_z_mixed(ρx, dB, dA, n)


dA = Beta(2, 3)
dB = Binomial(2, 0.2)
dC = Binomial(20, 0.2)

rho_z_mixed(-0.7, dA, dB)

# Mixed =======================================================================
# Values from Table 4, Col 1
@test -0.890 ≈ rho_z_mixed(-0.7, dB, dA) atol=0.005
@test -0.632 ≈ rho_z_mixed(-0.5, dB, dA) atol=0.005
@test -0.377 ≈ rho_z_mixed(-0.3, dB, dA) atol=0.005
@test  0.366 ≈ rho_z_mixed( 0.3, dB, dA) atol=0.005
@test  0.603 ≈ rho_z_mixed( 0.5, dB, dA) atol=0.005
@test  0.945 ≈ rho_z_mixed( 0.8, dB, dA) atol=0.005

# Values from Table 4, Col 2
@test -0.928 ≈ rho_z_mixed(-0.9, dC, dA) atol=0.005
@test -0.618 ≈ rho_z_mixed(-0.6, dC, dA) atol=0.005
@test -0.309 ≈ rho_z_mixed(-0.3, dC, dA) atol=0.005
@test  0.308 ≈ rho_z_mixed( 0.3, dC, dA) atol=0.005
@test  0.613 ≈ rho_z_mixed( 0.6, dC, dA) atol=0.005
@test  0.916 ≈ rho_z_mixed( 0.9, dC, dA) atol=0.005
