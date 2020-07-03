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

function z2x(margin, x)
    quantile.(margin, cdf.(Normal(0,1), x))
end

function aHe2x(margin, Z, n)
    a = get_coefs(margin, n)
    [sum([a[r+1]*He(Zi, r) for r = 0:n]) for Zi in Z]
end

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

margins = [
    (Beta(2, 3)),
    (Normal(12, 2)),
    (Exponential(3.14))
]

Z = randn(Float64, 10000)
n = 7

X1_1 = z2x(margins[3], Z)
X1_2 = aHe2x(margins[3], Z, n)
δ = X1_1 .- X1_2

mean(δ)
std(δ)

@time z2x(margins[1], Z)
@time aHe2x(margins[1], Z, n)
