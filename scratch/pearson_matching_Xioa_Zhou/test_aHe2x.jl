include("pc_cont_cont.jl")

function aHe2x(margin, Z, n)
    a = get_coefs(margin, n)
    [sum([a[r+1]*He(Zi, r) for r = 0:n]) for Zi in Z]
end

Z = randn(Float64, 10000)
n = 7

X1_1 = z2x(margins[3], Z)
X1_2 = aHe2x(margins[3], Z, n)
δ = X1_1 .- X1_2

mean(δ)
std(δ)

@time z2x(margins[1], Z)
@time aHe2x(margins[1], Z, n)
