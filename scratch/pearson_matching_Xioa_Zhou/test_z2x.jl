Z = randn(Float64, 100)
n = 11
a = get_coefs(margins[1], n)

X1_1 = z2x(margins[1], Z)
X1_2 = [sum([a[r+1]*He(Zi, r) for r = 0:n]) for Zi in Z]

δ = X1_1 .- X1_2

mean(δ)
std(δ)

@time z2x(margins[1], Z)
@time [sum([a[r+1]*He(Zi, r) for r = 0:n]) for Zi in Z]
