include("common.jl")

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

get_coefs(margins[1], 15)


"""
Estimate the input correlation coefficient ρ_z given the marginal CDFs of
(Xᵢ, Xⱼ), and the desired correlation coefficient ρ_x
"""
function rho_z_cc(ρx, marginᵢ, marginⱼ; n::Int=3)
    μᵢ = mean(marginᵢ)
    σᵢ = std(marginᵢ)
    μⱼ = mean(marginⱼ)
    σⱼ = std(marginⱼ)

    # Eq (25)
    k = 0:1:n
    a = get_coefs(marginᵢ, n)
    b = get_coefs(marginⱼ, n)
    # return (μᵢ, σᵢ, μⱼ, σⱼ, a, b)

    # Eq (22)
    c1 = -μᵢ*μⱼ
    c2 =  1 / (σᵢ*σⱼ)
    kab = factorial.(k).*a.*b
    ρx_l = c2*(c1 + sum((-1).^k .* kab))
    ρx_u = c2*(c1 + sum(kab))
    return (ρx_l, ρx_u)

    coef = [c2*a[k+1]*b[k+1]*factorial(k) for k = 1:n]
    coef = [c1*c2 - c2*a[1]*b[1] - ρx; coef]
    # return coef

    return (x -> Polynomial(coef)(x), (ρx_l, ρx_u), c1, c2)
end


function rho_z_bounds_cc(marginᵢ, marginⱼ; n::Int=3)
    μᵢ = mean(marginᵢ)
    σᵢ = std(marginᵢ)
    μⱼ = mean(marginⱼ)
    σⱼ = std(marginⱼ)

    # Eq (25)
    k = 0:1:n
    a = get_coefs(marginᵢ, n)
    b = get_coefs(marginⱼ, n)
    # return (μᵢ, σᵢ, μⱼ, σⱼ, a, b)

    # Eq (22)
    c1 = -μᵢ*μⱼ
    c2 =  1 / (σᵢ*σⱼ)
    kab = factorial.(k).*a.*b
    ρx_l = c1*c2 + c2*sum((-1).^k .* kab)
    ρx_u = c1*c2 + c2*sum(kab)
    (ρx_l, ρx_u)
end

rho_z_bounds_cc(margins[1], margins[1])
rho_z_bounds_cc(margins[1], margins[2])
rho_z_bounds_cc(margins[1], margins[3])
rho_z_bounds_cc(margins[2], margins[2])
rho_z_bounds_cc(margins[2], margins[3])
rho_z_bounds_cc(margins[3], margins[3])
