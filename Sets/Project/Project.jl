using LinearAlgebra, Statistics, Printf, Distributions, Plots


function calibrate_firm(alpha::Float64, r::Float64, L::Float64, w::Float64, i_to_y::Float64)
    K = (alpha - i_to_y) / (r * (1.0 - alpha))
    δ = i_to_y / (K * (1.0 - alpha))
    A = 1.0 / ((1.0 - alpha) * K^alpha)
    return (A, δ, K)
end


function tax_income(y::Float64, tau::Float64, lambda::Float64, ybar::Float64)
    return y - (1.0 - tau) * ((y / ybar)^(1.0 - lambda)) * ybar
end

function net_labor_income(y::Float64, tau::Float64, lambda::Float64; ybar::Float64=1.0)
    return (1.0 - tau) * ((y / ybar)^(1.0 - lambda)) * ybar
end


function u(c::Float64; sigma::Float64=2.0)
    if c <= 1e-14
        return -1e10
    else
        return (c^(1.0 - sigma)) / (1.0 - sigma)
    end
end


function build_asset_grid(φ::Float64, Amax::Float64, nA::Int)
    return collect(range(-φ, Amax, length=nA))
end


function tauchen(nz::Int, rho::Float64, sigma_eps::Float64; m::Float64=3.0)
    mu = -0.5 * sigma_eps^2 / (1 - rho^2)
    sigma_z = sigma_eps / sqrt(1 - rho^2)
    log_z_min = mu - m * sigma_z
    log_z_max = mu + m * sigma_z
    log_z_grid = collect(range(log_z_min, log_z_max, length=nz))
    d = (log_z_max - log_z_min) / (nz - 1)
    Π = zeros(nz, nz)
    for i in 1:nz
        for j in 1:nz
            if j == 1
                Π[i, j] = cdf(Normal(), (log_z_grid[j] - mu - rho*(log_z_grid[i]-mu) + d/2)/sigma_eps)
            elseif j == nz
                Π[i, j] = 1 - cdf(Normal(), (log_z_grid[j] - mu - rho*(log_z_grid[i]-mu) - d/2)/sigma_eps)
            else
                Π[i, j] = cdf(Normal(), (log_z_grid[j] - mu - rho*(log_z_grid[i]-mu) + d/2)/sigma_eps) -
                          cdf(Normal(), (log_z_grid[j] - mu - rho*(log_z_grid[i]-mu) - d/2)/sigma_eps)
            end
        end
    end
    zgrid = exp.(log_z_grid)
    return zgrid, Π
end


function solve_household!(
    r::Float64, w::Float64, tau::Float64, lambda::Float64,
    Agrid::Vector{Float64}, zgrid::Vector{Float64}, Π::Matrix{Float64},
    β::Float64, σ::Float64; ybar=1.0
)
    nA = length(Agrid)
    nz = length(zgrid)
    net_y = [net_labor_income(w*z, tau, lambda; ybar=ybar) for z in zgrid]

    V    = zeros(nA, nz)
    Vnew = similar(V)
    polA = zeros(Int, nA, nz)
    maxiter = 1000
    tol = 1e-6
    dist = 1.0

    R = 1.0 + r
    one_step_return = Array{Float64}(undef, nA, nA, nz)
    for iz in 1:nz
        for ia in 1:nA
            a_now = Agrid[ia]
            for iap in 1:nA
                a_next = Agrid[iap]
                c = net_y[iz] + R*a_now - a_next
                one_step_return[iap, ia, iz] = u(c; sigma=σ)
            end
        end
    end

    iter = 0
    while dist > tol && iter < maxiter
        iter += 1
        EV = V * Π'
        for iz in 1:nz
            for ia in 1:nA
                best_val = -1e20
                best_iap = 1
                for iap in 1:nA
                    val = one_step_return[iap, ia, iz] + β * EV[iap, iz]
                    if val > best_val
                        best_val = val
                        best_iap = iap
                    end
                end
                Vnew[ia, iz] = best_val
                polA[ia, iz] = best_iap
            end
        end
        dist = maximum(abs.(Vnew .- V))
        V .= Vnew
    end

    return V, polA
end


function stationary_distribution(polA::Array{Int,2}, Π::Matrix{Float64})
    nA, nz = size(polA)
    Γ = fill(1.0/(nA*nz), nA, nz)
    tol = 1e-10
    maxiter = 10_000
    dist = 1.0
    iter = 0
    Γnew = similar(Γ)
    while dist > tol && iter < maxiter
        iter += 1
        Γnew .= 0.0
        for iz in 1:nz
            for ia in 1:nA
                iap = polA[ia, iz]
                for iznext in 1:nz
                    Γnew[iap, iznext] += Γ[ia, iz] * Π[iz, iznext]
                end
            end
        end
        dist = maximum(abs.(Γnew .- Γ))
        Γ .= Γnew
    end
    return Γ
end


function compute_aggregate_assets(Γ::Matrix{Float64}, Agrid::Vector{Float64})
    nA, nz = size(Γ)
    aggA = 0.0
    for iz in 1:nz, ia in 1:nA
        aggA += Γ[ia, iz] * Agrid[ia]
    end
    return aggA
end


function compute_aggregate_consumption(Γ::Matrix{Float64}, Agrid::Vector{Float64},
                                       polA::Array{Int,2}, net_y::Vector{Float64}, r::Float64)
    nA, nz = size(Γ)
    aggC = 0.0
    R = 1.0 + r
    for iz in 1:nz, ia in 1:nA
        c = net_y[iz] + R*Agrid[ia] - Agrid[polA[ia, iz]]
        aggC += Γ[ia, iz] * c
    end
    return aggC
end


function compute_government_spending(Γ::Matrix{Float64}, zgrid::Vector{Float64},
                                     tau::Float64, lambda::Float64; w::Float64=1.0, ybar::Float64=1.0)
    nA, nz = size(Γ)
    G = 0.0
    for iz in 1:nz
        y = w * zgrid[iz]
        tax_rev = tax_income(y, tau, lambda, ybar)
        mass_z = sum(Γ[:, iz])
        G += tax_rev * mass_z
    end
    return G
end


function compute_average_labor(Γ::Matrix{Float64}, zgrid::Vector{Float64})
    nA, nz = size(Γ)
    avg_z = 0.0
    for iz in 1:nz
        mass_z = sum(Γ[:, iz])
        avg_z += mass_z * zgrid[iz]
    end
    return avg_z
end




function find_beta_lambda0(; beta_min=0.90, beta_max=0.99, tol_beta=1e-4, max_iter=50,
                           r=0.04, w=1.0, tau=0.3125, lambda=0.0,
                           alpha=0.36, i_to_y=0.20, L=1.0, φ=0.0,
                           nA=50, Amax=100.0, nz=5, sigma=2.0, ybar=1.0,
                           rho_z=0.9, sigma_eps=0.4)
    iter = 0
    beta_low = beta_min
    beta_high = beta_max
    beta_mid = (beta_low + beta_high) / 2
    diff = Inf
    while iter < max_iter && abs(diff) > tol_beta
        iter += 1
        (A, δ, K) = calibrate_firm(alpha, r, L, w, i_to_y)
        Agrid = build_asset_grid(φ, Amax, nA)
        zgrid, Π = tauchen(nz, rho_z, sigma_eps)
        V, polA = solve_household!(r, w, tau, lambda, Agrid, zgrid, Π, beta_mid, sigma; ybar=ybar)
        Γ = stationary_distribution(polA, Π)
        A_demand = compute_aggregate_assets(Γ, Agrid)
        diff = A_demand - K
        if diff > 0
            beta_high = beta_mid
        else
            beta_low = beta_mid
        end
        beta_mid = (beta_low + beta_high) / 2
        @printf("Iteration %d: beta = %.6f, asset market error = %.6e\n", iter, beta_mid, diff)
    end
    return beta_mid
end


function solve_equilibrium_lambda(new_lambda::Float64, beta::Float64;
    r_init=0.04, tau_init=0.3125, tol_eq=1e-4, max_iter=50,
    w_fixed=1.0, alpha=0.36, i_to_y=0.20, L=1.0, φ=0.0,
    nA=50, Amax=100.0, nz=5, sigma=2.0, ybar=1.0,
    rho_z=0.9, sigma_eps=0.4)
   
    r_current = r_init
    tau_current = tau_init
    step_r = 0.001
    step_tau = 0.001
    iter = 0
    A = 0.0
    δ = 0.0
    K = 0.0
    Y = 0.0
    V = zeros(nA, nz)
    polA = zeros(Int, nA, nz)
    Γ = zeros(nA, nz)
    zgrid = zeros(nz)
    Π = zeros(nz, nz)
    aggC = 0.0
    G = 0.0
    while iter < max_iter
        iter += 1
        (A, δ, K) = calibrate_firm(alpha, r_current, L, w_fixed, i_to_y)
        Y = A * K^alpha * L^(1.0 - alpha)
        Agrid = build_asset_grid(φ, Amax, nA)
        zgrid, Π = tauchen(nz, rho_z, sigma_eps)
        V, polA = solve_household!(r_current, w_fixed, tau_current, new_lambda, Agrid, zgrid, Π, beta, sigma; ybar=ybar)
        Γ = stationary_distribution(polA, Π)
        A_demand = compute_aggregate_assets(Γ, Agrid)
        aggC = compute_aggregate_consumption(Γ, Agrid, polA, [net_labor_income(w_fixed*z, tau_current, new_lambda; ybar=ybar) for z in zgrid], r_current)
        G = compute_government_spending(Γ, zgrid, tau_current, new_lambda; w=w_fixed, ybar=ybar)
        f1 = A_demand - K           
        f2 = (G / Y) - 0.2          
        @printf("Iter %d: r = %.4f, tau = %.4f, f1 = %.4e, f2 = %.4e\n", iter, r_current, tau_current, f1, f2)
        if abs(f1) < tol_eq && abs(f2) < tol_eq
            break
        end
        r_current -= step_r * f1
        tau_current -= step_tau * f2
    end
    result = Dict(
        :r => r_current,
        :tau => tau_current,
        :A => A,
        :δ => δ,
        :K => K,
        :Y => Y,
        :V => V,
        :polA => polA,
        :Gamma => Γ,
        :zgrid => zgrid,
        :Π => Π,
        :aggC => aggC,
        :G => G
    )
    return result
end



function gini_coefficient(values::Vector{Float64}, weights::Vector{Float64})
    inds = sortperm(values)
    x = values[inds]
    w = weights[inds]
    w = w / sum(w)
    cw = cumsum(w)
    cx = cumsum(x .* w)
    μ = sum(x .* w)
    B = sum(cx .* w)
    gini = 1 - 2 * B / μ
    return gini
end

function lorenz_curve(values::Vector{Float64}, weights::Vector{Float64})
    inds = sortperm(values)
    x = values[inds]
    w = weights[inds]
    w = w / sum(w)
    cumw = cumsum(w)
    cumx = cumsum(x .* w)
    μ = sum(x .* w)
    lorenz = cumx / μ
    return cumw, lorenz
end



function main()
   
    beta_target = find_beta_lambda0(r=0.04, w=1.0, tau=0.3125, lambda=0.0,
                                      alpha=0.36, i_to_y=0.20, L=1.0, φ=0.0,
                                      nA=50, Amax=100.0, nz=5, sigma=2.0, ybar=1.0,
                                      rho_z=0.9, sigma_eps=0.4)
    @printf("\nEquilibrium beta (λ = 0) found: %.6f\n\n", beta_target)
    
    
    (A0, δ0, K0) = calibrate_firm(0.36, 0.04, 1.0, 1.0, 0.20)
    Y0 = A0 * K0^0.36 * 1.0^(1.0 - 0.36)
    Agrid0 = build_asset_grid(0.0, 100.0, 50)
    zgrid0, Π0 = tauchen(5, 0.9, 0.4)
    V0, polA0 = solve_household!(0.04, 1.0, 0.3125, 0.0, Agrid0, zgrid0, Π0, beta_target, 2.0; ybar=1.0)
    Γ0 = stationary_distribution(polA0, Π0)
    A_demand0 = compute_aggregate_assets(Γ0, Agrid0)
    aggC0 = compute_aggregate_consumption(Γ0, Agrid0, polA0, [net_labor_income(1.0*z, 0.3125, 0.0; ybar=1.0) for z in zgrid0], 0.04)
    G0 = compute_government_spending(Γ0, zgrid0, 0.3125, 0.0; w=1.0, ybar=1.0)
    gov_ratio0 = G0 / Y0
    @printf("Economy (λ = 0):\n  K = %.4f, Aggregate asset demand = %.4f, G/Y = %.4f\n\n", K0, A_demand0, gov_ratio0)
    
    
    eq_result = solve_equilibrium_lambda(0.15, beta_target; r_init=0.04, tau_init=0.3125,
                                           w_fixed=1.0, alpha=0.36, i_to_y=0.20, L=1.0, φ=0.0,
                                           nA=50, Amax=100.0, nz=5, sigma=2.0, ybar=1.0,
                                           rho_z=0.9, sigma_eps=0.4)
    r15 = eq_result[:r]
    tau15 = eq_result[:tau]
    K15 = eq_result[:K]
    Y15 = eq_result[:Y]
    gov_ratio15 = eq_result[:G] / Y15
    @printf("Economy (λ = 0.15):\n  r = %.4f, tau = %.4f\n  K = %.4f, Y = %.4f, G/Y = %.4f\n\n", r15, tau15, K15, Y15, gov_ratio15)
    
    
    K_Y_ratio_0 = K0 / Y0
    K_Y_ratio_15 = K15 / Y15
    @printf("Capital-to-output ratio:\n  (λ = 0): %.4f\n  (λ = 0.15): %.4f\n\n", K_Y_ratio_0, K_Y_ratio_15)
    
    
    function extract_distribution(Γ::Matrix{Float64}, grid::Vector{Float64})
        nA, nz = size(Γ)
        vals = Float64[]
        weights = Float64[]
        for iz in 1:nz
            for ia in 1:nA
                push!(vals, grid[ia])
                push!(weights, Γ[ia, iz])
            end
        end
        return vals, weights
    end
    
    assets0, w_assets0 = extract_distribution(Γ0, Agrid0)
    assets15, w_assets15 = extract_distribution(eq_result[:Gamma], build_asset_grid(0.0, 100.0, 50))
    gini_assets0 = gini_coefficient(assets0, w_assets0)
    gini_assets15 = gini_coefficient(assets15, w_assets15)
    
    labor_income0 = [net_labor_income(1.0*z, 0.3125, 0.0; ybar=1.0) for z in zgrid0]
    mass0 = [sum(Γ0[:, iz]) for iz in 1:length(zgrid0)]
    gini_labor0 = gini_coefficient(labor_income0, mass0)
    
    labor_income15 = [net_labor_income(1.0*z, tau15, 0.15; ybar=1.0) for z in eq_result[:zgrid]]
    mass15 = [sum(eq_result[:Gamma][:, iz]) for iz in 1:length(eq_result[:zgrid])]
    gini_labor15 = gini_coefficient(labor_income15, mass15)
    
    @printf("Gini coefficients:\n")
    @printf("  After-tax labor income: λ = 0: %.4f, λ = 0.15: %.4f\n", gini_labor0, gini_labor15)
    @printf("  Assets:                λ = 0: %.4f, λ = 0.15: %.4f\n\n", gini_assets0, gini_assets15)
    
    
    median_z0 = Int(round(length(zgrid0)/2))
    median_z15 = Int(round(length(eq_result[:zgrid])/2))
    plt1 = plot(1:length(Agrid0), V0[:, median_z0], label="λ = 0", xlabel="Asset grid index", ylabel="Value", title="Value Function (median z)")
    plot!(plt1, 1:length(Agrid0), eq_result[:V][:, median_z15], label="λ = 0.15")
    savefig(plt1, "value_functions.png")
    
    plt2 = plot(1:length(Agrid0), Agrid0[polA0[:, median_z0]], label="λ = 0", xlabel="Asset index", ylabel="Next-period asset", title="Policy Function (median z)")
    grid15 = build_asset_grid(0.0, 100.0, 50)
    plot!(plt2, 1:length(Agrid0), grid15[eq_result[:polA][:, median_z15]], label="λ = 0.15")
    savefig(plt2, "policy_functions.png")
    
    plt3 = histogram(assets0, weights=w_assets0, label="λ = 0", xlabel="Assets", title="Marginal Distribution of Assets", alpha=0.5, nbins=20)
    histogram!(plt3, assets15, weights=w_assets15, label="λ = 0.15", alpha=0.5, nbins=20)
    savefig(plt3, "asset_distribution.png")
    
    lw0_x, lw0_y = lorenz_curve(labor_income0, mass0)
    lw15_x, lw15_y = lorenz_curve(labor_income15, mass15)
    plt4 = plot(lw0_x, lw0_y, label="λ = 0", xlabel="Cumulative share of population", ylabel="Cumulative share of income", title="Lorenz Curve: Labor Income")
    plot!(plt4, lw15_x, lw15_y, label="λ = 0.15")
    savefig(plt4, "lorenz_labor.png")
    
    la0_x, la0_y = lorenz_curve(assets0, w_assets0)
    la15_x, la15_y = lorenz_curve(assets15, w_assets15)
    plt5 = plot(la0_x, la0_y, label="λ = 0", xlabel="Cumulative share of population", ylabel="Cumulative share of assets", title="Lorenz Curve: Assets")
    plot!(plt5, la15_x, la15_y, label="λ = 0.15")
    savefig(plt5, "lorenz_assets.png")
    
    println("\nPlots saved: value_functions.png, policy_functions.png, asset_distribution.png, lorenz_labor.png, lorenz_assets.png")
end
