using Pkg
using QuadGK

#The agents maximize the following expected utility:
function expected_utility(consumption, beta, gamma)
    # Ensure beta is within the interval (0, 1)
    if beta <= 0 || beta >= 1
        error("The discount factor beta must be in the interval (0, 1).")
    end

    total_utility = 0.0
    discount_factor = 1.0
    for c in consumption
        total_utility += discount_factor * ((c^(1 - gamma) - 1) / (1 - gamma))
        discount_factor *= beta
    end
    return total_utility
end

#The period utility function
function utility(consumption, gamma)
    return (consumption^(1 - gamma) - 1) / (1 - gamma)
end

#The agents face the following budget constraint:
function budget_constraint(consumption, assets_next, income, tax, rate_of_return)
    return consumption + assets_next - (income - tax + (1 + rate_of_return) * assets_next)
end

#The post-tax labor income
function post_tax_income(y, τ, ȳ, λ)
    return (1 - τ) * y^(1 - λ) * ȳ^λ
end

#The pre tax labor income 
function pre_tax_income(z, w)
    return z * w
end

#The tax function 
function tax(y, τ, ȳ, λ)
    return y - (1 - τ) * ((y / ȳ)^(1 - λ)) * ȳ
end

#Post-tax labor income is a geo metric weighted average of the pre-tax labor income and the average labor income in the economy
function post_tax_income_weighted_average(y, τ, ȳ, λ)
    return (1 - τ) * ^(1 - λ) * ȳ^λ
end

#The idiosyncratic productivity level follows an AR(1) process
function productivity_ar1(z, ρ, ϵ, z̄)
    return exp(ρ * log(z) + (1 - ρ) * log(z̄) + ϵ)
end

#The borrowing constraint
# The borrowing constraint
function borrowing_constraint(a, ϕ)
    if ϕ <= 0
        error("ϕ must be greater than zero")
    end
    return a >= -ϕ
end

# The production function
function production(K, L, A, α)
    if α <= 0 || α >= 1
        error("The capital share α must be in the interval (0, 1).")
    end
    return A * K^α * L^(1 - α)
end

#Production function, profit maximization
function firms_problem(K, L, A, α, r, δ, w)
    F_K_L = A * K^α * L^(1 - α)
    profit = F_K_L - (r + δ) * K - w * L
    return profit
end

#Conditions
function rental_rate_of_capital(A, α, K, L, δ)
    return α * A * (K^(α - 1)) * (L^(1 - α)) - δ
end

function wage_rate(A, α, K, L)
    return (1 - α) * A * (K^α) * (L^(-α))
end

#The government uses the tax revenue to finance purchases of goods G
function integrate_tax_function(τ, ȳ, λ)
    G= quadgk(y -> tax_function(y, τ, ȳ, λ), 0, 1)
    return G
end

#The asset market clearing condition in this economy is given by:
function integrate_assets_next()
    K= quadgk(a_next -> a_next, 0, 1)
    return K
end

# The goods market clearing condition is:
function integrand(consumption, G, δ, K)
    return consumption + G + δ * K
end
function integrate_consumption_G_δ_K(G, δ, K)
    Production= quadgk(consumption -> integrand(consumption, G, δ, K), 0, 1)
    return Production
end

#The labor market clears when:
function integrand_z(z)
    return z
end
function integrate_z()
    integral, _ = quadgk(z -> integrand_z(z), 0, 1)
    return integral
end