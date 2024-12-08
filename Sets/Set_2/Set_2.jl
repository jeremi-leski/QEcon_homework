using ApproxFun, Calculus, DelimitedFiles, Example, FastGaussQuadrature, Interpolations, IterativeSolvers, LaTeXStrings, LinearAlgebra, NLsolve, Optim, Plots, PrettyTables, Roots, Statistics 

#Problem 3

#Step 1

# Inf?

function NPV(r,C)
    NPV_sum=0
    for t in 0:100
        C=randn(Float64)
        NPV_t=((1/(r+1))^t)*C
        NPV_sum += NPV_t
    end
return NPV_sum
end

NPV(0.2,2)


#Step 2

# r[1]?, not sure if nlsolve gives back right answers

function f_zero(r)
    return NPV(r,1)
end

quess = [1.0]
r0=nlsolve(r -> f_zero(r[1]), quess, ftol=1e-14, show_trace=true)


#Step 3

# Did I have to include formulas of NPV and f_zero into the internal_rate function or
# it was possible to run it without fully copying them, with just putting NPV(r, C) and f_zero(R) somehow?

function internal_rate(C)
    function NPV(r,C)
        NPV_sum=0
        for t in 0:100
            C=randn(Float64)
            NPV_t=((1/(r+1))^t)*C
            NPV_sum += NPV_t
        end
    return NPV_sum
    end
    function f_zero(r)
        return NPV(r,1)
    end
    quess = [0.0]
    r0=nlsolve(r -> f_zero(r[1]), quess, ftol=1e-14, show_trace=true)
return r0
end

internal_rate([1,2,3,4,5])



#Problem 4

#Step 3

function cost_function(x, α, σ, w1, w2, y)
    x1, x2 = x
    return w1 * x1 + w2 * x2
end

function ces_production(x1, x2, α, σ)
    if σ == 1
        return x1^α * x2^(1-α)
    else
        return ((α * x1^σ + (1-α) * x2^σ)^(1/σ))
    end
end

function solve_minimization(α, σ, w1, w2, y)
    objective(x) = cost_function(x, α, σ, w1, w2, y)
    constraints(x) = ces_production(x[1], x[2], α, σ) - y
    initial_guess = [1.0, 1.0]
    lower_bounds = [0.0, 0.0]
    upper_bounds = [Inf, Inf]
    result = optimize(objective, lower_bounds, upper_bounds, initial_guess, Fminbox())
    x_opt = result.minimizer
    return objective(x_opt), x_opt[1], x_opt[2]
end


cost, x1, x2 = solve_minimization(0.5, 0.8, 1.0, 1.0, 10.0)

