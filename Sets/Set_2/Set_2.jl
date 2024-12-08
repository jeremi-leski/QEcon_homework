using ApproxFun, Calculus, DelimitedFiles, Example, FastGaussQuadrature, Interpolations, IterativeSolvers, LaTeXStrings, LinearAlgebra, NLsolve, Optim, Plots, PrettyTables, Roots, Statistics 

#Problem 1
# function task_1(f, x0, α; tol=1e-6, maxiter=1000)

g(x) = f(x) + x

iterates = [x0]
residuals = []

for n in 1:maxiter
    g_x0 = g(x0)
    x1 = (1 - α) * g_x0 + α * x0
    abs_diff = abs(x1 - g_x0) 
    push!(residuals, abs(x1 - x0))
    push!(iterates, x1)

    if abs_diff < tol
        return (0, x1, g(x1), abs_diff, iterates, residuals)
    end

    x0 = x1
end

return (1, NaN, NaN, NaN, iterates, residuals)


function somerandomfunction(x)
return (x + 1)^(1/3) - x
end

x0 = 1.0
α= 0.0
tol = 1e-6
maxiter = 1000


result = task_1(somerandomfunction, x0, α; tol=tol, maxiter=maxiter)

println("Flag: ", flag)
println("Solution: ", solution)
println("Value at Solution: ", value_at_solution)
println("Absolute Difference: ", abs_diff)
println("Iterates: ", iterates)
println("Residuals: ", residuals)



###if x0 is near the fixed point the method converges faster.





#Problem 2

function exact_solution(α,β)
    x5=1
    x4=x5
    x3=x4
    x2=x3
    x1=x2-α-x4*(α-β)-x5*(β)
    return [x1,x2,x3,x4,x5]

    
    
end
exact_solution(4,7)
function backlslash_x(α,β) 
A=[1 -1 0 α-β β ; 0 1 -1 0 0 ; 0 0 1 -1 0 ; 0 0 0 1 -1 ; 0 0 0 0 1]
b=[α ; 0 ; 0 ; 0 ; 1]
x = A\ b
return x
end
backlslash_x(2,8) 

function relative_residual(α,β)
    A=[1 -1 0 α-β β ; 0 1 -1 0 0 ; 0 0 1 -1 0 ; 0 0 0 1 -1 ; 0 0 0 0 1]
b=[α ; 0 ; 0 ; 0 ; 1]
return norm(A*backlslash_x(α,β)-b)/norm(b)
end
relative_residual(4,5)

function condition_numerator(α,β)
    A=[1 -1 0 α-β β ; 0 1 -1 0 0 ; 0 0 1 -1 0 ; 0 0 0 1 -1 ; 0 0 0 0 1]
b=[α ; 0 ; 0 ; 0 ; 1]
 norm_A=norm(A,2)
 inv_norm_A=norm(inv(A))
 return norm_A*inv_norm_A
end

condition_numerator(3,16)


#Problem 3

#Step 1

# how would I calculate it to infinity? I wrote it to 100 for sample test.

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

# not sure why [1] in r[1], but it works?, not sure if nlsolve gives back right answers

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

#Steps 1,2
function task_4(x1, x2, α, σ)
    if σ == 1
        return (x1^α) * (x2^(1 - α))
    else
        CES= (α * (x1^((σ - 1) / σ)) + (1 - α) * (x2^((σ - 1) / σ))) ^ (σ / (σ - 1))
        return CES ^ (σ / (σ - 1))
    end
end


α = 0.5 
σ_values = [0.25, 1, 4]
x1_range = range(0.1, 10, length=100)
x2_range = range(0.1, 10, length=100)


x1 = repeat(x1_range, 1, length(x2_range))
x2 = repeat(x2_range', length(x1_range), 1)


contour(layout=(1, 3), size=(900, 300))

for (i, σ) in enumerate(σ_values)

    f_values = task_4.(x1, x2, α, σ)


    contour!(x1_range, x2_range, f_values, 
          xlabel="Input x1", ylabel="Input x2", 
          title="σ = $σ", color=:viridis, 
          levels=20, subplot=i, legend=false)
end


contour!(title="CES Production Function Contours (α = $α, Various σ)")


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

