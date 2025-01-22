using ApproxFun, Calculus, DelimitedFiles, Example, FastGaussQuadrature, Interpolations, IterativeSolvers, LaTeXStrings, LinearAlgebra, NLsolve, Optim, Plots, PrettyTables, Roots, Statistics


# There is now probability p that the worker will lose the job at
# the end of the period and become unemployed next period. If this
# happens, the worker gets to draw a new wage offer at the beginning
# of the next period. Let VU(w) be the value of being unemployed and
# receiving a job offer that pays wage w and VE(w) be the value of
# being employed and earning w.
# The Bellman equations are as follows:
# VU(w) = max{VE(w),c+β∑VU(w′) π(w′)}
# VE(w) = w+β [(1−p)VE(w)+p∑VU(w′) π(w′)]

#Create a plot that shows how the reservation wage w∗ changes with p.


β = 0.95  # Discount factor
c = 1.0   # Constant term in the Bellman equation

function π_1(w)
    return 1 / 10
end


function VU(w, VE_func)
    return max(VE_func(w), c + β * sum(π*(w′) * VE_func(w′) for w′ in 0:0.1:10))
end

function VE(w, p, VU_func)
    return w + β * ((1 - p) * VU_func(w) + p * sum(π*(w′) * VU_func(w′) for w′ in 0:0.1:10))
end

function reservation_wage(p)
    VE_values = Dict{Float64, Float64}()
    VU_values = Dict{Float64, Float64}()
    
    wage_grid = 0:0.1:10
    
    for w in wage_grid
        VE_values[w] = w
        VU_values[w] = c
    end
    
    for _ in 1:1000
        new_VE_values = Dict{Float64, Float64}()
        new_VU_values = Dict{Float64, Float64}()
        
        for w in wage_grid
            new_VE_values[w] = VE(w, p, w -> VU_values[w])
            new_VU_values[w] = VU(w, w -> VE_values[w])
        end
        
        VE_values = new_VE_values
        VU_values = new_VU_values
    end
    
    w_star = 0.0
    while VU_values[w_star] > VE_values[w_star]
        w_star += 0.01
    end
    return w_star
end

p_values = 0:0.01:1.0

w_star_values = [reservation_wage(p) for p in p_values]

plot(p_values, w_star_values, xlabel="Probability p", ylabel="Reservation Wage w*", title="Reservation Wage vs Probability p", legend=false)