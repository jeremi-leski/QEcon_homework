using ApproxFun, Calculus, DelimitedFiles, Example, FastGaussQuadrature, Interpolations, IterativeSolvers, LaTeXStrings, LinearAlgebra, NLsolve, Optim, Plots, PrettyTables, Roots, Statistics



# Step 1 

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

#Step 2

β = 0.95  # Discount factor
c = 1.0   # Constant term in the Bellman equation
w_max = 10.0  # Maximum possible wage

function π_w(w_max)
    return 1 / w_max  
end

function VU(w, VE_func)
    return max(VE_func(w), c + β * sum(π*(w′) * VE_func(w′) for w′ in 0:0.1:w_max))
end

function VE(w, p, VU_func)
    return w + β * ((1 - p) * VU_func(w) + p * sum(π*(w′) * VU_func(w′) for w′ in 0:0.1:w_max))
end

function reservation_wage(p)
    VE_values = Dict{Float64, Float64}()
    VU_values = Dict{Float64, Float64}()
    
    wage_grid = 0:0.1:w_max
    
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

q_values = [1 - w_star / w_max for w_star in w_star_values]

plot(p_values, w_star_values, xlabel="Probability p", ylabel="Reservation Wage w*", title="Reservation Wage vs Probability p", legend=false)

plot(p_values, q_values, xlabel="Probability p", ylabel="Probability q", title="Probability of Accepting a Job vs Probability p", legend=false)


#Step 3


β = 0.95  # Discount factor
c = 1.0   # Constant term in the Bellman equation
w_max = 10.0  # Maximum possible wage

function π_(w_max)
    return 1 /w_max 
end     

function VU(w, VE_func)
    return max(VE_func(w), c + β * sum(π*(w′) * VE_func(w′) for w′ in 0:0.1:w_max))
end

function VE(w, p, VU_func)
    return w + β * ((1 - p) * VU_func(w) + p * sum(π*(w′) * VU_func(w′) for w′ in 0:0.1:w_max))
end

function reservation_wage(p)
    VE_values = Dict{Float64, Float64}()
    VU_values = Dict{Float64, Float64}()
    
    wage_grid = 0:0.1:w_max
    
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

q_values = [1 - w_star / w_max for w_star in w_star_values]

expected_duration_values = [1 / q for q in q_values]

plot(p_values, w_star_values, xlabel="Probability p", ylabel="Reservation Wage w*", title="Reservation Wage vs Probability p", legend=false)

plot(p_values, q_values, xlabel="Probability p", ylabel="Probability q", title="Probability of Accepting a Job vs Probability p", legend=false)

plot(p_values, expected_duration_values, xlabel="Probability p", ylabel="Expected Duration of Unemployment", title="Expected Duration of Unemployment vs Probability p", legend=false)