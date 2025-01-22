using ApproxFun, Calculus, DelimitedFiles, Example, FastGaussQuadrature, Interpolations, IterativeSolvers, LaTeXStrings, LinearAlgebra, NLsolve, Optim, Plots, PrettyTables, Roots, Statistics, Parameters, Interpolations, Plots, LinearAlgebra, QuantEcon, Distributions, Random


#Task 1

X = 100.0 
p_min = 10.0 
p_max = 50.0 
f = 2.0 
q = 0.2 
nv = 50 
price_step = 0.1 

prices = p_min:price_step:p_max
prob_price = 1.0 / length(prices) 

v = zeros(nv + 1)
v_B = Dict()  
for n in nv👎1
     v_T = -f * n
    for p in prices
        v_B[(n, p)] = X - p - f * n
    end

    expected_v_B = sum(prob_price * v_B[(n, p)] for p in prices)
    v_A = -f + q * expected_v_B + (1 - q) * v[n + 1]

    
    v[n] = max(v_T, v_A)
end

println(v)


#Task 2

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



#Task 3

c = 1.0
β = 0.9
wages = 1:10
π = fill(1/10, 10)

function reservation_wage(p; max_iter=1000, tol=1e-6)
    V_U = fill(c / (1 - β), length(wages))  
    V_E = fill(0.0, length(wages))        
    w_star = 0                             

    for _ in 1:max_iter
        V_U_new = similar(V_U)
        V_E_new = similar(V_E)

      
        for i in 1:length(wages)
            w = wages[i]
            V_E_new[i] = w + β * ((1-p) * V_E[i] + p * dot(V_U, π))
            V_U_new[i] = max(V_E_new[i], c + β * dot(V_U, π))
        end

    
        if maximum(abs.(V_U_new .- V_U)) < tol
            V_U, V_E = V_U_new, V_E_new
            break
        end
        V_U, V_E = V_U_new, V_E_new
    end

   
    for i in 1:length(wages)
        if V_E[i] > c + β * dot(V_U, π)
            w_star = wages[i]
            break
        end
    end

    return w_star
end


ps = 0:0.1:1  
reservation_wages = [reservation_wage(p) for p in ps]
qs = [sum(π[i] for i in findall(w -> w >= reservation_wage(p), wages)) for p in ps]
durations = [1 / q for q in qs]

plot(ps, reservation_wages, xlabel="p", ylabel="Reservation Wage w*", title="Reservation Wage vs p", legend=false)

plot!(ps, qs, xlabel="p", ylabel="Job Acceptance Probability q", title="Job Acceptance Probability vs p", legend=false)

plot!(ps, durations, xlabel="p", ylabel="Expected Duration of Unemployment", title="Expected Duration vs p", legend=false)





#Task 4 

Z = [1, 2, 3] 
    X = [0, 1, 2, 3, 4, 5] 
    P = [
        0.5 0.3 0.2;
        0.2 0.7 0.1;
        0.3 0.3 0.4
    ] 
    
    
    function policy(X_t, Z_t)
        if Z_t == 1
            return 0
        elseif Z_t == 2
            return X_t
        elseif Z_t == 3
            return X_t + 1 <= 4 ? X_t + 1 : 3
        end
    end
    
   
    function joint_transition_matrix(X, Z, P, policy)
        nX = length(X)
        nZ = length(Z)
        
        
        nJoint = nX * nZ
        T = zeros(nJoint, nJoint)
    
        for z1 in 1:nZ, x1 in 1:nX, z2 in 1:nZ
            x_next = policy(X[x1], Z[z1])
            x2_idx = findall(x -> x == x_next, X)[1]
    
            
            T[(z1-1)*nX + x1, (z2-1)*nX + x2_idx] += P[z1, z2]
        end
    
        return T
    end
    
    T_joint = joint_transition_matrix(X, Z, P, policy)
    
    
    function stationary_distribution(T)
        eigvals, eigvecs = eigen(T')
        idx = argmax(real(eigvals))
        dist = real(eigvecs[:, idx])
        return dist / sum(dist) 
    end
    
    stationary_dist_joint = stationary_distribution(T_joint)
    
    
    function marginal_distribution_X(stationary_dist_joint, nX, nZ)
        marginal_X = zeros(nX)
    
        for x in 1:nX
            marginal_X[x] = sum(stationary_dist_joint[(z-1)*nX + x] for z in 1:nZ)
        end
    
        return marginal_X
    end
    
    marginal_X = marginal_distribution_X(stationary_dist_joint, length(X), length(Z))
    
   
    expected_X = sum(marginal_X .* X)
    
    
    println("Joint Transition Matrix:", T_joint)
    println("Stationary Distribution (Joint):", stationary_dist_joint)
    println("Marginal Distribution of X:", marginal_X)
    println("Expected Value of X:", expected_X)