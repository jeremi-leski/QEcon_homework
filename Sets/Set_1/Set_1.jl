using DelimitedFiles, Plots, Statistics

# Problem 1

function odd_or_even(n)
    if iseven(n) true 
        println("Even")

    elseif iseven(n) false
        else
        println("Odd")
    end
end

# Problem 2

function compare_three(a, b, c)
    if a>0 && b>0 && c>0 true 
        println("All numbers are positive")

    elseif a==0 && b==0 && c==0 true   
        println("All numbers are zero")

    else
        println("At least one number is not positive")

    end
end

# Problem 3

function my_factorial(n)
    output = factorial(n) 
        println("Output:" , output)
end

# Problem 4

function count_positives(arr)
    counter = 0  
    for num in arr
        if num > 0  
            counter += 1  
        end
    end
    println("Output: ", counter)  
end

# Problem 5

function plot_powers(n)
    x = -10:0.2:10  
    power_plots = plot()  

    for i in 1:n
        y = x .^ i   
        plot!(x, y, label="x^$i", linewidth=3, linestyle=:dash)  
    end

    xlabel!("x")
    ylabel!("y")
    title!("Powers of x")
    display(power_plots)  
end


# Problem 6
function standard_deviation(x)
    n=length(x)
    β1=sqrt((1/n-1)*sum((x.-(sum(x)/length(x))^2)))
  return β1
end
x=[1, 2, 5]
 standard_deviation(x)
 n
sum(x)

# Problem 7
data_homework = readdlm("data_analysis//datasets//dataset.csv",',',Float64) 
data_homework   
earnings=(data_homework[:,1])   
education=(data_homework[:,2])
hours_worked=(data_homework[:,3])
plot_1 = scatter(education,earnings ; legend=false, color=:blue, markersize = 5, opacity=0.7)
plot_2 = scatter(hours_worked,earnings ; legend=false, color=:blue, markersize = 5, opacity=0.7)
xaxis!(plot_1, "education")
yaxis!(plot_1, "earnings")
xaxis!(plot_2, "hours_worked")
yaxis!(plot_2, "earnings")
cor(earnings,education)
cor(earnings,hours_worked)



