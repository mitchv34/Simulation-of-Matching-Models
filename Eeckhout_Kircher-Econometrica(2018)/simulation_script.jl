using  ModelingToolkit, OrdinaryDiffEq, Distributions, Parameters
using Plots, ColorSchemes

#######################################################################################
############################## Auxiliary Functions ####################################
#######################################################################################

# This function will return all the relevant derivatives of the production function
function return_relevant_derivatives(F::Num, x::Num, y::Num, l::Num, r::Num)
    # Define the differential operators
    Dx = Differential(x)
    Dl = Differential(l)
    Dr = Differential(r)
    Dy = Differential(y)
    
	Fx  = expand_derivatives(     Dx(F)   ) 
	Flr = expand_derivatives( Dl( Dr(F) ) ) 
	Fxr = expand_derivatives( Dx( Dr(F) ) ) 
	Fxy = expand_derivatives( Dx( Dy(F) ) ) 
	Fyl = expand_derivatives( Dy( Dl(F) ) ) 

    return Fx, Fxy, Flr, Fxr, Fyl
end 

function convert_to_solution!(solution::Solution, sol::ODESolution, eval_profits, eval_wages)
    
    # Create a matrix of solution
	solution.x = sol.t
    solution.θ = sol[θ]
    solution.μ = sol[μ]

	# Calculate wages and add them to the solution :
	wages = map(i -> eval_wages(matrix_sol[i, :]), 1:nk)
	solution.w = wages
	# Calculate profits and put them in the matrix:
	profits = map(i -> eval_profits(matrix_sol[i, :]), 1:nk)
	solution.Π = profits

end

#######################################################################################
#######################################################################################
#######################################################################################

#######################################################################################
############################# Model Structures ########################################
#######################################################################################

# The following Structure is used to store the model
@with_kw struct Model
    
    # Variables and parameters
    vars::NamedTuple # Variables
    params::NamedTuple # Parameters
    
    # production function
    F::Num # The production function

    # Distribution of types
    type_dist :: NamedTuple # Distribution of types

    # Model solver
    solver::Solver # Model solver

end

# Structure containing the model solution
mutable struct Solution
    x::Vector{Float64} # Worker types
    μ::Vector{Float64} # Equilibrium Matching μ(x)
    θ::Vector{Float64} # Equilibrium Firm size θ(x)
    w::Vector{Float64} # Equilibrium wages w(x)
    Π::Vector{Float64} # Equilibrium profits Π(x, μ(x))    
end


mutable struct Solver

    # Structure containing the model solver

    ode_system::ODESystem # System of Ordinary Differential Equations

    # Solver parameters
    xspan::Vector{Float64} # Interval of integration
    inital_condition :: Vector{Pair{Num, Int64}} # Initial condition
    nk :: Int # Number of points in the grid

    # Solution
    solution::Solution # Solution of the model

    # Auxiliary functions
    eval_profits::(Vector{Float64} -> Float64) # Function to evaluate profits
    eval_wages::(Vector{Float64} -> Float64) # Function to evaluate wages
end


#######################################################################################
#######################################################################################
#######################################################################################

#######################################################################################
####################################### Model Initialization ##########################
#######################################################################################

# This function initializes the model
function model = initialize_model(param_values)
    
    # Defining the variables of the system
    @variables x, y, l, r # Independent variables
    @parameters ω_A, ω_B, σ_A # Parameters
    @variables μ(x), θ(x), w(x) # Dependent variables for solution
    
    # Define the production function

    vars = (x=x, y=y, l=l, r=r)
    params = (ω_A=ω_A, ω_B=ω_B, σ_A=σ_A)

    # define the part of the production function that depends on x and y
    A = ((ω_A * x^((σ_A - 1) / σ_A) + (1 - ω_A) * y^((σ_A - 1) / σ_A))^(σ_A / (σ_A - 1))) 
    # define symbolic expression for th epart that dependes on  l and r
    B = l^ω_B * r^(1-ω_B)

    # Crete a symbolic version of the production function
    F = A * B
    
    x_bounds = [1.0 100.0] # Bounds for the support of x Distribution
    y_bounds = [1.0 100.0] # Bounds for the support of y Distribution

    # Create the distribution objects
    wuni = Uniform(x_bounds[1], x_bounds[2]) 
    funi = Uniform(y_bounds[1], y_bounds[2])

    type_dist = (workers = wuni, firms = funi)
    
    # Define the assortativity of the model
    # TODO: Infer the assortativity from the production function
    assortativity = "positive"

    # Create the solver of the model
    H = 100 * (wuni.b - wuni.a)/(funi.b - funi.a)

    # Since we are in equilibrium:
    # the variable associated with firms (y) should be substituted with μ
    # l should be θ and r = 1.0 
    f = substitute(F, Dict( y => μ, l => θ, r => 1.0))	

    # Helper function to substitute in the derivatives
    sub(fun) = substitute(fun, Dict( y => μ, l => θ, r => 1.0))	
    # Getting relevant derivatives in terms of x, μ(x) and θ(x)
    Fx, Fxy, Flr, Fxr, Fyl = sub.(return_relevant_derivatives(F, x,y,l,r))

    # Functions to evaluate the profits of the firm 
    costs = θ * w # Costs of the firm
    Π = f - costs # Profit of the firm

    # System of ODE's that solves the model
	D = Differential(x)
	eqs = []
    if assortativity == "positive"
        push!(eqs, D(θ) ~ (H * Fyl - Fxr) / Flr)
        push!(eqs, D(μ) ~ H / θ)
    elseif assortativity == "negative"
        push!(eqs, D(θ) ~ -(H * Fyl + Fxr) / Flr)
        push!(eqs, D(μ) ~ -H / θ)
    end

    # Creating systems of ordinary differential equations
 	@named sys = ODESystem(eqs)
 
     # Construct functions to evaluate wages and profits
     eval_profits = eval( build_function( substitute(Π , p), [x, θ, μ, w]) )
     Dθ = Differential(θ)
     eval_wages = eval( build_function(substitute(expand_derivatives(Dθ(f)), p), [x,θ,μ]) )

     nk = 6000 # Number of points in the grid

    # Initialize empty solution
    empty_sol = zeros(nk, 5)
    solution = Solution(empty_sol...)
    # Initialize initial conditions
    x0 = y_bounds[ ( assortativity == "positive") ? 1 : 2 ]
    initial_condition = [x0, nothing] # Initial condition for θ(x) will be defined later
    # Define the solver of the model
    solver = Solver(sys, x_bounds, initial_condition, nk, solution, eval_profits, eval_wages)
    
    # Define the model
    model = Model(vars, params, F, type_dist, assortativity)
end

#######################################################################################
############################# Functions To Solve the Model ############################
#######################################################################################


# This function solves the IVP given an initial contions and parameter values
function solve_IVP(sys, u0, xspan, p; saveat=[])
	# Defining the ODEProblem object
	prob = ODEProblem(sys, u0, xspan, p, jac=true) 
	# Solving the system
	if length(saveat) == 0
		sol = solve(prob,Tsit5())
	else
		sol = solve(prob,Tsit5(), saveat=saveat)
	end
	return sol
end


# This function iteratively update the initial guess until convergence
# This is a shooting method to solve the IVP
function shooting_method(solver, param_values; tol=1e-5, max_iter = 100, saveat=[])

    sys = solver.ode_system
    xspan = solver.xspan
	u0 = solver.initial_condition


	# Define an initial guess for firm size upper and lower bound
	firm_size_lower = 0
	firm_size_upper = 10000
	u0[2][2] = (firm_size_lower + firm_size_upper)/2

	err = 100
	# Possitive assortative means that μ(1) = 1 we will use μ(100) = 100 to check convergence
	# This should be made more general in a .jl file.
    
    sol = solve_IVP(sys, u0, xspan, param_values)

	n = 0
	println("θ(1) = $guess_firm_size")
	while (abs(err) > tol) & (n < max_iter)
		n += 1
		sol = solve_IVP(sys, u0, xspan, param_values)
		err = sol(100)[2] - 100
		if err < 0
			print("err < 0 ")
			firm_size_upper = guess_firm_size
			guess_firm_size = (firm_size_lower + firm_size_upper)/2
			println("next θ(1) = $guess_firm_size")
		else
			print("err > 0 ")
			firm_size_lower = guess_firm_size
			guess_firm_size = (firm_size_lower + firm_size_upper)/2
			println("next θ(1) = $guess_firm_size")
		end
		u0[2][2] = guess_firm_size
		# println("Iteration: $n, err=$err, θ(1) = $guess_firm_size, μ(100) =$(sol(100)[2])")
	end
	
    if length(saveat) == 0
        return sol
    else
        sol = solve_IVP(sys, u0, xspan, param_values, saveat=saveat)
    end
end

# This function solves the model using the shooting method
function Solution!(param_values::Dict())
    
    model = initialize_model(param_values)

    save_range = range(model.x_bounds[1], stop=model.x_bounds[2], length = model.nk)

    sol = shooting_method(model.solver, param_values; saveat=save_range)

    eval_profits = model.eval_profits
    eval_wages = model.eval_wages
    convert_to_solution!(model.solver.solution,  sol, eval_profits, eval_wages)
end

