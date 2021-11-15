using  ModelingToolkit, OrdinaryDiffEq, Distributions, Parameters
using Plots, ColorSchemes

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

# The following Structure is used to store the model
@with_kw struct Model
    
    # Variables and parameters
    vars::NamedTuple # Variables
    params::NamedTuple # Parameters
    
    # production function
    F::Num # The production function

    # Distribution of types


end


# This function initializes the model
function model = initialize_model(vars, params)
    # Define the production function
    @variables x, y, l, r
    @parameters ω_A, ω_B, σ_A

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

    # Define the model
    model = Model(vars, params, F)
end


