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

end

@variables x, y, l, r
vars = (x=x, y=y, l=l, r=r)
model = Model(vars)

