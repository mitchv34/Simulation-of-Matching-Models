using NLsolve
using Plots
using Parameters
using DataFrames
using Term # You dont really need this package
using Latexify

mutable struct Params

    s₁::Float64  
	s₂::Float64
	p ::Float64
	b ::Float64
	β ::Float64
	δ ::Float64
	c ::Float64
	r ::Float64

    m::Function

    function Params(s₂::Float64)
        s₁ = 1;
        p = 2/3;
        b = 0.1;
        β = 0.5;
        δ = 0.2;
        c = 0.3;
        r = 0.05;

        m = (θ) -> 2√θ

        new(s₁, s₂, p, b, β, δ, c, r, m)
    end # constructor

end # Parameters

# Wage equation in equilibirum (takes value of unemployment as given ⟹ no need for s as argument)
w(y, rU, params) = params.β*(y-params.c) + (1-params.β) * rU

struct Equilibrium4 # There are missing variables 

    type::String

    s₁  ::Float64
    s₂  ::Float64
    θ   ::Float64
	m_θ ::Float64
	u   ::Float64
    γ   ::Float64
	ϕ   ::Float64
    w₁₁ ::Float64
    w₂₁ ::Float64
    w₂₂ ::Float64
    #Y   ::Float64

end


function plot_eq(list_eq::Array{Equilibrium4}; path="./Albrecht_Vroman_2002/figures/")
    
    eq_df = eq_to_df(list_eq)

    title = list_eq[1].type

    println("Ploting: $title")

    cols = names(eq_df);

    x = cols[2]
    ys = cols[3:end]

    for y in ys
        println("\t Generating plot for $(y)")
        plot(eq_df[:, x], eq_df[:, y], title = title, lw = 2, legend = false)
        xlabel!("$(latexify(x))")
        ylabel!("$(latexify(y))")
        savefig(path*"$(title) $(y)")
    end

end

function eq_to_df(list_eq::Array{Equilibrium4})

    s₁_list = []; s₂_list = []; θ_list = []; m_θ_list = [];
    u_list = []; γ_list = []; ϕ_list = [];
    w₁₁_list = []; w₂₁_list = []; w₂₂_list = [];# Y_list = [];

    n_eq = length(list_eq)

    for i ∈ 1:n_eq
        push!(s₁_list, list_eq[i].s₁); push!(s₂_list, list_eq[i].s₂);
        push!(θ_list, list_eq[i].θ); push!(m_θ_list, list_eq[i].m_θ);
        push!(u_list, list_eq[i].u); push!(γ_list, list_eq[i].γ); push!(ϕ_list, list_eq[i].ϕ);
        push!(w₁₁_list, list_eq[i].w₁₁); push!(w₂₁_list, list_eq[i].w₂₁); 
        push!(w₂₂_list, list_eq[i].w₂₂);# push!(Y_list, list_eq[i].Y, );
    end



    return DataFrame(
        [  
            :s₁  => s₁_list,
            :s₂  => s₂_list,
            :θ   => θ_list,
            :m_θ => m_θ_list,
            :u   => u_list,
            :γ   => γ_list,
            :ϕ   => ϕ_list,
            :w₁₁ => w₁₁_list,
            :w₂₁ => w₂₁_list,
            :w₂₂ => w₂₂_list,
            # :Y   => Y_list,
        ]
    )
end

function cs_equilibirum(s_2::Float64)

    params = Params(s_2)
    @unpack s₁ ,s₂ ,p  ,b  ,β  ,δ  ,c  ,r, m = params
    
    # Solve for the value of theta
    # define the equation that would give us theta
    function θ!(F, x) 
        F[1] = c - ( m(x[1]) * (1 - β )*(s₁ - c - b)) / ( x[1] * (r + δ + β*m(x[1])) ) ;
    end
    #  solve for theta
    θ_eq = nlsolve(θ!, [0.3]).zero[1];
    m_θ_eq = m(θ_eq);

    # solve for gamma
    #   Write the expression for phi in terms of gamma (using m_theta_eq)
    ϕ(γ) =( (1-γ)*p*m_θ_eq + (p - γ) * δ )/ (γ * m_θ_eq * (1 - p));

    # Write the equation that would give us gamma
    function γ!(F, x)
        F[1] = x[1]*(r+δ)*(s₁-c-b)-(1-x[1])*(s₂ - s₁)*(r+δ+m_θ_eq*ϕ(x[1])*β);
    end

    γ_eq = nlsolve(γ!, [0.3]).zero[1];
    ϕ_eq = ϕ(γ_eq);

    # Solve for u in equilibrium
    u_eq = δ * (1 - p) / ((1 - γ_eq)*(δ + m_θ_eq));
    
    # Value of unemployment for low skilled workers
    rUs₁ = (b*(r + δ) + m_θ_eq*β*ϕ_eq*(s₁ - c) )/(r + δ + m_θ_eq*ϕ_eq*β)
    # Value of unemployment for hihg-skilled workers
    rUs₂ = (b*(r+δ) + m_θ_eq * β *( ϕ_eq * s₁ + (1-ϕ_eq)s₂ - c))/(r+δ+m_θ_eq*β)
    
    # Calculate the value of wages 
    w₁₁ = w(s₁, rUs₁, params)
    w₂₁ = w(s₁, rUs₂, params)
    w₂₂ = w(s₂, rUs₂, params)

    # TODO: Calculate Agregate Production

    return Equilibrium4("Cross-Skill Matching", s₁ ,s₂, θ_eq, m_θ_eq, u_eq, γ_eq, ϕ_eq, w₁₁, w₂₁, w₂₂)

end


function es_equilibirum(s_2::Float64)

    params = Params(s_2)
    @unpack s₁ ,s₂ ,p  ,b  ,β  ,δ  ,c  ,r, m = params
    
    # Equation for ϕ in terms of θ and γ
    ϕ(θ,γ) = (p * (1-γ) * m(θ) + δ*(p-γ)) / (m(θ)*(γ+p-2*γ*p))
    
    # Equation for u in terms of  θ and γ
    u(θ,γ) = ( δ*(γ+p-2*γ*p)) / ( γ*(1-γ)*(m(θ)+2δ))
    
    # System of equations to obtain equilibrium values of θ and γ 
    function f!(F, x)
        θ, γ = x
        F[1] = (m(θ)/θ)*γ*( ((1-β)*(s₁-c-b)) / (r+δ+m(θ)*ϕ(θ,γ)*β)) - c
        F[2] = (m(θ)/θ)*(1-γ)*( ((1-β)*(s₂-c-b)) / (r+δ+m(θ)*(1-ϕ(θ,γ))*β)) - c
    end
    
    sol = nlsolve(f!, [0.3,0.4]);
	θ_eq, γ_eq = sol.zero
	m_θ_eq = m(θ_eq)
	ϕ_eq = ϕ(θ_eq, γ_eq)
	u_eq = u(θ_eq, γ_eq)
    
    # Value of unemployment for low skilled workers
    rUs₁ = (b*(r + δ) + m(θ_eq)*β*ϕ_eq*(s₁ - c) )/(r + δ + m(θ_eq)*ϕ_eq*β)
    # Value of unemployment for hihg-skilled workers
    rUs₂ = (b*(r + δ) + m(θ_eq)*β*(1 - ϕ_eq)* (s₂ - c) )/( r + δ + m(θ_eq)*β )
    
    # Calculate the value of wages 
    w₁₁ = w(s₁, rUs₁, params)
    w₂₁ = NaN
    w₂₂ = w(s₂, rUs₂, params)

    # TODO: Calculate Agregate Production

    return Equilibrium4("Ex-Post Segmentation", s₁ ,s₂, θ_eq, m_θ_eq, u_eq, γ_eq, ϕ_eq, w₁₁, w₂₁, w₂₂)

end

eq_list = es_equilibirum.(1.30:0.01:2.0)

plot_eq(eq_list)
-
