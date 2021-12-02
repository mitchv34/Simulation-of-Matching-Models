# println(pwd())
println("Loadig Packages")
using LaTeXStrings, Plots
include("simulation_functions.jl")

theme(:vibrant) 
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

println("Initializing Model")
assortativity = "positive"
model = initialize_model(assortativity)

ω_A, ω_B, σ_A = model.params # unpack parameters into variables

mutable struct Solutions_for_Plots
   # Variation of ω_A
   ω_A::Dict
   # Variation of σ_A
   σ_A::Dict
   # Variation of ω_B
   ω_B::Dict
end

solution_plot = Solutions_for_Plots(Dict(), Dict(), Dict())

param_var = Dict( :ω_A => [0.5, 0.25, 0.75], # parameter variation
                  :ω_B => [0.5, 0.25, 0.75], # parameter variation
                  :σ_A => [0.1, 0.3, 0.99]) # parameter variation



μ, θ, w = model.solver.vars
x = model.vars[:x]
Dθ = Differential(θ)



println("Running Simulations")
for i ∈ 1:3
   
   # # pack parameter values into a dictionary
   param_values_ω_A = Dict(ω_A => param_var[:ω_A][i], ω_B => 0.5, σ_A => 0.9)
   param_values_σ_A = Dict(ω_A => 0.3, ω_B => 0.6, σ_A => param_var[:σ_A][i])
   param_values_ω_B = Dict(ω_A => 0.5, ω_B => param_var[:ω_B][i], σ_A => 0.9)
   
   Solution!(model, param_values_ω_A)
   mat = hcat(model.solver.solution.x, model.solver.solution.θ, model.solver.solution.μ)
   # eval_wages = eval( build_function( substitute( expand_derivatives( Dθ( model.f ) ), param_values_ω_A ) , [x,θ,μ] ) )
   # model.solver.solution.w  = map(i -> eval_wages(mat), 1:model.solver.nk)
   solution_plot.ω_A[param_var[:ω_A][i]] = copy(model.solver.solution)
   
   Solution!(model, param_values_σ_A)
   mat = hcat(model.solver.solution.x, model.solver.solution.θ, model.solver.solution.μ)
   # eval_wages = eval( build_function( substitute( expand_derivatives( Dθ( model.f ) ), param_values_σ_A ) , [x,θ,μ] ) )
   # model.solver.solution.w  = map(i -> eval_wages(mat), 1:model.solver.nk)
   solution_plot.σ_A[param_var[:σ_A][i]] = copy(model.solver.solution)
   
   Solution!(model, param_values_ω_B)
   mat = hcat(model.solver.solution.x, model.solver.solution.θ, model.solver.solution.μ)
   # eval_wages = eval( build_function( substitute( expand_derivatives( Dθ( model.f ) ), param_values_ω_B ) , [x,θ,μ] ) )
   # model.solver.solution.w  = map(i -> eval_wages(mat), 1:model.solver.nk)
   solution_plot.ω_B[param_var[:ω_B][i]] = copy(model.solver.solution)

end


println("Plotting Results")

solver_vars = model.solver.vars
labels = Dict(:ω_A => "\\omega_A", :ω_B => "\\omega_B", :σ_A => "\\sigma_A")

for j in 1:3
   field = fieldnames(Solutions_for_Plots)[j]

   s = getfield( solution_plot, field)

   p_μ = plot()
   p_θ_1 = plot()
   p_θ_2 = plot()
   p_w = plot()
   
   lab = labels[field]
   scal = 2
   for i∈1:3
      println("Plotting: " , field , " = " , x)
      x = param_var[field][i]
      
      # Plot μ, θ(x), θ(μ(x))
      legend_loc =  (field == :σ_A) ? :topleft : :topright
      p_μ = plot!(p_μ, s[x].x, s[x].μ, title=L"\mu(x)", xlabel = L"x", label="",  thickness_scaling = scal)
      p_θ_1 = plot!(p_θ_1, s[x].x, s[x].θ, title=L"\theta(x)", xlabel = L"x", label=L"%$(lab) = %$x",  thickness_scaling = scal)
      ylims!(0,400 - 150 * (field == :σ_A))
      p_θ_2 = plot!(p_θ_2, s[x].μ, s[x].θ, title=L"\theta(\mu(x))", xlabel = L"\mu(x)", label="",  thickness_scaling = scal)
      ylims!(0,400 - 150 * (field == :σ_A))

      
      p_μ.series_list[end][:linecolor] = "#d1244f"
      p_θ_1.series_list[end][:linecolor] = "#d1244f"
      p_θ_2.series_list[end][:linecolor] = "#d1244f"
      p_μ.series_list[end][:linewidth] = 3
      p_θ_1.series_list[end][:linewidth] = 3
      p_θ_2.series_list[end][:linewidth] = 3
      # p_θ_1.series_list[end][:linecolor] = "#d1244f"
      # p_θ_2.series_list[end][:linecolor] = "#d1244f"
      # p_w = plot(matrix_sol[:, 1], matrix_sol[:,3], title="w(x)", xlabel = "x", label="")
      # p_Π = plot(matrix_sol[:, 1], matrix_sol[:,3], title="Π(x)", xlabel = "x", label="")
      plot(p_μ, p_θ_1, p_θ_2, layout=(1,3), size=(1200,450), legend = legend_loc)  
      savefig("./Eeckhout_Kircher-Econometrica(2018)/document/figures/plot_$(assortativity)_$(field)_$(i).pdf")
      p_μ.series_list[end][:linealpha] = 0.4
      p_θ_1.series_list[end][:linealpha] = 0.4
      p_θ_2.series_list[end][:linealpha] = 0.4
      p_μ.series_list[end][:linecolor] = "grey"
      p_θ_1.series_list[end][:linecolor] = "grey"
      p_θ_2.series_list[end][:linecolor] = "grey"
      p_μ.series_list[end][:linewidth] = 1
      p_θ_1.series_list[end][:linewidth] = 1
      p_θ_2.series_list[end][:linewidth] = 1
      
      # Plot w(x)
      p_w = plot!(p_w, s[x].x, s[x].w, title=L"w(x)", xlabel = L"x",  label=L"%$(lab) = %$x" )
   end
   
   colors = ["#66C2A5", "#FC8D62", "#8DA0CB"]

   for i ∈ 1:3
      p_μ.series_list[i][:linealpha] = 0.9
      p_θ_1.series_list[i][:linealpha] = 0.9
      p_θ_2.series_list[i][:linealpha] = 0.9
      p_w.series_list[i][:linealpha] = 0.9
      p_μ.series_list[i][:linecolor] = colors[i]
      p_θ_1.series_list[i][:linecolor] = colors[i]
      p_θ_2.series_list[i][:linecolor] = colors[i]
      p_w.series_list[i][:linecolor] = colors[i]
      p_μ.series_list[i][:linewidth] = 1
      p_θ_1.series_list[i][:linewidth] = 1
      p_θ_2.series_list[i][:linewidth] = 1
      p_w.series_list[i][:linewidth] = 1
   end

   
   plot(p_μ, p_θ_1, p_θ_2, layout=(1,3), size=(1200,420), thickness_scaling = 1)  
   
   savefig("./Eeckhout_Kircher-Econometrica(2018)/document/figures/plot_$(assortativity)_$(field).pdf")
   
   plot(p_w, size = (600, 600))
   savefig("./Eeckhout_Kircher-Econometrica(2018)/document/figures/plot_$(assortativity)_$(field)_wages.pdf")

end  