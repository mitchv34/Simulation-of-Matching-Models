using LaTeXStrings
include("simulation_functions.jl")

theme(:vibrant) 
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

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

param_var = [0.5, 0.25, 0.75] # parameter variation

for val ∈ param_var
   
   param_values_ω_A = Dict(ω_A => val, ω_B => 0.5, σ_A => 0.9)
   param_values_σ_A = Dict(ω_A => 0.3, ω_B => 0.6, σ_A => val)
   param_values_ω_B = Dict(ω_A => 0.5, ω_B => val, σ_A => 0.9)
   # pack parameter values into a dictionary
   
   Solution!(model, param_values_ω_A)
   solution_plot.ω_A[val] = model.solver.sol

   Solution!(model, param_values_σ_A)
   solution_plot.σ_A[val] = model.solver.sol

   Solution!(model, param_values_ω_B)
   solution_plot.ω_B[val] = model.solver.sol

end

solver_vars = model.solver.vars

labels = Dict(:ω_A => "\\omega_A", :ω_B => "\\omega_B", :σ_A => "\\sigma_A")
for j in 1:3
   field = fieldnames(Solutions_for_Plots)[j]

   s = getfield( solution_plot, field)

   p_μ = plot()
   p_θ_1 = plot()
   p_θ_2 = plot()
   
   lab = labels[field]
   scal = 2
   for i∈1:3

      x = param_var[i]
      p_μ = plot!(p_μ, s[x].t, s[x][solver_vars[:μ]], title=L"\mu(x)", xlabel = L"x", label="",  thickness_scaling = scal)
      p_θ_1 = plot!(p_θ_1, s[x].t, s[x][solver_vars[:θ]], title=L"\theta(x)", xlabel = L"x", label=L"%$(lab) = %$x",  thickness_scaling = scal)
      ylims!(0,400)
      p_θ_2 = plot!(p_θ_2, s[x][solver_vars[:μ]], s[x][solver_vars[:θ]], title=L"\theta(\mu(x))", xlabel = L"\mu(x)", label="",  thickness_scaling = scal)
      ylims!(0,400)
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
      plot(p_μ, p_θ_1, p_θ_2, layout=(1,3), size=(1200,450))  
      savefig("Eeckhout_Kircher-Econometrica(2018)/document/figures/plot_$(assortativity)_$(field)_$(i).pdf")
      p_μ.series_list[end][:linealpha] = 0.4
      p_θ_1.series_list[end][:linealpha] = 0.4
      p_θ_2.series_list[end][:linealpha] = 0.4
      p_μ.series_list[end][:linecolor] = "grey"
      p_θ_1.series_list[end][:linecolor] = "grey"
      p_θ_2.series_list[end][:linecolor] = "grey"
      p_μ.series_list[end][:linewidth] = 1
      p_θ_1.series_list[end][:linewidth] = 1
      p_θ_2.series_list[end][:linewidth] = 1
   end

   colors = ["#66C2A5", "#FC8D62", "#8DA0CB"]

   for i ∈ 1:3
      p_μ.series_list[i][:linealpha] = 0.9
      p_θ_1.series_list[i][:linealpha] = 0.9
      p_θ_2.series_list[i][:linealpha] = 0.9
      p_μ.series_list[i][:linecolor] = colors[i]
      p_θ_1.series_list[i][:linecolor] = colors[i]
      p_θ_2.series_list[i][:linecolor] = colors[i]
      p_μ.series_list[end][:linewidth] = 2.5
      p_θ_1.series_list[end][:linewidth] = 2.5
      p_θ_2.series_list[end][:linewidth] = 2.5
   end

   
   plot(p_μ, p_θ_1, p_θ_2, layout=(1,3), size=(1200,420))  
   
   savefig("Eeckhout_Kircher-Econometrica(2018)/document/figures/plot_$(assortativity)_$(field).pdf")

end  