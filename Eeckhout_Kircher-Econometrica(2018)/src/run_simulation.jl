include("simulation_functions.jl")

model = initialize_model()

ω_A, ω_B, σ_A = model.params

param_values = Dict(ω_A => 0.25, ω_B => 0.5, σ_A => 0.5)


model.solver.nk

Solution!(model::Model, param_values)

sys_vars = model.solver.vars

model.solver.initial_condition[sys_vars[:θ]] = 34.0

model.solver.initial_condition

sol = model.solver.solution

p_μ = plot(matrix_sol[:, 1], matrix_sol[:,3], title="μ(x)", xlabel = "x", label="")
	ylims!(0,200)
	plot!([0, 100], [100, 100], label="100")
p_θ = plot(matrix_sol[:, 1], matrix_sol[:,2], title="θ(x)", xlabel = "x", label="")
p_w = plot(matrix_sol[:, 1], matrix_sol[:,3], title="w(x)", xlabel = "x", label="")
p_Π = plot(matrix_sol[:, 1], matrix_sol[:,3], title="Π(x)", xlabel = "x", label="")
plots = [p_μ, p_θ, p_w, p_Π]