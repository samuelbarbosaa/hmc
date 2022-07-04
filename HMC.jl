### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 8e871532-e991-4c21-997f-9228d007153f
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.add("ForwardDiff")
	Pkg.add("Plots")
end;

# ╔═╡ 0047b944-96da-4515-9c9f-cd962c034337
# ╠═╡ show_logs = false
begin
	using Random, Distributions, LinearAlgebra, ForwardDiff, Plots;
	Random.seed!(123)
end;

# ╔═╡ 8da87152-fa10-11ec-0f98-557dd3908a29
md"""
# Hamiltonian Monte Carlo
"""

# ╔═╡ 153713bb-a384-4af7-87f2-ac3c9ae6c09b
md"""
Função $\pi(\theta) = \exp[\frac{1}{2} (\theta_1^2 \theta_2^2 + \theta_1^2 + \theta_2^2 - 8\theta_1 -8\theta_2) ]$:
"""

# ╔═╡ 9e0308ca-317b-4e48-ab4d-f0187581ef47
function π(θ)
	exp(-1/2 *  (prod(θ)^2 + sum(θ.^2) - 8 * sum(θ)) )
end

# ╔═╡ bc6289bd-1aaa-4537-a5ac-a1789a17bae9
begin
	x = -0.5:0.01:6.2;
	y = -0.5:0.01:6.2;
	contour(x, y, (x,y) -> π((x,y)))
end

# ╔═╡ 679b6cb2-72e4-4fda-8928-2a207f817f59
begin
	Z = zeros(size(x)[1], size(y)[1])
	for (i, u) in enumerate(x)
		for (j, v) in enumerate(y)
			Z[i,j] = π([u,v])
		end
	end
	plot(x, sum(Z, dims=2) / sum(Z))
end

# ╔═╡ afb8e8e0-e2f3-4c1b-816f-d7c67c5932b9
function U(θ)
	-log(π(θ))
end

# ╔═╡ df158b2f-7c19-4aa5-906b-b3e58e7e263b
function dU(U, θ)
	ForwardDiff.gradient(U, θ)
end

# ╔═╡ f64603ce-4965-460a-89b5-4f46daa52934
function K(ρ, M)
	1/2 * ρ' * M^-1 * ρ
end

# ╔═╡ 1204af10-2d1f-4f8f-9474-bdd2e09180d3
function leapfrog_step(θ, ρ, U, dU, M, L, ϵ)
	ρ = ρ + ϵ/2 * dU(U, θ)
	θ = θ + ϵ * M^-1 * ρ
	ρ = ρ + ϵ/2 * dU(U, θ)
	return(θ, ρ)
end

# ╔═╡ df203af8-d84b-4a45-a894-4b7d11d20a2e
function metropolis_step(θ, ρ, θₜ, ρₜ, U, K, M)
	r = exp(U(θ) - U(θₜ) + K(ρ, M) - K(ρₜ, M))
	u = rand(Uniform())
	if u < r
		θ = θₜ
		ρ = -ρₜ
	end
	return(θ, ρ)
end

# ╔═╡ efad6391-13f4-4b7f-b306-a5b3bfb23742
function HMC(θ, U, dU, K, M, L, ϵ, num_iter)
	d = size(θ, 1)
	Θ = zeros(num_iter, d)
	Ρ = zeros(num_iter, d)
	
	for t = 1:num_iter
		
		# Amostra ρ
		Σ = cholesky(M).L
		P = MvNormal(Σ)
		ρ = rand(P)
		
		# Inicializa estado no passo t
		ρₜ = ρ
		θₜ = θ
		
		# Leapfrog steps
		for l = 1:L
			θₜ, ρₜ = leapfrog_step(θₜ, ρₜ, U, dU, M, L, ϵ)
		end
		
		# Metropolis step
		θ, ρ = metropolis_step(θ, ρ, θₜ, ρₜ, U, K, M)

		# Registra valores amostrados
		Θ[t,:] = θ
		Ρ[t,:] = ρ
		
	end
	return((Θ, Ρ))
end

# ╔═╡ 4f25207a-2f1d-47d2-8798-5fbe104aecda
begin
	d = 2
	θ₀ = rand(d)
	M = I(d)
	L = 10
	ϵ = 0.01
	burn_in = 10000
	N = 10000
	num_iter = burn_in + N
end;

# ╔═╡ 6b6cef15-d0bc-4a00-b82a-60e4875caad3
(A, B) = HMC(θ₀, U, dU, K, M, L, ϵ, num_iter);

# ╔═╡ c3657dc4-3f74-4efb-9bc0-9bbb69894bd4
histogram(A[:, 1])

# ╔═╡ f1fe72ee-76bf-40c3-9fd1-1f0e10426dcb
plot(1:num_iter, A[:, 1])

# ╔═╡ Cell order:
# ╟─8da87152-fa10-11ec-0f98-557dd3908a29
# ╠═8e871532-e991-4c21-997f-9228d007153f
# ╠═0047b944-96da-4515-9c9f-cd962c034337
# ╟─153713bb-a384-4af7-87f2-ac3c9ae6c09b
# ╠═9e0308ca-317b-4e48-ab4d-f0187581ef47
# ╠═bc6289bd-1aaa-4537-a5ac-a1789a17bae9
# ╠═679b6cb2-72e4-4fda-8928-2a207f817f59
# ╠═afb8e8e0-e2f3-4c1b-816f-d7c67c5932b9
# ╠═df158b2f-7c19-4aa5-906b-b3e58e7e263b
# ╠═f64603ce-4965-460a-89b5-4f46daa52934
# ╠═1204af10-2d1f-4f8f-9474-bdd2e09180d3
# ╠═df203af8-d84b-4a45-a894-4b7d11d20a2e
# ╠═efad6391-13f4-4b7f-b306-a5b3bfb23742
# ╠═4f25207a-2f1d-47d2-8798-5fbe104aecda
# ╠═6b6cef15-d0bc-4a00-b82a-60e4875caad3
# ╠═c3657dc4-3f74-4efb-9bc0-9bbb69894bd4
# ╠═f1fe72ee-76bf-40c3-9fd1-1f0e10426dcb
