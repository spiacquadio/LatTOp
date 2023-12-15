include("assembly.jl")
include("fileio.jl")

using FiniteDifferences
using Printf
using Plots

# Setup
quad = getQuadratureRules(Quad1)			# Get quadrature rules for 1st order quadrelaterals
shape = getShapeAndDerivatives(Quad1, quad)	# Get shape functions for 1st order quads

# Import mesh. (ncon = node connectivity, pos = nodal positions, boundaries = maps boundary name to boundary elements, boundaryNodes = maps boundary name to boundary nodes)
ncon, pos, boundaries, boundaryNodes = importMsh("meshes/ThermalDomain.msh")

# Initialise domain
K,f = getSparsityPattern(ncon)					# Gets sparcity pattern of global matrix and forcing array
Ke = zeros(4,4,length(boundaries["Domain"]))	# Initialise element matrices (used in sensitivity analysis)
HF = getWeightFactorMatrix(ncon[boundaries["Domain"]], pos, 0.005)	# Weight factor matrix (filter matrix to filter sensitivities)
neighbours = getNeighbours(ncon[boundaries["Domain"]])				# Get element neighbours


#=
	Topology Optimization function
	Optimises topology using the SIMP method

	This function might need to be optimised because the bisection part can take a 
		significant amount of time. The solve time is relatively low. But it seems like
		Julia doesnt like for loops (kinda like MATLAB) so the assembly time is quite high.
		However, I dont really want to use a library (yet) because by doing it manually,
		i can save all the information that is needed, and this could potentially be expanded
		to use an AD library with more ease. But we can discuss this at a later time.

	Inputs:
		- k (temporary variable which will be used for porosity of cell)
		- volFrac (desired volume fraction)

	Outputs:
		- relDense (array of element relative densities)
		- c (domain compliance)
=#
function optimiseTopology(k, volFrac)
	itr = 0; p = 1; pmax = 4; c = 0.0; u = []; change = 1.0
	nElements = size(ncon[boundaries["Domain"]])[1]

	# Preallocate some arrays
	relDense = volFrac * ones(Float64, nElements)
	relDenseNew = 0.0 * ones(Float64, nElements)
	sensitivity = Vector{Float64}(undef, nElements)

	while (itr < 150 && change > 1e-4)
		# Main FE solve
		assembleTime = @elapsed begin
			globalAssembly!(k, relDense, p, ncon[boundaries["Domain"]], pos, Ke, K, f, quad, shape)
			handleNeumannBCs!(ncon[boundaries["Heat"]], pos, 100.0, f)
		end
		solveTime = @elapsed begin
			u = K \ f	# Solve for temperature
		end

		# Compute element sensitivities (Computes ∂c/∂e)
		sensitivityTime = @elapsed begin
			# Increase penalisation parameter after 10 iterations to improve convergence
			if itr > 10 && p < pmax
				p *= 1.02
			end

			# Note objective function is min: c = T' * f; s.t: ∫ρ_e = γ * f (where ρ_e is the element density, γ is the volume fraction)
			# Compute element sensitivities for each element. ∂c/∂ρ = p * ρ^(p-1) * u_e' * Ke_e * u_e
			Threads.@threads for e = 1:nElements
				sensitivity[e] = -(p * relDense[e]^(p - 1)) * u[ncon[boundaries["Domain"]][e]]' * Ke[:,:,e] * u[ncon[boundaries["Domain"]][e]]
			end
			sensitivity = HF * sensitivity	# Filter sensitivities
		end

		# Bisection method to find optimum Lagrange Multiplier
		# This is essentially a basic root finder
		# Refer to the paper (https://doi.org/10.3390/app11073175) for equation numbers
		# Currently this uses the 'non-generalised' version as I only have 1 constraint
		bisectionTime = @elapsed begin
			lambda1 = 0.0; lambda2 = 1.0e9; 
			move = 0.2	# This parameter controls how much an element density can change. Too large and a local optimum is found too quickly which isnt close to the global optimum
			while (lambda2 - lambda1) / (lambda1 + lambda2) > 1e-4	# While root not found
				mid = 0.5 * (lambda1 + lambda2)	# Bisect
				B = -sensitivity / mid			# Equation 6
				xb = relDense .* sqrt.(B)		# Equation 7 (part 1)

				# Equation 7 (fully implimented)
				relDenseNew = max.(0.0, max.(relDense.-move), min.(1.0, min.(relDense.+move, xb)))

				# This checks the ∂L/∂λ = ∫ρ_e - γ = 0 condition and updates the bisection accordingly
				if sum(relDenseNew) - (volFrac * nElements) > 0
					lambda1 = mid
				else
					lambda2 = mid
				end
			end
		end

		# Iterate
		itr += 1
		c = u' * f										# Compute compliance
		change = abs(maximum(relDenseNew - relDense))	# Compute the maximum change in relative densities
		relDense = relDenseNew							# Update relDense
		
		Printf.@printf "Itr: %i - Conv: %.2e - Comp: %.2e [TIMERS] assemble = %.4f, solve = %.4f, sens = %.4f, bisection = %.4f\n" itr change c assembleTime solveTime sensitivityTime bisectionTime
	end

	export2vtk("Results/test-opt.vtu", pos, ncon, u, relDense)
	return relDense, c
end;


# Simple solve
function solve(k)
	relDense = ones(length(ncon[boundaries["Domain"]]))
	globalAssembly!(k, relDense, 1, ncon[boundaries["Domain"]], pos, Ke, K, f, quad, shape)
	handleDirichletBCs!(boundaryNodes["Traction"], 0.0, K, f)
	handleNeumannBCs!(ncon[boundaries["Heat"]], pos, 10.0, f)
	u = K \ f	# Solve for temperature

	export2vtk("Results/test.vtu", pos, ncon, u)
end;