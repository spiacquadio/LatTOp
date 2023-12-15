include("assembly.jl")
include("processing.jl")
include("fileio.jl")

using FiniteDifferences
using Printf

nu = 0.2

quad = getQuadratureRules()
shape = getShapeAndDerivatives(quad)

ncon, pos, boundaries, boundaryNodes = importMsh("meshes/CombinedDomain.msh")
K,f = getSparsityPattern(ncon)
Ke = zeros(8,8,length(boundaries["Domain"]))
HF = getWeightFactorMatrix(ncon[boundaries["Domain"]], pos, 0.005)

# Topology Optimization
function optimiseTopology(E, volFrac)
	itr = 0; p = 1; pmax = 3; c = 0.0; u = []; change = 1.0
	
	nElements = size(ncon[boundaries["Domain"]])[1]
	relDense = volFrac * ones(Float64, nElements)
	relDenseNew = 0.0 * ones(Float64, nElements)
	sensitivity = Vector{Float64}(undef, nElements)
	while (itr < 100 && change > 1e-4)
		# Main FE solve
		assembleTime = @elapsed begin
			globalAssembly!(E, relDense.^p, ncon[boundaries["Domain"]], pos, Ke, K, f, quad, shape)
			handleDirichletBCs!(boundaryNodes["Left"], both, 0.0, K, f)
			handleNeumannBCs!(ncon[boundaries["Traction"]], pos, [0.0, -20.0], f)
		end
		solveTime = @elapsed begin
			u = K \ f	# Solve for temperature
		end

		# Compute element sensitivities (Computes ∂c/∂e)
		sensitivityTime = @elapsed begin
			if itr > 10 && p < pmax
				p *= 1.02
			end
			tmp = -(p .* relDense.^(p - 1))
			Threads.@threads for e = 1:nElements
				uElement = zeros(8)
				for i = 1:4
					uElement[2*i] = u[2*ncon[boundaries["Domain"]][e][i]]
					uElement[2*i-1] = u[2*ncon[boundaries["Domain"]][e][i]-1]
				end
				sensitivity[e] = tmp[e] * uElement' * Ke[:,:,e] * uElement
			end
			sensitivity = HF * sensitivity
		end

		# Bisection method to find optimum Lagrange Multiplier
		bisectionTime = @elapsed begin
			lambda1 = 0.0; lambda2 = 1.0e9; move = 0.2
			while (lambda2 - lambda1) / (lambda1 + lambda2) > 1e-4
				mid = 0.5 * (lambda1 + lambda2)
				
				B = -sensitivity / mid
				xb = relDense .* sqrt.(B)

				relDenseNew = max.(0.01, max.(relDense.-move), min.(1.0, min.(relDense.+move, xb)))

				if sum(relDenseNew) - (volFrac * nElements) > 0
					lambda1 = mid
				else
					lambda2 = mid
				end
			end
		end

		# Iterate
		itr += 1
		c = u' * f	# Compute compliance
		change = abs(maximum(relDenseNew - relDense))
		relDense = relDenseNew

		Printf.@printf "Itr: %i - Conv: %.2e - Comp: %.2e [TIMERS] assembly = %.4f, solve = %.4f, sens = %.4f, bisection = %.4f\n" itr change c assembleTime solveTime sensitivityTime bisectionTime
	end

	export2vtk("Results/test-opt_solid.vtu", pos, ncon, u, relDense)
	return relDense, c
end;

# Simple Solve
function solve(E)
	nElements = size(ncon[boundaries["Domain"]])[1]
	relDense = ones(Float64, nElements)
	globalAssembly!(E, relDense, ncon[boundaries["Domain"]], pos, Ke, K, f, quad, shape)
	handleDirichletBCs!(boundaryNodes["Left"], both, 0.0, K, f)
	handleNeumannBCs!(ncon[boundaries["Traction"]], pos, [0.0, -20.0], f)
	u = K \ f
	export2vtk("Results/test_solid.vtu", pos, ncon, u)
	return u
end;
