include("assembly.jl")
include("generateGrid.jl")
include("fileio.jl")

using FiniteDifferences
using Printf

quad = getQuadratureRules(Quad1)
shape = getShapeAndDerivatives(Quad1, quad)

ncon, pos, boundaries, boundaryNodes = importMsh("meshes/SimpleBar.msh")

function computeCompliance(k)
	K,f = globalAssembly(k, ncon, pos, quad, shape)
	handleDirichletBCs!(boundaryNodes["Top"], 0.0, K, f)
	handleDirichletBCs!(boundaryNodes["Left"], 0.0, K, f)
	handleDirichletBCs!(boundaryNodes["Right"], 0.0, K, f)
	handleDirichletBCs!(boundaryNodes["Bottom"], 0.0, K, f)
	
	u = K \ f
	return u' * f
end


# Find minimum
function findOptimum(guess, tol, alpha = 0.1)
	while true
		g = FiniteDifferences.grad(FiniteDifferences.central_fdm(3,1), computeCompliance, guess)[1]
		@printf "Value: %.4f, %.4f - Gradient: %.4f, %.4f - α: %.2e - Tolerance: %.2e\n" guess[1] guess[2] g[1] g[2] alpha tol

		if abs(g[1]) < tol && abs(g[2]) < tol
			break;
		else
			guess = guess .- alpha .* g
		end
	end

	return guess
end;