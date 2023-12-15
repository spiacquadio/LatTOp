include("assembly.jl")
include("generateGrid.jl")
include("processing.jl")
include("fileio.jl")

using FiniteDifferences
using Printf

nu = 0.2

quad = getQuadratureRules()
shape = getShapeAndDerivatives(quad)
ncon, pos, boundaries, boundaryNodes = importMsh("meshes/SimpleBar.msh")


function computeCompliance(E)
	K,f = globalAssembly(E, ncon, pos, quad, shape)
	handleDirichletBCs!(boundaryNodes["Left"], both, 0.0, K, f)
	handleNeumannBCs!(ncon[boundaries["Right"]], pos, [0.0, 10.0], f)
	u = K \ f
	return u' * f
end

function solve(E)
	K,f = globalAssembly(E, ncon, pos, quad, shape)
	handleDirichletBCs!(boundaryNodes["Left"], both, 0.0, K, f)
	handleNeumannBCs!(ncon[boundaries["Right"]], pos, [0.0, 10.0], f)
	u = K \ f
	export2vtk("test.vtu", pos, ncon, u)
	return u
end


#export2vtk("test.vtu", pos, ncon, u)

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