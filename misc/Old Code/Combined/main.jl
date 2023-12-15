include("assembly.jl")
include("fileio.jl")

using FiniteDifferences, FiniteDiff
using Printf

nu = 0.2

quad = getQuadratureRules()
shape = getShapeAndDerivatives(quad)
ncon, pos, boundaries, boundaryNodes = importMsh("meshes/SimpleBar.msh")


function computeCompliance(x)
	KS,fS,KH,fH = globalAssembly(x, ncon, pos, quad, shape)
	handleDirichletBCsSolid!(boundaryNodes["Left"], both, 0.0, KS, fS)
	handleNeumannBCsSolid!(ncon[boundaries["Right"]], pos, [0.0, 10.0], fS)

	handleDirichletBCsHeat!(boundaryNodes["Top"], 0.0, KH, fH)
	handleDirichletBCsHeat!(boundaryNodes["Left"], 0.0, KH, fH)
	handleDirichletBCsHeat!(boundaryNodes["Right"], 0.0, KH, fH)
	handleDirichletBCsHeat!(boundaryNodes["Bottom"], 0.0, KH, fH)

	d = KS \ fS
	T = KH \ fH
	return d' * fS, T' * fH
end;

function computeComplianceSolid(x)
	KS,fS = globalAssembly(x, ncon, pos, quad, shape, "solid")
	handleDirichletBCsSolid!(boundaryNodes["Left"], both, 0.0, KS, fS)
	handleNeumannBCsSolid!(ncon[boundaries["Right"]], pos, [0.0, 10.0], fS)
	d = KS \ fS
	return d' * fS
end

function computeComplianceHeat(x)
	_,_,KH,fH = globalAssembly(x, ncon, pos, quad, shape, "heat")
	handleDirichletBCsHeat!(boundaryNodes["Top"], 0.0, KH, fH)
	handleDirichletBCsHeat!(boundaryNodes["Left"], 0.0, KH, fH)
	handleDirichletBCsHeat!(boundaryNodes["Right"], 0.0, KH, fH)
	handleDirichletBCsHeat!(boundaryNodes["Bottom"], 0.0, KH, fH)
	T = KH \ fH
	return T' * fH
end;

function solve(x)
	KS,fS,KH,fH = globalAssembly(x, ncon, pos, quad, shape)
	handleDirichletBCsSolid!(boundaryNodes["Left"], both, 0.0, KS, fS)
	handleNeumannBCsSolid!(ncon[boundaries["Right"]], pos, [0.0, 10.0], fS)

	handleDirichletBCsHeat!(boundaryNodes["Top"], 0.0, KH, fH)
	handleDirichletBCsHeat!(boundaryNodes["Left"], 0.0, KH, fH)
	handleDirichletBCsHeat!(boundaryNodes["Right"], 0.0, KH, fH)
	handleDirichletBCsHeat!(boundaryNodes["Bottom"], 0.0, KH, fH)

	d = KS \ fS
	T = KH \ fH
	export2vtk("combinedDisp.vtu", pos, ncon, d)
	export2vtk("combinedTemp.vtu", pos, ncon, T)

	return;
end;

function generateParetoFront(nPoints = 10)
	# Get minimum of individual problems
	println("Finding minimum of heat problem")
	xHeat = findOptimum(computeComplianceHeat, [5.0,6.0], 1e-4, 0.75, false)
	cHeatMin = computeComplianceHeat(xHeat)
	cSolidMax = computeComplianceSolid(xHeat)
	Printf.@printf "Optimum of heat problem at [%.4f %.4f]\n" xHeat[1] xHeat[2]

	println("Finding minimum of solid problem")
	xSolid = findOptimum(computeComplianceSolid, [2.0,9.0], 1e-4, 0.75, false)
	cSolidMin = computeComplianceSolid(xSolid)
	cHeatMax = computeComplianceHeat(xSolid)
	Printf.@printf "Optimum of solid problem at [%.4f %.4f]\n" xSolid[1] xSolid[2]

	Printf.@printf "Range of heat problem = [%.4f %.4f] - Range of solid problem = [%.4f %.4f]\n" cHeatMin cHeatMax cSolidMin cSolidMax

	# Constrained Optimisation
	paretoFront = zeros(nPoints, 2)
	constraintArray = LinRange(cSolidMin, cSolidMax, nPoints);
	guess = [xSolid.-0.1; -1.0]

	println("Generating pareto front")
	paretoFront[1,:] = [cHeatMax, cSolidMin]
	paretoFront[end,:] = [cHeatMin, cSolidMax]
	for i = 2:nPoints-1
		c = constraintArray[i]
		L(x) = computeComplianceHeat(x[1:2]) - x[3] * (computeComplianceSolid(x[1:2]) - c)
		guess = newtonMethod(L, guess, 1e-4, false)
		paretoFront[i,:] = [computeComplianceHeat(guess[1:2]) computeComplianceSolid(guess[1:2])] 
		
		Printf.@printf "Pair %i: %.4f, %.4f\n\n" i paretoFront[i,1] paretoFront[i,2]
	end

	return paretoFront;
end;

# Find minimum using gradient decent
function findOptimum(funcptr, guess, tol, alpha = 0.1, suppress::Bool = true)
	while true
		g = FiniteDifferences.grad(FiniteDifferences.central_fdm(3,1), funcptr, guess)[1]
		if !suppress
			@printf "Value: %.4f, %.4f - Gradient: %.4f, %.4f - α: %.2e - Tolerance: %.2e\n" guess[1] guess[2] g[1] g[2] alpha tol
		end

		if abs(g[1]) < tol && abs(g[2]) < tol
			break;
		else
			guess = guess .- alpha .* g
		end
	end

	return guess
end;

# Find optimum using newton-raphson method (for finding saddle points)
function newtonMethod(funcptr, guess, tol, suppress::Bool = true)
	while true
		g = FiniteDifferences.grad(FiniteDifferences.central_fdm(3,1), funcptr, guess)[1]
		h = FiniteDiff.finite_difference_hessian(funcptr, guess)
		update = h \ g
		if !suppress
			@printf "Value: %.4f, %.4f, %.4f - Update: %.4f, %.4f, %.4f\n" guess[1] guess[2] guess[3] update[1] update[2] update[3]
		end

		if all(abs.(update) .< tol)
			break;
		else
			guess = guess .- update
		end
	end

	return guess
end;