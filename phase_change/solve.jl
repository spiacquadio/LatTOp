include("assembly.jl")
include("../src/fileio.jl")
using Printf
using SparseArrays
using ReverseDiff

## Main solve function
function solve(endtime, dt, initialTemp, itrsPerSave = 1)
	temperature = initialTemp .* ones(nNodes);
	time = 0.0; saveCounter = 1;

	# Save initial conditions and init pvd
	pvdpath = "res/results_pc/test_Paraffin_$dt.pvd";
	initExport2pvd(pvdpath, pos, ncon, time, temperature);

	# Assemble constant matrices
	K,_ = getSparsityPatternHeat(ncon[boundaries["Domain"]]);
	C,_ = getSparsityPatternHeat(ncon[boundaries["Domain"]]);
	F = Array{Float64}(undef, nNodes);
	globalAssembly!(ncon[boundaries["Domain"]], pos, K,C,F, quad, shape);
	handleNeumannBCs!(ncon[boundaries["Heat"]], pos, 10000.0, F);

	# Main transient iteration loop
	while time < endtime
		temperature = solveTimestep(dt, temperature, K,C,F);
		time += dt;
		Printf.@printf "Time: %.2f, Temperature Range (w.r.t T_m): %.2f - %.2f\n" time minimum(temperature)-tm maximum(temperature)-tm

		# Save timestep
		if saveCounter == itrsPerSave || time >= endtime
			export2pvd(pvdpath, pos, ncon, time, temperature);
			saveCounter = 0;
		end

		saveCounter += 1;
	end

	return temperature;
end

## Solves a single timestep
function solveTimestep(dt, prevTemp, K,C,F)
	global newTemp = prevTemp;
	res = 1.0; tol = 1e-4; maxItrs = 10000; itr = 0;

	# Solve linear system
	if maximum(prevTemp) < tm-tr || minimum(prevTemp) > tm+tr
		residual = (F .* dt) .+ (C * prevTemp);
		jacobian = (K .* dt) .+ C;
		newTemp = sparse(jacobian) \ residual;

		# Ensure solution is in linear region
		if maximum(newTemp) < tm - tr || minimum(newTemp) > tm+tr
			return newTemp
		end
	end
	
	# Preallocate matrices
	dL, L = getSparsityPatternHeat(ncon[boundaries["Domain"]]);
	Lprev = getLatentHeatVector(ncon[boundaries["Domain"]], elements, prevTemp);

	# Precompute vector and matrix
	global vecPart = (F .* dt) .+ (C * prevTemp) .+ Lprev;
	global matPart = (C .+ K .* dt);
	while res > tol && itr < maxItrs 
		iterateLatentHeat!(ncon[boundaries["Domain"]], elements, dL, L, newTemp);
		residual = (vecPart .- L) .- (matPart) * newTemp;
		jacobian = matPart .+ dL;
		global direction = jacobian \ residual;

		stepsize, m = linesearchBacktrack(newTemp, 0.9, direction);
		newTemp =  newTemp .+ stepsize .* direction;
		res = sqrt(residualNorm2(newTemp));
		itr += 1;
		Printf.@printf "Newton Iteration %i: res = %.2e, α = %.2f, m = %.2e\n" itr res stepsize m
	end;
	return newTemp;
end;

function residualNorm2(x)
	L = getLatentHeatVector(ncon[boundaries["Domain"]], elements, x);
	residual = (vecPart .- L) .- (matPart) * x;
	return residual' * residual;
end;

function linesearchBacktrack(x, stepsize, direction)
	decay = 0.8;	# decays the stepsize
	beta = 0.01;	# Used for Armijo condition
	diff = DiffResults.DiffResult(0.0, Vector{Float64}(undef, nNodes));	# ReverseDiff diffresult object

	itr = 0;
	maxItrs = 4;
	
	m = 0;
	ReverseDiff.gradient!(diff, residualNorm2, x);
	fx = diff.value;
	gradfx = DiffResults.gradient(diff);
	m = direction' * gradfx;
	while itr < maxItrs
		fs = residualNorm2(x + stepsize.*direction);
		if fs > fx + beta * stepsize * m
			stepsize *= decay;
		else
			break;
		end
		itr += 1;
	end

	return stepsize, m;
end;

function dampedNewton(x)

end;