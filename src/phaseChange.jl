include("assembly.jl")
include("fileio.jl")
using Printf
using SparseArrays
using Ferrite
using ReverseDiff
using Distributed
## Main solve function
@everywhere function solve(endtime, dt, initialTemp, reldense, savename, itrsPerSave = 1)
	temperature = initialTemp .* ones(nNodes);
	global timeElapsed = 0.0; 
	saveCounter = 1; 
    global p = 1.0; global useLattice = true;

	# Save initial conditions and init pvd
	if itrsPerSave > 0
		pvdpath = "res/results_pc/";
		while handleOverwrite(pvdpath * savename * ".pvd") == true
			println("You are about to overwrite an existing file. Would you like to continue? (y/n)");
			answer = readline();
			if answer == "y"
				break;
			elseif answer == "n"
				println("Enter a new save name");
				savename = readline();
			else
				println("Unrecognised key");
			end
		end
		savepath = pvdpath * savename * ".pvd"
		initExport2pvd(savepath, pos, ncon, timeElapsed, temperature, reldense);
	end


	# Assemble constant matrices
	K,_ = getSparsityPatternHeat(ncon[boundaries["Domain"]]);
	C,_ = getSparsityPatternHeat(ncon[boundaries["Domain"]]);
	F = Array{Float64}(undef, nNodes);
	globalAssemblyPC!(ncon[boundaries["Domain"]], pos, reldense, K,C,F, quad, shape);
	thermalBCs!(K, F);

	# Main transient iteration loop
	deltaTemp = zeros(nNodes);
	while timeElapsed < endtime
		prevTemp = temperature;
		temperature = solveTimestep(dt, temperature, deltaTemp, reldense, K,C,F);
		deltaTemp = temperature .- prevTemp;
		timeElapsed += dt;
		Printf.@printf "Time: %.2f, Temperature Range (w.r.t T_m): %.2f - %.2f\n" timeElapsed minimum(temperature)-tm maximum(temperature)-tm

		# Save timestep
		if saveCounter == itrsPerSave || (timeElapsed >= endtime && itrsPerSave > 0)
			export2pvd(savepath, pos, ncon, timeElapsed, temperature, reldense);
			saveCounter = 0;
		end

		saveCounter += 1;
	end

	return temperature, temperature' * F;
end


## Main solve with gradients
@everywhere function solveWithGradient(endTime, dt, initialTemp, reldense)
	temperature = initialTemp .* ones(nNodes);
	global timeElapsed = 0.0; 
	saveCounter = 1; 
    global p = 1.0; global useLattice = true;

	K,_ = getSparsityPatternHeat(ncon[boundaries["Domain"]]);
	C,_ = getSparsityPatternHeat(ncon[boundaries["Domain"]]);
	F = Array{Float64}(undef, nNodes);
	globalAssemblyPC!(ncon[boundaries["Domain"]], pos, reldense, K,C,F, quad, shape);
	handleNeumannBCs!(ncon[boundaries["Heat"]], pos, 16000.0, F, 1);

	dKe = zeros(4,4, nElements);
	dCe = zeros(4,4, nElements);
	getDerivatives(dKe, dCe, reldense);

	timeElapsed = 0.0;
	dT = zeros(nNodes,nElements);
	temperature = initialTemp .* ones(nNodes);
	while timeElapsed < endTime
		temperature,dT = solveTimestepGradient(dt, temperature, reldense, K,C,F, dKe,dCe, dT);
		
		timeElapsed += dt;
		Printf.@printf "Time: %.2f, Temperature Range (w.r.t T_m): %.2f - %.2f\n" timeElapsed minimum(temperature)-tm maximum(temperature)-tm

		saveCounter += 1;
	end

	return temperature, temperature' * F, dT' * F;
end;

@everywhere function getDerivatives(dKe, dCe, reldense)
	@distributed for e = 1:nElements
		elePoints = pos[ncon[boundaries["Domain"][e]],:];
		_,_,_, dK, dC = assembleElementdPC(elePoints, reldense[e], quad, shape)

		dKe[:,:,e] .= dK;
		dCe[:,:,e] .= dC;
	end
end;

@everywhere function assembleElementdPC(element, reldense, quad, shape)
	Ke = zeros(typeof(reldense), 4,4);
	Ce = zeros(typeof(reldense), 4,4);
	dKe = zeros(typeof(reldense), 4,4);
	dCe = zeros(typeof(reldense), 4,4);
	fe = zeros(4);
	
	D = DiffResults.DiffResult(Matrix{Float64}(undef, 2,2), Matrix{Float64}(undef, 2,2));
	cp = DiffResults.DiffResult(0.0, [0.0]);
	rho = DiffResults.DiffResult(0.0, [0.0]);
	ForwardDiff.derivative!(D, getKmatrix, reldense);
	@distributed for q = 1:quad.nQuadPoints
		detJ, dNdX = getJacobian(q, shape, element);
		integrator = quad.weights[q] * detJ

		ForwardDiff.gradient!(cp, getCp, [reldense]);
		ForwardDiff.gradient!(rho, getRho, [reldense]);

		fe += zeros(4);
		Ke += integrator .* (dNdX' * D.value * dNdX);
		Ce += (integrator * rho.value * cp.value) .* (shape.N[q,:] * shape.N[q,:]');
		
		dKe += integrator .* (dNdX' * DiffResults.derivative(D) * dNdX);
		dCe += integrator * (rho.value * DiffResults.derivative(cp)[1] + 
			cp.value * DiffResults.derivative(rho)[1]) .* (shape.N[q,:] * shape.N[q,:]');
	end

	return Ke, Ce, fe, dKe, dCe;
end;

@everywhere function computeLinearGradient(ncon, A, C, b, dKe, dCe, dt, T, dTprev)
	dT = zeros(nNodes, nElements);
	dB = zeros(nNodes, nElements);
	invA = inv(Matrix(A)); x = A \ b;
	@distributed for i = 1:nElements
		conn = ncon[i];
		Asparse = dt.*dKe[:,:,i] .+ dCe[:,:,i];

		dT[:,i] = -invA[conn,:]' * (Asparse * x[conn]);
		dB[conn,i] = dCe[:,:,i] * T[conn];
	end

	return dT .+ (invA*(dB .+ C*dTprev));
end;

@everywhere function computeNonlinearGradient(reldense, ncon, A, J, C, r, dKe, dCe, dLPrev, dt, T, Tn, dTprev, dTi)
	dT = zeros(nNodes, nElements);
	dR = zeros(nNodes, nElements);
	invJ = inv(Matrix(J)); x = J \ r;

	dL,ddLe = computeLatentHeatGradient(ncon, elements, reldense, T, dTi, true);
	ddL,_ = getSparsityPatternHeat(ncon);
	assembler = start_assemble(ddL);
	@distributed for i = 1:nElements
		fill!(ddL, 0.0)
		@distributed for e = 1:nElements
			assemble!(assembler, ncon[e], ddLe[:,:,i,e]);
		end
		conn = ncon[i];
		dA = dt.*dKe[:,:,i] .+ dCe[:,:,i];
		assemble!(assembler, conn, dA);

		dT[:,i] = -invJ * (ddL * x);
		dR[conn,i] = dCe[:,:,i] * Tn[conn] - dA * T[conn];
	end
	
	dR += (C * dTprev - dL + dLPrev - A * dTi)
	return dT .+ (invJ * dR);
end;

@everywhere function computeLatentHeatGradient(ncon, elements, reldense, T, dT, matrix::Bool)
	dL = zeros(nNodes, nElements);
	ddLe = [];
	if matrix == true
		ddLe = zeros(4,4, nElements, nElements);
	end

	L = DiffResults.DiffResult(0.0, [0.0]);
	rho = DiffResults.DiffResult(0.0, [0.0]);
	
	@distributed for e = 1:nElements
		conn = ncon[e];

		ForwardDiff.gradient!(L, getL, [reldense[e]]);
		ForwardDiff.gradient!(rho, getRho, [reldense[e]]);
		@distributed for i = 1:4
			pf, dpf = getPhaseFraction(T[conn[i]]);
			ddpf = getdds(T[conn[i]]);
			dL[conn[i],e] += elements[e].fe[i] * (L.value*pf*DiffResults.derivative(rho)[1] +
				rho.value*pf*DiffResults.derivative(L)[1]);

			elements[e].fe[i]
			dL[conn[i],:] += elements[e].fe[i] * rho.value*L.value*dpf .* dT[conn[i],:];
			
			if matrix == true
				ddLe[i,:,e,e] += elements[e].Me[i,:] .* (L.value*dpf*DiffResults.derivative(rho)[1] + 
					rho.value*dpf*DiffResults.derivative(L)[1]);
					@distributed for j = 1:4
					ddLe[i,j,:,e] += elements[e].Me[i,j] * rho.value*L.value*ddpf .* dT[conn[i],:];
				end
			end
		end
	end
	return dL,ddLe;
end;

@everywhere function solveTimestepGradient(dt, prevTemp, reldense, K,C,F, dKe,dCe, dTprev)
	newTemp = prevTemp;	# Predict the new temp base on previous solutions
	res = 1.0; tol = 1e-4; maxItrs = 10000; itr = 0;

	# Solve linear system
	if isInLinearRange(prevTemp)
		vecPart = (F .* dt) .+ (C * prevTemp);
		matPart = (K .* dt) .+ C;
		linTemp = matPart \ vecPart;

		# Ensure solution is in linear region
		if isInLinearRange(prevTemp)
			gradTemp = computeLinearGradient(ncon[boundaries["Domain"]], matPart, C, vecPart, dKe, dCe, dt, prevTemp, dTprev);
			return linTemp, gradTemp;
		end
	end

	# Preallocate matrices
	dL, L = getSparsityPatternHeat(ncon[boundaries["Domain"]]);
	Lprev = getLatentHeatVector(ncon[boundaries["Domain"]], reldense, elements, prevTemp);

	vecPart = (F .* dt) .+ (C * prevTemp) .+ Lprev;	# Made global for linesearch (ReverseDiff limitation)
	matPart = (C .+ K .* dt)							# Made global for linesearch (ReverseDiff limitation)

	# Perform newton iterations
	gradTemp = dTprev;
	dLprev,_ = computeLatentHeatGradient(ncon[boundaries["Domain"]], elements, reldense, prevTemp, dTprev, false);
	while res > tol && itr < maxItrs 
		iterateLatentHeat!(ncon[boundaries["Domain"]], reldense, elements, dL, L, newTemp);
		residual = (vecPart .- L) .- (matPart) * newTemp;
		jacobian = matPart .+ dL;
		direction = jacobian \ residual;

		gradTemp += computeNonlinearGradient(reldense, ncon[boundaries["Domain"]], matPart, jacobian, C, residual, dKe, dCe, dLprev, dt, newTemp, prevTemp, dTprev, gradTemp);
		
		stepsize = 1.0;
		newTemp =  newTemp .+ stepsize .* direction;
		res = sqrt(residual'*residual);
		itr += 1;
		Printf.@printf "Newton Iteration %i: res = %.2e, alpha = %.2f\n" itr res stepsize
	end;
	return newTemp, gradTemp;
end;
@everywhere function isInLinearRange(T)
	return (maximum(T) < tm-tr || minimum(T) > tm+tr);
end;


## Solves a single timestep
@everywhere function solveTimestep(dt, prevTemp, deltaTemp, reldense, K,C,F)
	# useLineSearch = true;
	global newTemp = prevTemp# + 0.5.*deltaTemp;	# Predict the new temp base on previous solutions
	global Tprev = prevTemp;
	res = 1.0; tol = 1e-4; maxItrs = 10000; itr = 0;

	# Solve linear system
	if maximum(prevTemp) < tm-tr || minimum(prevTemp) > tm+tr
		residual = (F .* dt) .+ (C * prevTemp);
		jacobian = (K .* dt) .+ C;
		newTemp = jacobian \ residual;

		# Ensure solution is in linear region
		if maximum(newTemp) < tm - tr || minimum(newTemp) > tm+tr
			return newTemp
		end
	end

	# Preallocate matrices
	dL, L = getSparsityPatternHeat(ncon[boundaries["Domain"]]);
	Lprev = getLatentHeatVector(ncon[boundaries["Domain"]], reldense, elements, prevTemp);

	global vecPart = (F .* dt) .+ (C * prevTemp) .+ Lprev;	# Made global for linesearch (ReverseDiff limitation)
	global matPart = (C .+ K .* dt)							# Made global for linesearch (ReverseDiff limitation)

	# Perform newton iterations
	while res > tol && itr < maxItrs 
		iterateLatentHeat!(ncon[boundaries["Domain"]], reldense, elements, dL, L, newTemp);
		residual = (vecPart .- L) .- (matPart) * newTemp;
		jacobian = matPart .+ dL;
		direction = jacobian \ residual;

		# stepsize, m = linesearchBacktrack(newTemp, 0.9, direction);
		m = 0.0; stepsize = 1.0;
		newTemp =  newTemp .+ stepsize .* direction;
		res = sqrt(residual'*residual);
		itr += 1;
		Printf.@printf "Newton Iteration %i: res = %.2e, alpha = %.2f, m = %.2e\n" itr res stepsize m
	end;
	return newTemp;
end;


## ----- Line Search Functions ----- ##
@everywhere function residualNorm2(x)
	L = getLatentHeatVector(ncon[boundaries["Domain"]], reldense, elements, x);
	residual = (vecPart .- L) .- (matPart) * x;
	return residual' * residual;
end;

@everywhere function linesearchBacktrack(x, stepsize, direction)
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


## ----- Helper functions ----- ##
@everywhere function handleOverwrite(savename)
	if isfile(savename)
		return true;
	end
	return false;
end;