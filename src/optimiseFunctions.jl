include("phaseChange.jl")
using Statistics
#=
	These functions are boundary condition helpers which are used in the topology Optimization
	function. This way the optimiseTopology function doesn't have to be edited
=#
@everywhere function thermalBCs!(K, f)
	handleNeumannBCs!(ncon[boundaries["Heat"]], pos, 25000.0, f, 1)
end
@everywhere function solidBCs!(K, f)
	handleDirichletBCs!(boundaryNodes["Left"], xyDof, [0.0, 0.0], K, f, 2)
	handleNeumannBCs!(ncon[boundaries["TractionBottom"]], pos, [0.0, -10000.0], f, 2)
	# f[2*4+1] = -1.0
end

#=
	Computes the compliance vs compliance pareto front
=#
@everywhere function generateParetoFront(volfrac)
	nPoints = 21;
	weights = range(0.0,1.0, nPoints);
	
	c = zeros(nPoints, 2);
	 for i = 1:nPoints
		weight = weights[i];
		println("[Pareto Generation] Iteration: $i, Weight $weight");
		reldense, _, _ = optimiseTopology(volfrac, weights[i], "", false);
		_, cT = solveHeat(reldense, "", false);
		_, cS = solveSolid(reldense, "", false);
		c[i,:] = [cT, cS];
	end

	return c;
end;

#=
	Determines the pareto optimal structure for a given 
		utopia point
=#
@everywhere function getParetoOptimalStructure(utopia, volfrac)
	reldense, cTmin = optimiseTopologyThermal(volfrac, "", false);
	_, cSmax = solveSolid(reldense, "", false);

	reldense, cSmin = optimiseTopologySolid(volfrac, "", false);
	_, cTmax = solveHeat(reldense, "", false);

	# Used to normalise
	cSrange = (cSmax - cSmin);
	cTrange = (cTmax - cTmin);
	#cTnorm=cTrange/cTmax;
	#cSnorm=cSrange/cSmax;
	# Initialise parameters
	guess = getInitialGuess([cTmax,cSmin]./[cTrange, cSrange], [cTmin,cSmax]./[cTrange, cSrange], utopia./[cTrange, cSrange]);
	#guess = getInitialGuess([cTrange, cSrange]./[cTmax,cSmin], [cTrange, cSrange]./[cTmin,cSmax], utopia./[cTnorm, cSnorm]);#changed by Stefano
	println("Initial guess is $guess");

	reldense, _, _ = optimiseTopology(volfrac, guess, "", false);
	_, cT = solveHeat(reldense, "", false);
	_, cS = solveSolid(reldense, "", false);

	a = 0.0; b = guess; c = 1.0;
	ga = goalFunction(cTmax, cSmin, utopia, cTrange, cSrange);
	gb = goalFunction(cT, cS, utopia, cTrange, cSrange);
	gc = goalFunction(cTmin, cSmax, utopia, cTrange, cSrange);
	parabola = [a,b,c, ga,gb,gc];

	# Find optimal point
	change = 1.0; progress = guess; itr = 0; xOpt = 0.0;
	while change > 0.01
		itr += 1;
		prev = [a,b,c];
		guess = 0.5 * (ga*(b^2-c^2) + gb*(c^2-a^2) + gc*(a^2-b^2)) / (ga*(b-c) + gb*(c-a) + gc*(a-b));

		# Bisection if outside limits
		if guess < 0
			a = 0;
			b = (b - a) / 2;
			ga = gb;
		elseif guess > 1
			b = (c - b) / 2;
			c = 1;
			gc = gb;

		else
			reldense, _, _ = optimiseTopology(volfrac, guess, "", false);
			_, cT = solveHeat(reldense, "", false);
			_, cS = solveSolid(reldense, "", false);
			gx = goalFunction(cT, cS, utopia, cTrange, cSrange);
			if guess > b
				if gx > gb
					c = guess;
					gc = gx;
				else
					a = b;
					b = guess;
					ga = gb;
					gb = gx;
				end
			elseif guess < b
				if gx > gb
					a = guess;
					ga = gx;
				else
					c = b;
					b = guess;
					gc = gb;
					gb = gx;
				end
			end
		end

		g = [ga,gb,gc];
		xv = [a,b,c];
		parabola = [parabola;; xv;g]

		xOpt = xv[findmin(g)[2]];
		change = maximum(abs.(prev .- xOpt));
		progress = [progress; xOpt];

		Printf.@printf "Itr: %i - Change: %.2f - Error: %.2e - Guess: %.2f - Compliance: %.6e, %.6e\n" itr change minimum(g) xOpt cT cS
	end
	return reldense, cT, cS, xOpt, progress, parabola;
end;

@everywhere function goalFunction(cT,cS, utopia, cTrange,cSrange, p=3)
	err = (abs((cT - utopia[1])/cTrange)^p + abs((cS - utopia[2])/cSrange)^p)^(1.0/p);
end;

@everywhere function getInitialGuess(p1, p2, px)
	vec = p1 - p2;
	veclength = sqrt(vec[1]^2 + vec[2]^2);
	vec = vec ./ veclength;

	vecpx = p1 - px;
	distance = vecpx[1]*vec[1] + vecpx[2]*vec[2];

	return distance / veclength;
end;

#=
	Topology Optimization function
	Optimises topology using the SIMP method

	Inputs:
		- k (temporary variable which will be used for porosity of cell)
		- volFrac (desired volume fraction)

	Outputs:
		- relDense (array of element relative densities)
		- c (domain compliance)
=#
@everywhere function optimiseTopology(volFrac, weight, savename = "", save = true)
	global p = 1.0
	itr = 0; pmax = 3; change = 1.0

	w1 = weight;
	w2 = 1.0 - weight;

	# GOCM Stuff
	lambda = 1.0;
	constraint = 0.0;
	first = true; f0 = []

	# Preallocate some arrays
	u_T = []; c_T = 0.0;
	u_S = []; c_S = 0.0;
	relDense = volFrac * ones(Float64, nElements)

	while (itr < 50)
		# Main FE solve
		FEAtime = @elapsed begin
			# Solve thermal problem
			global elementType = 0
			globalAssembly!(relDense, ncon[boundaries["Domain"]], Ke_T, K_T, f_T, quad, shape, 1)
			thermalBCs!(K_T, f_T)
			u_T = K_T \ f_T	# Solve for temperature
			c_T = u_T' * f_T

			# Solve solid-elastic problem
			global elementType = 1
			globalAssembly!(relDense, ncon[boundaries["Domain"]], Ke_S, K_S, f_S, quad, shape, 2)
			solidBCs!(K_S, f_S)
			u_S = K_S \ f_S	# Solve for temperature
			c_S = u_S' * f_S
		end


		# Compute element sensitivities (Computes dc/de)
		sensitivityTime = @elapsed begin
			# Increase penalisation parameter after 10 iterations to improve convergence
			if itr > 10 && p < pmax
				p *= 1.02
			end
			sensitivity_T = sensitivityAnalysis(ncon[boundaries["Domain"]], relDense, p, u_T, Ke_T, HF, nElements, 1)
			sensitivity_S = sensitivityAnalysis(ncon[boundaries["Domain"]], relDense, p, u_S, Ke_S, HF, nElements, 2)
			
			if first == true
				f0 = [norm(sensitivity_T), norm(sensitivity_S)]
			end
			sensitivity_T = sensitivity_T ./ f0[1];
			sensitivity_S = sensitivity_S ./ f0[2];
			sensitivity = (w1 .* sensitivity_T) .+ (w2 .* sensitivity_S)
		end

		# Update element densities 
		bisectionTime = @elapsed begin
			relDenseNew = updateDensity(sensitivity, relDense, volFrac, nElements)
			# constraint, lambda, relDenseNew = updateDensityGOCM(constraint, sensitivity, relDense, volFrac, nElements, lambda)
		end

		# Iterate
		itr += 1
		change = maximum(abs.(relDense - relDenseNew))
		relDense = relDenseNew	# Update relDense
		
		Printf.@printf "Itr: %i - Conv: %.2e - Compliance: %.2e, %.2e - Constraint: %.2e [TIMERS] FEA = %.4f, sens = %.4f, bisection = %.4f\n" itr change c_T c_S (sum(relDense)/(volFrac*nElements)) FEAtime sensitivityTime bisectionTime
	end

	if save == false
		return relDense, c_T, c_S;
	end

	if length(savename) == 0
		println("Enter save name");
		savename = readline();
	end

	while handleOverwrite(savename) == true
		println("You are about to overwrite an existing file. Would you like to continue? (y/n)");
		println("If you would like to cancel press 'c'");
		answer = readline();
		if answer == "y"
			println("Saving")
			break;
		elseif answer == "n"
			println("Enter a new save name");
			savename = readline();
		elseif answer == "c"
			println("Cancelling");
			save = false;
			break;
		else
			println("Unrecognised key");
		end
	end

	if save == true
		export2vtk(savepath * savename * ".vtu", pos, ncon, u_T, relDense)
		export2vtk(savepath * savename * ".vtu", pos, ncon, u_S, relDense)
	end
	return relDense, c_T, c_S
end;

@everywhere function optimiseTopologyThermal(volFrac, savename = "", save = true)
	global p = 1.0;
	itr = 0; pmax = 3.0; change = 1.0

	# GOCM Stuff
	lambda = 1.0;
	constraint = 0.0;
	first = true; f0 = []

	# Preallocate some arrays
	u_T = []; c_T = 0.0;
	relDense_T = volFrac * ones(Float64, nElements)

	while (itr < 50)
		# Main FE solve
		FEAtime = @elapsed begin
			# Solve thermal problem
			global elementType = 0
			globalAssembly!(relDense_T, ncon[boundaries["Domain"]], Ke_T, K_T, f_T, quad, shape, 1)
			thermalBCs!(K_T, f_T)
			u_T = K_T \ f_T	# Solve for temperature
		end


		# Compute element sensitivities (Computes dc/dde)
		sensitivityTime = @elapsed begin
			# Increase penalisation parameter after 10 iterations to improve convergence
			if itr > 10 && p < pmax
				p *= 1.02
			end
			sensitivity_T = sensitivityAnalysis(ncon[boundaries["Domain"]], relDense_T, p, u_T, Ke_T, HF, nElements, 1)
			
			# if first == true
			# 	f0 = u_T' * f_T
			# end
			# sensitivity_T = sensitivity_T ./ f0
		end

		# Update element densities 
		bisectionTime = @elapsed begin
			relDenseNew_T = updateDensity(sensitivity_T, relDense_T, volFrac, nElements)
			# constraint, lambda, relDenseNew_T = updateDensityGOCM(constraint, sensitivity_T, relDense_T, volFrac, nElements, lambda)
		end

		# Iterate
		itr += 1
		change = maximum(abs.(relDense_T - relDenseNew_T))
		c_T = u_T' * f_T			# Compute compliance
		relDense_T = relDenseNew_T	# Update relDense
		
		Printf.@printf "Itr: %i - Conv: %.2e - Compliance: %.2e - Constraint: %.2e [TIMERS] FEA = %.4f, sens = %.4f, bisection = %.4f\n" itr change c_T (sum(relDense_T)/(volFrac*nElements)) FEAtime sensitivityTime bisectionTime
	end

	if save == false
		return relDense_T, c_T;
	end

	if length(savename) == 0
		println("Enter save name");
		savename = readline();
	end

	while handleOverwrite(savename) == true
		println("You are about to overwrite an existing file. Would you like to continue? (y/n)");
		println("If you would like to cancel press 'c'");
		answer = readline();
		if answer == "y"
			println("Saving")
			break;
		elseif answer == "n"
			println("Enter a new save name");
			savename = readline();
		elseif answer == "c"
			println("Cancelling");
			save = false;
			break;
		else
			println("Unrecognised key");
		end
	end

	if save == true
		export2vtk(savepath * savename * ".vtu", pos, ncon, u_T, relDense_T)
	end
	return relDense_T, c_T
end;

@everywhere function optimiseTopologySolid(volFrac, savename = "", save = true)
	global p = 1.0; 
	itr = 0; pmax = 3.0; change = 1.0

	# GOCM stuff
	lambda = 1.0;
	constraint = 0.0;
	first = true; f0 = []

	# Preallocate some arrays
	u_S = []; c_S = 0.0;
	relDense_S = volFrac * ones(Float64, nElements)

	while (itr < 50)
		# Main FE solve
		FEAtime = @elapsed begin
			# Solve solid-elastic problem
			global elementType = 1;
			globalAssembly!(relDense_S, ncon[boundaries["Domain"]], Ke_S, K_S, f_S, quad, shape, 2)
			solidBCs!(K_S, f_S)
			u_S = K_S \ f_S	# Solve for temperature
		end


		# Compute element sensitivities (Computes dc/de)
		sensitivityTime = @elapsed begin
			# Increase penalisation parameter after 10 iterations to improve convergence
			if itr > 10 && p < pmax
				p *= 1.02
			end
			sensitivity_S = sensitivityAnalysis(ncon[boundaries["Domain"]], relDense_S, p, u_S, Ke_S, HF, nElements, 2)
			
			if first == true
				# f0 = u_S' * f_S
				f0 = norm(sensitivity_S);
			end
			sensitivity_S = sensitivity_S ./ f0
		end

		# Update element densities 
		bisectionTime = @elapsed begin
			relDenseNew_S = updateDensity(sensitivity_S, relDense_S, volFrac, nElements)
			# constraint, lambda, relDenseNew_S = updateDensityGOCM(constraint, sensitivity_S, relDense_S, volFrac, nElements, lambda)
		end

		# Iterate
		itr += 1
		change = maximum(abs.(relDense_S - relDenseNew_S))
		c_S = u_S' * f_S	# Compute compliance
		relDense_S = relDenseNew_S	# Update relDense
		
		Printf.@printf "Itr: %i - Conv: %.2e - Compliance: %.2e - Constraint: %.2e [TIMERS] FEA = %.4f, sens = %.4f, bisection = %.4f\n" itr change c_S (sum(relDense_S)/(volFrac*nElements)) FEAtime sensitivityTime bisectionTime
	end

	if save == false
		return relDense_S, c_S;
	end

	if length(savename) == 0
		println("Enter save name");
		savename = readline();
	end

	while handleOverwrite(savename) == true
		println("You are about to overwrite an existing file. Would you like to continue? (y/n)");
		println("If you would like to cancel press 'c'");
		answer = readline();
		if answer == "y"
			println("Saving")
			break;
		elseif answer == "n"
			println("Enter a new save name");
			savename = readline();
		elseif answer == "c"
			println("Cancelling");
			save = false;
			break;
		else
			println("Unrecognised key");
		end
	end

	if save == true
		export2vtk(savepath * savename * ".vtu", pos, ncon, u_S, relDense_S)
	end

	return relDense_S, c_S
end;

@everywhere function optimiseTopologyPhaseChange(volFrac, endTime, savename = "", save = true)
	global p = 1.0; 
	itr = 0; pmax = 3.0; change = 1.0

	# GOCM stuff
	lambda = 1.0;
	constraint = 0.0;
	first = true; f0 = []

	# Preallocate some arrays
	u_PC = []; c_PC = 0.0;
	relDense_PC = volFrac * ones(Float64, nElements)

	while (itr < 50)
		# Main FE solve
		FEAtime = @elapsed begin
			# Solve phase-change problem
			dt = 1.0e-3;
			initialTemp = 18.0+273.15;
			u_PC, c_PC, sensitivity_PC = solveWithGradient(endTime, dt, initialTemp, relDense_PC);
		end


		# Filter sensitivities
		sensitivityTime = @elapsed begin
			# Increase penalisation parameter after 10 iterations to improve convergence
			if itr > 10 && p < pmax
				p *= 1.02
			end

			# if first == true
			# 	# f0 = u_S' * f_S
			# 	f0 = norm(sensitivity_PC);
			# end
			# sensitivity_PC = sensitivity_PC ./ f0

			rowSum = zeros(size(HF)[1]);
			for i = 1:size(HF)[1]
				rowSum[i] = sum(HF[i,:]);
			end
			sensitivity_PC = HF * (relDense_PC .* sensitivity_PC) ./ rowSum ./ max.(1e-3, relDense_PC)
		end

		# Update element densities 
		bisectionTime = @elapsed begin
			relDenseNew_PC = updateDensity(sensitivity_PC, relDense_PC, volFrac, nElements)
		end

		# Iterate
		itr += 1
		change = maximum(abs.(relDense_PC - relDenseNew_PC))
		relDense_PC = relDenseNew_PC	# Update relDense
		
		Printf.@printf "Itr: %i - Conv: %.2e - Compliance: %.2e - Constraint: %.2e [TIMERS] FEA = %.4f, sens = %.4f, bisection = %.4f\n" itr change c_PC (sum(relDense_PC)/(volFrac*nElements)) FEAtime sensitivityTime bisectionTime
	end

	if save == false
		return relDense_PC, c_PC;
	end

	if length(savename) == 0
		println("Enter save name");
		savename = readline();
	end

	while handleOverwrite(savename) == true
		println("You are about to overwrite an existing file. Would you like to continue? (y/n)");
		println("If you would like to cancel press 'c'");
		answer = readline();
		if answer == "y"
			println("Saving")
			break;
		elseif answer == "n"
			println("Enter a new save name");
			savename = readline();
		elseif answer == "c"
			println("Cancelling");
			save = false;
			break;
		else
			println("Unrecognised key");
		end
	end

	if save == true
		export2vtk(savepath * savename * ".vtu", pos, ncon, u_PC, relDense_PC)
	end

	return relDense_PC, c_PC
end;

@everywhere function solveSolid(relDense_S, savename = "", save = true)
	global p = 1.0;

	global elementType = 1;
	globalAssembly!(relDense_S, ncon[boundaries["Domain"]], Ke_S, K_S, f_S, quad, shape, 2)
	solidBCs!(K_S, f_S)
	u_S = K_S \ f_S

	if save == false
		return u_S, u_S'*f_S;
	end

	if length(savename) == 0
		println("Enter save name or press 'c' to skip");
		savename = readline();
		if savename == "c"
			save = false;
		end
	end

	while handleOverwrite(savename) == true
		println("You are about to overwrite an existing file. Would you like to continue? (y/n)");
		println("If you would like to cancel press 'c'");
		answer = readline();
		if answer == "y"
			println("Saving")
			break;
		elseif answer == "n"
			println("Enter a new save name");
			savename = readline();
		elseif answer == "c"
			println("Cancelling");
			save = false;
			break;
		else
			println("Unrecognised key");
		end
	end

	if save == true
		export2vtk(savepath * savename * ".vtu", pos, ncon, u_S, relDense_S)
	end

	return u_S, u_S'*f_S;
end;

@everywhere function solveHeat(relDense_T, savename = "", save = true)
	global p = 1.0;

	global elementType = 0;
	globalAssembly!(relDense_T, ncon[boundaries["Domain"]], Ke_T, K_T, f_T, quad, shape, 1)
	thermalBCs!(K_T, f_T)
	u_T = K_T \ f_T

	if save == false
		return u_T, u_T'*f_T;
	end

	if length(savename) == 0
		println("Enter save name");
		savename = readline();
	end

	while handleOverwrite(savename) == true
		println("You are about to overwrite an existing file. Would you like to continue? (y/n)");
		println("If you would like to cancel press 'c'");
		answer = readline();
		if answer == "y"
			println("Saving")
			break;
		elseif answer == "n"
			println("Enter a new save name");
			savename = readline();
		elseif answer == "c"
			println("Cancelling");
			save = false;
			break;
		else
			println("Unrecognised key");
		end
	end

	if save == true
		export2vtk(savepath * savename * ".vtu", pos, ncon, u_T, relDense_T)
	end

	return u_T, u_T'*f_T;
end;


#===== Helper Functions =====#
@everywhere function handleOverwrite(savename)
	if isfile(savepath * savename * ".vtu")
		return true;
	end
	return false;
end;

#optimiseTopologyThermal(0.3, "", false)