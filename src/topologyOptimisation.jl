include("latticeProperties.jl")
using ForwardDiff

@everywhere function sensitivityAnalysis(ncon, relDense, p, u, Ke, HF, nElements, nDofs)
	# Note objective function is min: c = T' * f; s.t: int p_e = gamma * f (where p_e is the element density, gamma is the volume fraction)
	# Compute element sensitivities for each element. dc/dp
	sensitivity = Vector{Float64}(undef, nElements)
	Threads.@threads for e = 1:nElements
		uElement = zeros(nDofs * 4)
		for i = 1:4
			for j = 1:nDofs
				uElement[nDofs*i - (nDofs - j)] = u[nDofs*ncon[e][i] - (nDofs - j)]
			end
		end
		sensitivity[e] = -uElement' * (Ke[:,:,e]) * uElement
	end
	rowSum = zeros(size(HF)[1]);
	for i = 1:size(HF)[1]
		rowSum[i] = sum(HF[i,:]);
	end

	return HF * (relDense .* sensitivity) ./ rowSum ./ max.(1e-3, relDense)	# Filter sensitivities
end


@everywhere function updateDensity(sensitivity, relDense, volFrac, nElements)
	# Bisection method to find optimum Lagrange Multiplier
	# This is essentially a basic root finder
	# Refer to the paper (https://doi.org/10.3390/app11073175) for equation numbers
	# Currently this uses the 'non-generalised' version as I only have 1 constraint
	relDenseNew = []
	lambda1 = 0.0; lambda2 = 1.0e9; 
	move = 0.2	# This parameter controls how much an element density can change. Too large and a local optimum is found too quickly which isnt close to the global optimum
	while (lambda2 - lambda1) / (lambda1 + lambda2) > 1e-4	# While root not found
		mid = 0.5 * (lambda1 + lambda2)	# Bisect
		B = -sensitivity / mid			# Equation 6
		# xb = relDense .* B.^0.4	# Equation 7 (part 1)
		xb = relDense .* (1.0 .+ (B .- 1.0)./(1.0 * 2.0));

		# Equation 7 (fully implimented)
		relDenseNew = max.(0.01, max.(relDense.-move, min.(1.0, min.(relDense.+move, xb))))

		# This checks the dL/dlambda = int p_e - gamma = 0 condition and updates the bisection accordingly
		if sum(relDenseNew) - (volFrac * nElements) > 0
			lambda1 = mid
		else
			lambda2 = mid
		end
	end

	return relDenseNew
end

#=
	This GOCM function updates the density using the generalised version of
	the above function.
	
	The constraint is not necessarily satisfied in all iterations.
	In its current form, only the relative volume constraint is implimented.
	This will be generalised to be able to include number of constraints assuming
		their derivatives are known. For example, a minimum displacement can be 
		included if the derivative can be found. This will require more work
=#
@everywhere function updateDensityGOCM(	# Inputs
		gl, 				# Previous constraint 
		sensitivity, 		# Derivative of objective function w.r.t relDense
		relDense, 			# Relative density of all elements
		volFrac, 			# Desired volume fraction
		nElements, 			# Number of elements
		lambda)				# Lagrange multiplier(s)


	esp = 0.05; move = 0.1;
	g = sum(relDense) / (nElements * volFrac) - 1;	# Constraint
	dg = g - gl										# Change in constraint from Previous iteration

	tol = 1e-10
	# Determine update parameter (essentially controls how fast the lagrange multiplier can update)
	if (g > tol && dg > tol) || (g < -tol && dg < -tol);		p0 = 1.0;
	elseif (g > tol && dg > -esp) || (g < -tol && dg < esp);	p0 = 0.5;
	else;														p0 = 0.0;
	end

	
	lambda = lambda * (1.0 + p0 * (g + dg));	# Update lagrange multiplier
	D = -sensitivity ./ (lambda / nElements);	# Determine scale factor
	return g, lambda, max.(0.001, max.(relDense.-move, min.(1.0, min.(relDense.+move, relDense .* abs.(D) .^0.4))));	# return new constraint, lambda, and updated relative density
end