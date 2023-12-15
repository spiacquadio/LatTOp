using LinearAlgebra
using Ferrite

# Material properties
k = 10.0;		# Conductivity	[W/m.K]
rho = 910.0;	# Density		[kg/m3]
Cp = 2000.0;	# Specific heat (at constant pressure)	[J/kg.K]
LH = 190000.0;	# Latent heat of melting				[J/kg]

tm = 273.15 + 29;	# Melting temperature
tr = 10.0;			# Half Mushy zone temperature range (T_mushy = [tm-tr : tm+tr])


## Gets the fraction of the phase based on its temperature
function getPhaseFraction(temp)
	a = 15.0;
	if (temp < tm - tr)
		return 0.0, 0.0;
	elseif (temp > tm + tr)
		return 1.0, 0.0;
	else
		s = 1 / (1 + exp(-a * (temp-tm)));
		ds = (a * exp(-a * (temp-tm))) / (1 + exp(-a * (temp-tm)))^2;
		return s, ds;
	end
end;

## Gets only the Latent Heat vector
function getLatentHeatVector(ncon, elements, temp)
	L = Vector{typeof(temp[1])}(undef,nNodes);
	fill!(L, 0.0);
	for e = 1:nElements
		Le = Vector{typeof(temp[1])}(undef,4);		# Init element matrix
		fill!(Le, 0.0);
		pf, _ = getPhaseFraction(temp[e]);			# Get phase fraction
		Le .= (rho * LH * pf) .* elements[e].fe;	# Compute element latent heat vector

		# Assemble to global latent heat vector
		elementMap = ncon[e];
		L[elementMap] += Le;
	end

	return L;
end;

## Iterates the latent heat vector given a new temperature
function iterateLatentHeat!(ncon, elements, dL, L, temp)
	# Reset matrix and vector
	fill!(L, 0.0);
	fill!(dL, 0.0);

	assembler = start_assemble(dL, L);
	for e = 1:nElements
		Le = zeros(4);
		dLe = zeros(4,4);

		# Compute matrices
		pf, dpf = getPhaseFraction(temp[e]);
		Le .= (rho * LH * pf) .* elements[e].fe;
		dLe .= (rho * LH * dpf) .* elements[e].Me;

		# Assemble to global matrix and vector
		elementMap = ncon[e]
		assemble!(assembler, elementMap, dLe, Le);
	end
end;


## Gets the constant matrices and vectors
function assembleElementPC(element, quad, shape)
	Ke = zeros(4,4);
	Ce = zeros(4,4);
	fe = zeros(4);
	
	for q = 1:quad.nQuadPoints
		detJ, dNdX = getJacobian(q, shape, element);
		integrator = quad.weights[q] * detJ

		dNdX1 = dNdX[1,:]	# for ease of access	
		dNdX2 = dNdX[2,:]	# for ease of access

		for i = 1:4
			fe[i] += 0.0;
			for j = 1:4
				Ke[i,j] += integrator * k * (dNdX1[i] * dNdX1[j] + dNdX2[i] * dNdX2[j]);
				Ce[i,j] += integrator * rho * Cp * (shape.N[q,i] * shape.N[q,j]);
			end
		end
	end

	return Ke, Ce, fe;
end;

## Assembles the constant global matrices and vector
function globalAssembly!(ncon, coords, K, C, F, quad, shape)
	fill!(K, 0.0); fill!(C, 0.0); fill!(F, 0.0);

	assemblerK = start_assemble(K);
	assemblerC = start_assemble(C);

	for e = 1:nElements
		elePoints = coords[ncon[e],:]
		Ke, Ce, fe = assembleElementPC(elePoints, quad, shape)

		elementMap = ncon[e]
		assemble!(assemblerK, elementMap, Ke);
		assemble!(assemblerC, elementMap, Ce);
		F[elementMap] += fe;
	end
end

function handleNeumannBCs!(
		constraints, 		# Node connectivity of the constrained elements
		pos, 				# Nodal coordinates
		value, 				# Neumann BC value (size = nNodalDofs x 1)
		f)					# Global RHS

	nElements = size(constraints)[1]
	fe = zeros(2)
	for e = 1:nElements
		idx = constraints[e]
		elePoints = pos[idx,1:2]
		eleLength = norm(elePoints[2,:] - elePoints[1,:])	# weight * det(J) = eleLength (2 * L/2 = L)

		N = eleLength * 0.5
		for i = 1:2
			fe[i] = value * N
		end
		
		for i = 1:2
			f[idx[i]] += fe[i]
		end
	end
end;

function getSparsityPatternHeat(ncon)
	nElements = size(ncon)[1]
	nDofs = maximum(maximum.(ncon))

	# Create temporary arrays to store which elements are non-zero
	indexI = [Vector{Int64}(undef,0) for _ = 1:nDofs]	# Non-zeros in the rows
	indexJ = [Vector{Int64}(undef,0) for _ = 1:nDofs]	# Non-zeros in the columns
	for e = 1:nElements
		element = ncon[e]
		for r in element
			for c in element
				indexI[r] = [indexI[r]; r]	# Add the non-zeros to the relevant array
				indexJ[r] = [indexJ[r]; c]	# Add the non-zeros to the relevant array
			end
		end
	end

	indexI = collect(Iterators.flatten(indexI))	# Flatten vectors into 1D arrays for constructing sparse matrix
	indexJ = collect(Iterators.flatten(indexJ))	# Flatten vectors into 1D arrays for constructing sparse matrix		
	values = zeros(Float64, length(indexI))		# Create a values array of 0's
	return sparse(indexI, indexJ, values), Vector{Float64}(undef, nDofs)	# Return the sparce matrix and full RHS
end;