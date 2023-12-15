include("shapeFunctions.jl")
include("latticeProperties.jl")

using Ferrite
using SparseArrays
using LinearAlgebra
using Printf
using Distributed
#=
	This is a "wrapper" function to easily select which element assembly 
	function to use. The "elementType" variable is in the global namespace
=#
 function assembleElement(	# Inputs
	relDense, 				# Element relative density
	element, 				# Element positions
	quad::QuadratureRules, 	# Gauss quarature rules
	shape::ShapeFunctions2D)	# Shape functions
	
	# If heat equation
	if elementType == 0
		return assembleElementHeat!(relDense, element, quad, shape);
	
	# If structural problem 
	elseif elementType == 1
		return assembleElementSolid(relDense, element, quad, shape);
	end
end

#=
	This function returns the 2D material matrix for a structural problem.
	
	It will be extended in the future to use the material properties of the 
	selected unit cell and the input will be the porosity of the cell.
=#
 function getDmatrix(reldense)
	if useLattice == true
		r = determineRadius(reldense2volfrac(reldense^p));
		D = getElasticityMatrix(r);
	else
		nu = 0.3;
		Eeval = 70.0e9 * reldense^p;
		Eeval3 = 70.0e9 * reldense^p;
		D = (1 / (1 - nu^2)) * [Eeval; Eeval*nu; 0;; Eeval3*nu; Eeval3; 0;; 0; 0; Eeval*0.5*(1-nu)]
	end
	return D
end


#=
	This is a typical finite element assembly function. It is specifically for 
	4 noded finite elements. It assembles the element matrix for a structual
	problem.

	This function returns the finite element matrix and forcing array (Ke, fe)
	as well as KeSave which is used in the sensitivity analysis in the optimise 
	topology function
=#
 function assembleElementSolid(# Inputs
	relDense, 				# Element relative density
	element, 				# Element index
	quad::QuadratureRules, 	# Gauss quarature rules
	shape::ShapeFunctions2D)	# Element shape functions


	# Initialise
	fe = zeros(8);
	Ke = zeros(typeof(relDense),8,8);
	dKe = zeros(typeof(relDense), 8,8);
	D = DiffResults.DiffResult(Matrix{Float64}(undef, 3,3), Matrix{Float64}(undef, 3,3));

	# Compute Elemental Matrix and Vector
	B = zeros(3,8);
	ForwardDiff.derivative!(D, getDmatrix, relDense);

	# Loop over quadrature points
	for q = 1:quad.nQuadPoints
		detJ, dNdX = getJacobian(q, shape, element)
		getBmatrix!(B,dNdX[1,:], dNdX[2,:])

		integrator = quad.weights[q] * detJ

		Ke += integrator .* (B' * D.value * B) 
		dKe += integrator .* (B' * DiffResults.derivative(D) * B)
	end

	return Ke, dKe, fe;
end;


#=
	This Function will be modified in the future to get the ETC of the element
	taking into account which unit cell it belongs to
=#
 function getKmatrix(reldense)
	if useLattice == true
		r = determineRadius(reldense2volfrac(reldense^p));
		D = getConductivityMatrix(r);
	else
		kMax = 160.0;
		kMin = 25e-3;
		kEval = (kMax - kMin) * reldense^p + kMin
		D = [kEval; 0;; 0; kEval];
	end
	return D;
end;


#=
	This is a typical finite element assembly function. It is specifically for 
	4 noded finite elements. It assembles the element matrix using a conductivity
	and convection matrix. 

	The reason the convection matrix (Kec) is include is because for the thermal
	topology optimisation to work correctly, the elements need to loss heat to 
	the environment. The paper was based on natural convection 

	This function returns the finite element matrix and forcing array (Ke, fe)
	as well as KeSave which is used in the sensitivity analysis in the optimise 
	topology function
=#
 function assembleElementHeat!(	# Inputs
		relDense, 				# relative density array (length = nElements)
		element, 				# Element node positions (size 4x2)
		quad::QuadratureRules, 	# Quadrature rules
		shape::ShapeFunctions2D)	# Shape functions

		
	# Initialise
	fe = zeros(4);
	Ke = zeros(typeof(relDense),4,4);
	dKe = zeros(typeof(relDense), 4,4);
	D = DiffResults.DiffResult(Matrix{Float64}(undef, 2,2), Matrix{Float64}(undef, 2,2));
		
	# Get material properties	
	h = 2.0;
	ForwardDiff.derivative!(D, getKmatrix, relDense);

	# Compute Elemental Matrix and Vector
	for q = 1:quad.nQuadPoints
		detJ, dNdX = getJacobian(q, shape, element)	# Get detminant of jacobian and shape function gradients

		integrator = quad.weights[q] * detJ
		
		fe += (integrator * 25.0 * h) .* shape.N[q,:];
		Ke += integrator .* ((dNdX' * D.value * dNdX) + (shape.N[q,:] * shape.N[q,:]'));
		dKe += integrator .* (dNdX' * DiffResults.derivative(D) * dNdX);
	end

	return Ke, dKe, fe;
end;


#=
	
=#
 function assembleElementPC(element, reldense, quad, shape)
	Ke = zeros(4,4);
	Ce = zeros(4,4);
	fe = zeros(4);
	
	D = getKmatrix(reldense)
	for q = 1:quad.nQuadPoints
		detJ, dNdX = getJacobian(q, shape, element);
		integrator = quad.weights[q] * detJ

		rho, cp, _ = getRhoCpL(reldense);

		fe .+= zeros(4);
		Ke .+= integrator .* (dNdX' * D * dNdX);
		Ce .+= (integrator * rho * cp) .* (shape.N[q,:] * shape.N[q,:]');
	end

	return Ke, Ce, fe;
end;

#=

=#
 function globalAssemblyPC!(ncon, coords, reldense, K, C, F, quad, shape)
	fill!(K, 0.0); fill!(C, 0.0); fill!(F, 0.0);
	
	assemblerK = start_assemble(K);
	assemblerC = start_assemble(C);

	 for e = 1:nElements
		elePoints = coords[ncon[e],:]
		Ke, Ce, fe = assembleElementPC(elePoints, reldense[e], quad, shape)

		elementMap = ncon[e]
		assemble!(assemblerK, elementMap, Ke);
		assemble!(assemblerC, elementMap, Ce);
		F[elementMap] += fe;
	end
end

#=

=#
## Gets the fraction of the phase based on its temperature
 function getPhaseFraction(temp)
	a = 0.5;
	if (temp < tm - tr)
		return 0.0, 0.0;
	elseif (temp > tm + tr)
		return 1.0, 0.0;
	else
		s = 1 / (1 + exp(-a * (temp-tm))); #liquid fraction
		ds = (a * exp(-a * (temp-tm))) / (1 + exp(-a * (temp-tm)))^2;
		return s, ds;
	end
	f=sum(s)./nElements;
end;
 function getdds(temp)
	a = 0.5;
	dds = (2*a^2*exp(-2*a*(temp - tm)))/(exp(-a*(temp - tm)) + 1)^3 - (a^2*exp(-a*(temp - tm)))/(exp(-a*(temp - tm)) + 1)^2;
	return dds;
end;

## Gets only the Latent Heat vector
 function getLatentHeatVector(ncon, reldense, elements, temp)
	L = Vector{typeof(temp[1])}(undef,nNodes);
	fill!(L, 0.0);
	 for e = 1:nElements
		rho, _, LH =  getRhoCpL(reldense[e])
		Le = zeros(typeof(reldense[1]),4);								# Init element matrix
		 for i = 1:4
			pf, _ = getPhaseFraction(temp[ncon[e][i]]);			# Get phase fraction
			Le[i] += (rho * LH * pf) .* elements[e].fe[i];		# Compute element latent heat vector
		end

		# Assemble to global latent heat vector
		elementMap = ncon[e];
		L[elementMap] += Le;
	end

	return L;
end;

## Iterates the latent heat vector given a new temperature
 function iterateLatentHeat!(ncon, reldense, elements, dL, L, temp)
	# Reset matrix and vector
	fill!(L, 0.0);
	fill!(dL, 0.0);

	assembler = start_assemble(dL, L);
	 for e = 1:nElements
		rho, _, LH = getRhoCpL(reldense[e]);
		Le = zeros(typeof(reldense[e]), 4);
		dLe = zeros(typeof(reldense[e]), 4,4);
		 for i = 1:4
			pf, dpf = getPhaseFraction(temp[ncon[e][i]]);
			Le[i] = (rho * LH * pf) .* elements[e].fe[i];
			dLe[i,:] = (rho * LH * dpf) .* elements[e].Me[i,:];
		end

		# Assemble to global matrix and vector
		elementMap = ncon[e]
		assemble!(assembler, elementMap, dLe, Le);
	end
end;


#= 
	This is a typical global assemble for finite element matrices.
=# 
 function globalAssembly!(		# Inputs
		relDense, 				# relative density array (length = nElements)
		ncon, 					# Node connectivity array (size = nElements x nElementNodes)
		KeSave, 				# Will be returned for sensitivity analysis
		K::SparseMatrixCSC, 	# Global stiffness matrix
		f::Vector{Float64}, 	# Global RHS
		quad::QuadratureRules, 	# Gauss quadrature rules
		shape::ShapeFunctions2D,	# Element shape functions
		nNodalDofs)				# Number of nodal degrees of freedom (1 for heat 2 for solid)


	nElements = size(ncon)[1]
	
	# Use Ferrite assembler for efficiency
	assembler = start_assemble(K, f);
	
	# Initialise system matrix/vector to 0
	fill!(f, 0.0);
	fill!(K, 0.0);

	# Assemble Linear System
	 for e = 1:nElements
		element = pos[ncon[e],1:2];
		Ke, dKe, fe = assembleElement(relDense[e], element, quad, shape);
		KeSave[:,:,e] .= dKe;
	
		# Maps the element connectivity to the global matrix indices
		elementMap = zeros(Int64, nNodalDofs * shape.nNodes)
		 for i = 1:shape.nNodes
			 for j = 1:nNodalDofs
				elementMap[nNodalDofs * i - (nNodalDofs - j)] = nNodalDofs * ncon[e][i] - (nNodalDofs - j)
			end
		end

		# Assemble element contributions to global matrix
		assemble!(assembler, elementMap, Ke, fe);
	end
end;

#=
	This function handles the Dirichlet BCs by setting the diagonal of the 
	stiffness matrix to 1 and the row to 0. It then sets f[dirichletNodes] 
	to the value supplied.

	The DOF enum is for selecting which DOF to constrain
=#
@enum DOF begin
	xDof = 1 << 0;
	yDof = 1 << 1;
	zDof = 1 << 2;
	all = (1 << 0) | (1 << 1) | (1 << 2);

	xyDof = (1 << 0) | (1 << 1);
	yzDof = (1 << 1) | (1 << 2);
	xzDof = (1 << 2) | (1 << 0);

end
function handleDirichletBCs!(	# Inputs
		constraints::Vector, 	# The nodes to be constrained
		dof::DOF,				# Constrained DOF
		value, 					# Dirichlet values (size = nNodalDofs x 1)
		K, 						# Global stiffness matrix
		f,						# Global RHS
		nNodalDofs)				# Number of nodal DOFs

	nNodes = size(K)[1]

	# Loop over each constraints
	for c in constraints
		d1 = nNodalDofs * c - (nNodalDofs - 1)	# If 1d, this evaluates to c, else 2*c-1 else 3*c-2
		d2 = nNodalDofs * c - (nNodalDofs - 2)	# If 1d, this is irrelevant, else 2*c else 3*c-1
		d3 = nNodalDofs * c - (nNodalDofs - 3)	# If 1d, this is irrelevant, else this is irrelevant else 3*c
		for i = 1:nNodes
			if (Int(dof) & Int(xDof)) != 0; K[d1, i] = 0.0	end
			if (Int(dof) & Int(yDof)) != 0; K[d2, i] = 0.0; end
			if (Int(dof) & Int(zDof)) != 0; K[d3, i] = 0.0;	end
		end
		if (Int(dof) & Int(xDof)) != 0
			K[d1, d1] = 1.0
			f[d1] = value[1]
		end
		if (Int(dof) & Int(yDof)) != 0
			K[d2, d2] = 1.0
			f[d2] = value[2]
		end
		if (Int(dof) & Int(zDof)) != 0
			K[d3, d3] = 1.0
			f[d3] = value[3]
		end
	end
end;

#=
	This function handles the Neumann BCs by using a linear finite elements
	Im not 100% sure its implimented correctly, but it seems to give the correct 
	results when I compared to other libraries
=#
function handleNeumannBCs!(
		constraints, 		# Node connectivity of the constrained elements
		pos, 				# Nodal coordinates
		value, 				# Neumann BC value (size = nNodalDofs x 1)
		f,					# Global RHS
		nNodalDofs)			# Number of nodal DOFs


	nElements = size(constraints)[1]
	fe = zeros(2 * nNodalDofs)
	for e = 1:nElements
		idx = constraints[e]
		elePoints = pos[idx,1:2]
		eleLength = norm(elePoints[2,:] - elePoints[1,:])	# weight * det(J) = eleLength (2 * L/2 = L)

		N = eleLength * 0.5
		elementMap = zeros(Int32, 2 * nNodalDofs)
		for i = 1:2
			for j = 1:nNodalDofs
				elementMap[nNodalDofs * i - (nNodalDofs - j)] = nNodalDofs * idx[i] - (nNodalDofs - j)
				fe[nNodalDofs * i - (nNodalDofs - j)] = value[j] * N
			end
		end
		
		f[elementMap] += fe
	end
end;

#=
	This function gets the sparcity pattern of the global matrix based on
	the element connectivity and then returns a sparse matrix.
=#
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
function getSparsityPatternSolid(ncon)
	nElements = size(ncon)[1]
	nDofs = 2 * maximum(maximum.(ncon))

	indexI = [Vector{Int64}(undef,0) for _ = 1:nDofs]
	indexJ = [Vector{Int64}(undef,0) for _ = 1:nDofs]
	for e = 1:nElements
		element = ncon[e]
		for r in element
			for c in element
				indexI[2*r] = [indexI[2*r]; [2*r;2*r]]
				indexJ[2*r] = [indexJ[2*r]; [2*c;2*c-1]]
				indexI[2*r-1] = [indexI[2*r-1]; [2*r-1;2*r-1]]
				indexJ[2*r-1] = [indexJ[2*r-1]; [2*c;2*c-1]]
			end
		end
	end
	indexI = collect(Iterators.flatten(indexI))
	indexJ = collect(Iterators.flatten(indexJ))
	values = zeros(Float64, length(indexI))
	return sparse(indexI, indexJ, values), Vector{Float64}(undef, nDofs)
end;


# This isnt really relevant
function getNeighbours(ncon)
	nNodes = maximum(maximum.(ncon))
	nElements = size(ncon)[1]
	connectedNodes = [Vector{Int64}(undef,0) for _ = 1:nNodes]	# Stores which elements are connected to node i
	
	for e = 1:nElements
		nodes = ncon[e]
		for node in nodes
			connectedNodes[node] = [connectedNodes[node]; e]
		end
	end

	neighbours = [Vector{Int64}(undef,0) for _ = 1:nElements]	# Stores which elements neighbour each element

	for n = 1:nNodes
		for element in connectedNodes[n]
			indices = findall(!=(element), connectedNodes[n])
			for index in indices
				neighbours[element] = [neighbours[element]; connectedNodes[n][index]]
			end
		end
	end

	for e = 1:nElements
		neighbours[e] = unique(neighbours[e])
	end
	return neighbours
end


#=
	This computes the weight factor matrix for filtering the element 
	sensitivities.

	This is constructed by finding the distance between each element, and
	storing the distance in a sparse matrix. The filter radius defines 
	how far sensititivy is filtered. Only the positive r - distance values
	are stored

	This distances are the distance between element centroids
=#
function getWeightFactorMatrix(ncon, pos, radius::Float64)
	nElements = size(ncon)[1]
	centroids = zeros(Float64, (nElements, 2))

	# Compute element centroids
	for e = 1:nElements
		element = pos[ncon[e],1:2]
		centroid = sum(element, dims=1) ./ 4
		centroids[e,:] = centroid
	end

	# Compute weight factor matrix
	indexI = [Vector{Int64}(undef,0) for _ = 1:nElements]
	indexJ = [Vector{Int64}(undef,0) for _ = 1:nElements]
	values = [Vector{Float64}(undef,0) for _ = 1:nElements]
	Threads.@threads for i = 1:nElements
		for j = 1:nElements
			distance = norm(centroids[i,:] - centroids[j,:])
			if distance < radius
				indexI[i] = [indexI[i]; i]
				indexJ[i] = [indexJ[i]; j]
				values[i] = [values[i]; radius - distance]
			end
		end
	end

	indexI = collect(Iterators.flatten(indexI))
	indexJ = collect(Iterators.flatten(indexJ))
	values = collect(Iterators.flatten(values))
	return sparse(indexI, indexJ, values)
end;