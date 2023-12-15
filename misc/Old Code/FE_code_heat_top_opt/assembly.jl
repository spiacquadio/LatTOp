include("shapeFunctions.jl")

using SparseArrays
using LinearAlgebra

#=
	This Function will be modified in the future to get the ETC of the element
	taking into account which unit cell it belongs to
=#
conductivity(x) = 25e-3 + (200 - 25e-3) * x

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
function assembleElement!(		# Inputs
		k, 						# temporary
		relDense, 				# relative density array (length = nElements)
		p, 						# Penalisation factor
		element, 				# Element node positions (size 4x2)
		e, 						# Element index
		Ke, 					# Element stiffness matrix
		fe, 					# Element forcing array
		KeSave, 				# Returns to topology optimisation for sensitivty analysis
		quad::QuadratureRules, 	# Quadrature rules
		shape::ShapeFunctions)	# Shape functions


	# Initialise matrix and array to 0
	nEleNodes = size(Ke)[1]
	Keh = zeros(4,4)
	Kec = zeros(4,4)
	for i = 1:nEleNodes
		fe[i] = 0.0
		for j = 1:nEleNodes
			KeSave[i,j,e] = 0.0
			Ke[i,j] = 0.0
		end
	end

	# Determine suitable convection coefficient
	# This is temporary code I was experiementing with. It doesnt do much because the
	#	elements all become either 0 or 1.
	h = 0.0
	if relDense[e] < 0.6
		h = 0.25
	elseif relDense[e] < 0.7
		h = 5.0
	else
		h = 0.1
	end

	# Compute Elemental Matrix and Vector
	for q = 1:quad.nQuadPoints
		detJ, dNdX = getJacobian(q, shape, element)	# Get detminant of jacobian and shape function gradients
		dNdX1 = dNdX[1,:]	# for ease of access	
		dNdX2 = dNdX[2,:]	# for ease of access

		integrator = quad.weights[q] * detJ
		for i = 1:nEleNodes
			fe[i] += 25.0 * h * integrator * shape.N[q,i]	# Convective boundary condition on all elements (ambient temp 25°C, h is the convection coefficient)
			for j = 1:nEleNodes
				Kec[i,j] += integrator * (dNdX1[i] * dNdX1[j] + dNdX2[i] * dNdX2[j])	# Assemble conductivity Matrix
				Keh[i,j] += h * integrator * (shape.N[i] * shape.N[j])					# Assemble convection matrix
				KeSave[i,j,e] += Kec[i,j]												# Only conduction is required for exnsitivity analysis
				Ke[i,j] += conductivity(relDense[e]^p) * Kec[i,j] + Keh[i,j]			# Combine matrices
			end
		end
	end
end;

#= 
	This is a typical global assemble for finite element matrices. 
	This could be generalised.
=# 
function globalAssembly!(		# Inputs
		k,						# Will be used to pass in unit cell variables (porosity of each cell most likely)
		relDense, 				# relative density array (length = nElements)
		p, 						# Penalisation factor
		ncon, 					# Node connectivity array (size = nElements x nElementNodes)
		coords, 				# Coords of each node (size = nNodes x 2)
		KeSave, 				# Will be returned for sensitivity analysis
		K::SparseMatrixCSC, 	# Global stiffness matrix
		f::Vector{Float64}, 	# Global RHS
		quad::QuadratureRules, 	# Gauss quadrature rules
		shape::ShapeFunctions)	# Element shape functions


	nNodes = size(coords)[1]
	nElements = size(ncon)[1]
	nDofs = nNodes

	# Initialise global forcing array (RHS) to 0
	for i = 1:nDofs
		f[i] = 0.0
	end

	# Initialise non-zeros in sparse matrix to 0
	for i = 1:length(K.nzval)
		K.nzval[i] = 0.0
	end

	# Assemble
	fe = Vector{Float64}(undef, 4)		# Allocate memory for fe vector
	Ke = Matrix{Float64}(undef, 4, 4)	# Allocate memory for Ke matrix
	for e = 1:nElements
		elePoints = coords[ncon[e],:]
		assembleElement!(k, relDense, p, elePoints, e, Ke, fe, KeSave, quad, shape)
		for i = 1:4
			f[ncon[e][i]] += fe[i]
			for j = 1:4
				K[ncon[e][i],ncon[e][j]] += Ke[i,j]
			end
		end
	end
end;

#=
	This function handles the Dirichlet BCs by setting the diagonal of the 
	stiffness matrix to 1 and the row to 0. It then sets f[dirichletNodes] 
	to the value supplied

	This function could be generalised
=#
function handleDirichletBCs!(	# Inputs
		constraints::Vector, 	# The nodes to be constrained
		value::Float64, 		# Dirichlet value
		K, 						# Global stiffness matrix
		f)						# Global RHS

	nNodes = size(K)[1]

	# Loop over each constraints
	for c in constraints
		for i = 1:nNodes
			K[c,i] = 0.0
		end
		K[c,c] = 1.0
		f[c] = value
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
		value, 				# Neumann BC value
		f)					# Global RHS


	nElements = size(constraints)[1]
	fe = zeros(2)
	for e = 1:nElements
		idx = constraints[e]
		elePoints = pos[idx,1:2]
		eleLength = norm(elePoints[2,:] - elePoints[1,:])	# weight * det(J) = eleLength (2 * L/2 = L)

		N1 = eleLength * 0.5
		N2 = eleLength * 0.5
		fe[1] = value * N1
		fe[2] = value * N2

		for i = 1:2
			f[idx[i]] += fe[i]
		end
	end
end;

#=
	This function gets the sparcity pattern of the global matrix based on
	the element connectivity and then returns a sparse matrix.
=#
function getSparsityPattern(ncon)
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
				#@show (element, connectedNodes[n][index])
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