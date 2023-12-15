include("shapeFunctions.jl")

using SparseArrays

function getDmatrix(E)
	Eeval = 1
	D = (Eeval / (1 - nu^2)) * [1; nu; 0;; nu; 1; 0;; 0; 0; 0.5*(1-nu)]
	return D
end

function assembleElement!(E, relDense, element, e, Ke, fe, KeSave, quad::QuadratureRules, shape::ShapeFunctions)
	# Initialise
	nEleNodes = size(Ke)[1]
	for i = 1:nEleNodes
		fe[i] = 0.0
		for j = 1:nEleNodes
			Ke[i,j] = 0.0
			KeSave[i,j,e]
		end
	end

	# Compute Elemental Matrix and Vector
	B = zeros(3,8)
	for q = 1:quad.nQuadPoints
		detJ, dNdX = getJacobian(q, shape, element)
		getBmatrix!(B,dNdX[1,:], dNdX[2,:])

		integrator = quad.weights[q] * detJ
		for i = 1:nEleNodes
			fe[i] += 0
			for j = 1:nEleNodes
				KeSave[i,j,e] += integrator * B[:,i]' * getDmatrix(E) * B[:,j]
				Ke[i,j] += relDense[e] * KeSave[i,j,e]
			end
		end
	end
end;

function globalAssembly!(E, relDense, ncon, coords, KeSave, K::SparseMatrixCSC, f::Vector{Float64}, quad::QuadratureRules, shape::ShapeFunctions)
	nNodes = size(coords)[1]
	nElements = size(ncon)[1]
	nDofs = 2 * nNodes

	# Initialise
	for i = 1:nDofs
		f[i] = 0.0
	end
	for i = 1:length(K.nzval)
		K.nzval[i] = 0.0
	end

	# Assemble
	Ke = Matrix{Float64}(undef, 8,8)
	fe = Vector{Float64}(undef, 8)
	for e = 1:nElements
		if length(ncon[e]) != 4
			continue;
		end

		elePoints = coords[ncon[e],1:2]
		assembleElement!(E, relDense, elePoints, e, Ke, fe, KeSave, quad, shape)
		for i = 1:length(ncon[e])
			f[2*ncon[e][i]-1] += fe[2*i-1]
			f[2*ncon[e][i]] += fe[2*i]
			for j = 1:length(ncon[e])
				K[2*ncon[e][i]-1,2*ncon[e][j]-1] += Ke[2*i-1,2*j-1]
				K[2*ncon[e][i]-1,2*ncon[e][j]] += Ke[2*i-1,2*j]
				K[2*ncon[e][i],2*ncon[e][j]-1] += Ke[2*i,2*j-1]
				K[2*ncon[e][i],2*ncon[e][j]] += Ke[2*i,2*j]
			end
		end
	end
end;

@enum DegreeOfFreedom begin
	xDof = 1 << 0;
	yDof = 1 << 1;
	both = (1 << 0) | (1 << 1);
end

function handleDirichletBCs!(constraints::Vector, dof::DegreeOfFreedom, value, K, f)
	nNodes = size(K)[1]
	for c = 1:length(constraints)
		idx = constraints[c]
		for i = 1:nNodes
			if Int(dof) & Int(xDof) != 0; K[2*idx-1, i] = 0.0 	end
			if Int(dof) & Int(yDof) != 0; K[2*idx, i] = 0.0		end
		end
		if Int(dof) & Int(xDof) != 0
			K[2*idx-1, 2*idx-1] = 1.0 	
			f[2*idx-1] = value
		end
		if Int(dof) & Int(yDof) != 0
			K[2*idx, 2*idx] = 1.0
			f[2*idx] = value
		end
	end
end;

using LinearAlgebra
function handleNeumannBCs!(constraints, pos, values, f)
	nElements = size(constraints)[1]
	fe = zeros(4)
	for e = 1:nElements
		idx = constraints[e]
		elePoints = pos[idx,1:2]
		eleLength = norm(elePoints[2,:] - elePoints[1,:])	# weight * det(J) = eleLength (2 * L/2 = L)

		N1 = eleLength * 0.5
		N2 = eleLength * 0.5
		fe[1] = values[1] * N1
		fe[2] = values[2] * N1
		fe[3] = values[1] * N2
		fe[4] = values[2] * N2

		for i = 1:2
			f[2*idx[i]-1] += fe[2*i-1]
			f[2*idx[i]] += fe[i]
		end
	end
end;

function getSparsityPattern(ncon)
	nElements = size(ncon)[1]
	nDofs = 2 * maximum(maximum.(ncon))

	indexI = [Vector{Int64}(undef,0) for _ = 1:nDofs]
	indexJ = [Vector{Int64}(undef,0) for _ = 1:nDofs]
	for e = 1:nElements
		element = ncon[e]
		for r in element
			for c in element
				indexI[2*r] = [indexI[2*r]; 2*r]
				indexJ[2*r] = [indexJ[2*r]; 2*c]
				indexI[2*r-1] = [indexI[2*r-1]; 2*r-1]
				indexJ[2*r-1] = [indexJ[2*r-1]; 2*c-1]
			end
		end
	end

	indexI = collect(Iterators.flatten(indexI))
	indexJ = collect(Iterators.flatten(indexJ))
	values = zeros(Float64, length(indexI))
	return sparse(indexI, indexJ, values), Vector{Float64}(undef, nDofs)
end;

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