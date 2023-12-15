include("shapeFunctions.jl")

conductivity(x) = 1/((x[1]-1)^2 + 0.25 + (x[2]-2)^2)

function assembleElement!(k, element, Ke, fe, quad::QuadratureRules, shape::ShapeFunctions)
	# Initialise
	nEleNodes = size(Ke)[1]
	for i = 1:nEleNodes
		fe[i] = 0.0
		for j = 1:nEleNodes
			Ke[i,j] = 0.0
		end
	end
	# Compute Elemental Matrix and Vector
	for q = 1:quad.nQuadPoints
		detJ, dNdX = getJacobian(q, shape, element)
		dNdX1 = dNdX[1,:]
		dNdX2 = dNdX[2,:]

		integrator = quad.weights[q] * detJ
		for i = 1:nEleNodes
			fe[i] += integrator * shape.N[q,i]
			for j = 1:nEleNodes
				Ke[i,j] += conductivity(k) * integrator * (dNdX1[i] * dNdX1[j] + dNdX2[i] * dNdX2[j])
			end
		end
	end
end

function globalAssembly(k, ncon, coords, quad::QuadratureRules, shape::ShapeFunctions)
	nNodes = size(coords)[1]
	nElements = size(ncon)[1]
	nDofs = nNodes

	K = Matrix{Float64}(undef, nDofs, nDofs)
	f = Vector{Float64}(undef, nDofs)

	# Initialise
	for i = 1:nDofs
		f[i] = 0.0
		for j = 1:nDofs
			K[i,j] = 0.0
		end
	end

	# Assemble
	Ke = Matrix{Float64}(undef, 4,4)
	fe = Vector{Float64}(undef, 4)
	for e = 1:nElements
		if length(ncon[e]) != 4
			continue;
		end

		elePoints = coords[ncon[e],:]
		assembleElement!(k, elePoints, Ke, fe, quad, shape)
		for i = 1:4
			f[ncon[e][i]] += fe[i]
			for j = 1:4
				K[ncon[e][i],ncon[e][j]] += Ke[i,j]
			end
		end
	end

	return K,f
end

function handleDirichletBCs!(constraints::Vector, value::Float64, K, f)
	nNodes = size(K)[1]
	for c = 1:length(constraints)
		idx = constraints[c]
		for i = 1:nNodes
			K[idx,i] = 0.0
		end
		K[idx,idx] = 1.0
		f[idx] = value
	end
end

function getSparsityPattern(ncon)
	nElements = size(ncon)[1]
	nDofs = maximum(maximum(ncon))

	indexI = [Vector{Int64}(undef,0) for _ = 1:nDofs]
	indexJ = [Vector{Int64}(undef,0) for _ = 1:nDofs]
	for e = 1:nElements
		element = ncon[e]
		for r in element
			for c in element
				indexI[r] = [indexI[r]; r]
				indexJ[r] = [indexJ[r]; c]
			end
		end
	end

	indexI = collect(Iterators.flatten(indexI))
	indexJ = collect(Iterators.flatten(indexJ))
	values = zeros(Float64, length(indexI))
	return sparse(indexI, indexJ, values), Vector{Float64}(undef, nDofs)
end;