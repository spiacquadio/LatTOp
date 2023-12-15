include("shapeFunctions.jl")

conductivity(x) = 1/((x[1]-1)^2 + 0.25 + (x[2]-2)^2)
function getDmatrix(E)
	Eeval = 1e3 / ((E[1] - 5.0)^2 + 1.5 + (E[2] - 7.0)^2)
	D = (Eeval / (1 - nu^2)) * [1; nu; 0;; nu; 1; 0;; 0; 0; 0.5*(1-nu)]
	return D
end

function assembleElementSolid!(x, element, Ke, fe, quad::QuadratureRules, shape::ShapeFunctions)
	# Initialise
	nEleNodes = size(Ke)[1]
	for i = 1:nEleNodes
		fe[i] = 0.0
		for j = 1:nEleNodes
			Ke[i,j] = 0.0
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
				Ke[i,j] += integrator * B[:,i]' * getDmatrix(x) * B[:,j]
			end
		end
	end
end;

function assembleElementHeat!(x, element, Ke, fe, quad::QuadratureRules, shape::ShapeFunctions)
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
				Ke[i,j] += conductivity(x) * integrator * (dNdX1[i] * dNdX1[j] + dNdX2[i] * dNdX2[j])
			end
		end
	end
end;

function globalAssembly(x, ncon, coords, quad::QuadratureRules, shape::ShapeFunctions, solveFor::String = "both")
	nNodes = size(coords)[1]
	nElements = size(ncon)[1]
	nDofs = 2 * nNodes

	KS = []
	fS = []
	KH = []
	fH = []
	Kes = []
	fes = []
	Keh = []
	feh = []
	if solveFor == "solid" || solveFor == "both"
		KS = Matrix{Float64}(undef, nDofs, nDofs)
		fS = Vector{Float64}(undef, nDofs)
		Kes = Matrix{Float64}(undef, 8,8)
		fes = Vector{Float64}(undef, 8)
	end
	if solveFor == "heat" || solveFor == "both"
		KH = Matrix{Float64}(undef, nNodes, nNodes)
		fH = Vector{Float64}(undef, nNodes)
		Keh = Matrix{Float64}(undef, 4,4)
		feh = Vector{Float64}(undef, 4)
	end

	# Initialise
	for i = 1:nNodes
		if solveFor == "solid" || solveFor == "both"
			fS[2*i] = 0.0
			fS[2*i-1] = 0.0
		end
		if solveFor == "heat" || solveFor == "both"
			fH[i] = 0.0
		end

		for j = 1:nNodes
			if solveFor == "solid" || solveFor == "both"
				KS[2*i,2*i] = 0.0
				KS[2*i-1,2*i-1] = 0.0
				KS[2*i,2*i-1] = 0.0
				KS[2*i-1,2*i] = 0.0
			end
			if solveFor == "heat" || solveFor == "both"
				KH[i,j] = 0.0
			end
		end
	end

	# Assemble
	for e = 1:nElements
		if length(ncon[e]) != 4
			continue;
		end

		elePoints = coords[ncon[e],1:2]
		if solveFor == "solid" || solveFor == "both"
			assembleElementSolid!(x, elePoints, Kes, fes, quad, shape)
		end
		if solveFor == "heat" || solveFor == "both"
			assembleElementHeat!(x, elePoints, Keh, feh, quad, shape)
		end

		for i = 1:length(ncon[e])
			if solveFor == "solid" || solveFor == "both"
				fS[2*ncon[e][i]-1] += fes[2*i-1]
				fS[2*ncon[e][i]] += fes[2*i]
			end
			if solveFor == "heat" || solveFor == "both"
				fH[ncon[e][i]] += feh[i]
			end

			for j = 1:length(ncon[e])
				if solveFor == "solid" || solveFor == "both"
					KS[2*ncon[e][i]-1,2*ncon[e][j]-1] += Kes[2*i-1,2*j-1]
					KS[2*ncon[e][i]-1,2*ncon[e][j]] += Kes[2*i-1,2*j]
					KS[2*ncon[e][i],2*ncon[e][j]-1] += Kes[2*i,2*j-1]
					KS[2*ncon[e][i],2*ncon[e][j]] += Kes[2*i,2*j]
				end
				if solveFor == "heat" || solveFor == "both"
					KH[ncon[e][i],ncon[e][j]] += Keh[i,j]
				end
			end
		end
	end

	return KS, fS, KH, fH
end;

@enum DegreeOfFreedom begin
	xDof = 1 << 0;
	yDof = 1 << 1;
	both = (1 << 0) | (1 << 1);
end

function handleDirichletBCsSolid!(constraints::Vector, dof::DegreeOfFreedom, value, K, f)
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

function handleDirichletBCsHeat!(constraints::Vector, value::Float64, K, f)
	nNodes = size(K)[1]
	for c = 1:length(constraints)
		idx = constraints[c]
		for i = 1:nNodes
			K[idx,i] = 0.0
		end
		K[idx,idx] = 1.0
		f[idx] = value
	end
end;

using LinearAlgebra
function handleNeumannBCsSolid!(constraints, pos, values, f)
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