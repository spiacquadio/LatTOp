
using Distributed
#=
	This struct stores the finite element matrices which can be precomputed
		for efficiency. Their basic forms are computed and the material
		properties can be multiplied with these matrices in the assembly stage
=#
mutable struct Element
	Me::AbstractMatrix{Float64};	# Mass matrix (Me = int N'*N dOmega)
	Ke::AbstractMatrix{Float64};	# Stiffness matrix (Ke = int gradN*gradN' dOmega)
	fe::AbstractVector{Float64};	# Forcing array (fe = int N dOmega)

	Element(size::Int64) = new(zeros(size,size), zeros(size,size), zeros(size));
end;

#=
	This function initialises the finite element matrices.
=#
@everywhere function initElementMatrices(ncon, coords, quad, shape)
	elementMatrices = Vector{Element}(undef, nElements);
	@distributed for e = 1:nElements
		Me = zeros(4,4);
		Ke = zeros(4,4);
		fe = zeros(4);

		elePoints = coords[ncon[e],:]

		# Loop over quadrature points
		@distributed for q = 1:quad.nQuadPoints
			detJ, dNdX = getJacobian(q, shape, elePoints);	# Gets the determinant and shape function deriviatives
			dNdX1 = dNdX[1,:];	# For easy access
			dNdX2 = dNdX[2,:];	# For easy access
			integrator = quad.weights[q] * detJ;	# Integration constant

			# Loop over element nodes
			@distributed for i = 1:4
				fe[i] += integrator * shape.N[q,i];
				@distributed for j = 1:4
					Me[i,j] += integrator * shape.N[q,i] * shape.N[q,j];
					Ke[i,j] += integrator * (dNdX1[i] * dNdX1[j] + dNdX2[i] * dNdX2[j]);
				end
			end
		end

		# Assign matrices
		elementMatrices[e] = Element(4);
		elementMatrices[e].Me = Me;
		elementMatrices[e].Ke = Ke;
		elementMatrices[e].fe = fe;
	end

	return elementMatrices;
end;


# Struct to store quadrature rules
struct QuadratureRules
	nQuadPoints::Int16
	weights::Vector
	quadPoints::Matrix
end;

# Enum for selecting type of element
@enum ElementType begin
	Quad1 = 1;
	Lin1 = 2;
end;

# Gauss weughts and positions
@everywhere function getQuadratureRules(type::ElementType)
	if type == Quad1
		QuadratureRules(4, [1.0,1.0,1.0,1.0], 1/sqrt(3.0) * [-1.0;-1.0;;  1.0;-1.0;;  1.0;1.0;;  -1.0;1.0])
	elseif type == Lin1
		QuadratureRules(1, [2.0], [0.0;;])
	end
end;

# Struct to store the shape functions and derivatives
struct ShapeFunctions1D
	N::Matrix
	dNdXi::Matrix
end;
struct ShapeFunctions2D
	nNodes::Int32
	N::Matrix
	dNdXi1::Matrix
	dNdXi2::Matrix
end;
struct ShapeFunctions3D
	nNodes::Int32
	N::Matrix
	dNdXi1::Matrix
	dNdXi2::Matrix
	dNdXi3::Matrix
end;

# Gets the shape functions of the selected element
@everywhere function getShapeFunctions(type::ElementType, quad::QuadratureRules)
	if type == Quad1
		xi1 = quad.quadPoints[1,:]
		xi2 = quad.quadPoints[2,:]
		shapeFuncs = [	
			(1.0 .- xi1) .* (1.0 .- xi2)/4.0;; 
			(1.0 .+ xi1) .* (1.0 .- xi2)/4.0;;
			(1.0 .+ xi1) .* (1.0 .+ xi2)/4.0;;
			(1.0 .- xi1) .* (1.0 .+ xi2)/4.0]
		return shapeFuncs
	elseif type == Lin1
		xi = quad.quadPoints
		shapeFuncs = 0.5 * [1.0 .- xi;; xi .- 1.0]
		return shapeFuncs
	end
end;

# Gets the shape function gradients
@everywhere function getShapeDXi1(quad::QuadratureRules)
	xi2 = quad.quadPoints[2,:]
	shapeDxi1 = [
		(xi2 .- 1.0)/4.0;;
		(1.0 .- xi2)/4.0;;
		(1.0 .+ xi2)/4.0;;
	   -(1.0 .+ xi2)/4.0]
	return shapeDxi1
end;
@everywhere function getShapeDXi2(quad::QuadratureRules)
	xi1 = quad.quadPoints[1,:]
	shapeDxi1 = [
		(xi1 .- 1.0)/4.0;;
	   -(1.0 .+ xi1)/4.0;;
		(1.0 .+ xi1)/4.0;;
		(1.0 .- xi1)/4.0]
	return shapeDxi1
end;
@everywhere function getShapeDxi(quad::QuadratureRules)
	dNdXi = [-0.5;; 0.5]
	return dNdXi
end;
@everywhere function getBmatrix!(B::Matrix, dNdX::Vector, dNdY::Vector)
	for i = 1:4
		B[1,2*i-1]	= dNdX[i];
		B[2,2*i]	= dNdY[i];
		B[3,2*i-1]	= dNdY[i];
		B[3,2*i]	= dNdX[i];
	end
end


@everywhere function getShapeAndDerivatives(type::ElementType, quad::QuadratureRules)
	if type == Quad1
		ShapeFunctions2D(4, getShapeFunctions(type, quad), getShapeDXi1(quad), getShapeDXi2(quad))
	elseif type == Lin1
		ShapeFunctions1D(getShapeFunctions(type, quad), getShapeDxi(quad))
	end
end;

# Gets the jacobian, value of gradients at a gauss point
@everywhere function getJacobian(iQuad, shapeFuncs::ShapeFunctions2D, elePoints)
	dNdXi1 = shapeFuncs.dNdXi1[iQuad,:]'
	dNdXi2 = shapeFuncs.dNdXi2[iQuad,:]'
	jacobian = [
		dNdXi1 * elePoints[:,1]; dNdXi1 * elePoints[:,2];;
		dNdXi2 * elePoints[:,1]; dNdXi2 * elePoints[:,2]]

	detJ = (jacobian[1,1] * jacobian[2,2]) - (jacobian[2,1] * jacobian[1,2])
	dNdX = jacobian \ [dNdXi1; dNdXi2]
	return detJ, dNdX
end;