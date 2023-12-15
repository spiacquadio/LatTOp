using ForwardDiff
using ReverseDiff

include("../src/fileio.jl")
include("../src/shapeFunctions.jl")


quad = getQuadratureRules(Quad1);
shape = getShapeAndDerivatives(Quad1, quad);

ncon, pos, boundaries, boundaryNodes = importMsh("res/meshes/ThermalDomain.msh");
nNodes = size(pos)[1];
nElements = size(boundaries["Domain"])[1];
elements = initElementMatrices(ncon[boundaries["Domain"]], pos, quad, shape);
# Input flux 10-40kW