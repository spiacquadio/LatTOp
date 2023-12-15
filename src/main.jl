#=
	Currently this file acts as an input file for the "optimiseFunctions.jl" file.
	1) First run this file to load all the include functions
	2) Run the optimseFunctions.jl file
	3) Call the functions in the optimseFunctions.jl file from the terminal
=#

include("fileio.jl")
include("assembly.jl")
include("topologyOptimisation.jl")
include("optimiseFunctions.jl")
using Printf
using Plots
using Distributed
# Setup
quad = getQuadratureRules(Quad1);				# Get quadrature rules for 1st order quadrelaterals
shape = getShapeAndDerivatives(Quad1, quad);	# Get shape functions for 1st order quad

# Import mesh. (ncon = node connectivity, pos = nodal positions, boundaries = maps boundary name to boundary elements, boundaryNodes = maps boundary name to boundary nodes)
meshpath = "/Users/spiacqua/Documents/TESI/Topology_Optimisation/res/meshes/CombinedDomain_new_5mm.msh";
savepath = "/Users/spiacqua/Documents/TESI/Topology_Optimisation/res/results/";
ncon, pos, boundaries, boundaryNodes = importMsh(meshpath);

# Mesh stats
nElements = size(ncon[boundaries["Domain"]])[1];
nNodes = size(pos)[1];

# Initialise domain
K_T,f_T = getSparsityPatternHeat(ncon);		# Gets sparcity pattern of global matrix and forcing array		Heat problem
Ke_T = zeros(4,4,nElements);				# Initialise element matrices (used in sensitivity analysis)	Heat problem
K_S,f_S = getSparsityPatternSolid(ncon);	# Gets sparcity pattern of global matrix and forcing array		Solid problem
Ke_S = zeros(8,8,nElements);				# Initialise element matrices (used in sensitivity analysis)	Solid problem

HF = getWeightFactorMatrix(ncon[boundaries["Domain"]], pos, 6e-3);	# Weight factor matrix (filter matrix to filter sensitivities)
# neighbours = getNeighbours(ncon[boundaries["Domain"]]);					# Get element neighbours
elements = initElementMatrices(ncon[boundaries["Domain"]], pos, quad, shape);	# Precompute the element matrices

# Execute the program
# define unit cell and if lattices are used or if you should do a traditional topology optimisation
#useLattice=true;
#celltype = "bcc";#"f2ccz" or "bcc"
# optimal_weight=Float64[]; # initialise the vector for the optimal weights
# # Loop through different domain volume fractions and save the optimal weights, as well as the relative density field
# @distributed for i in range(1,stop=5,step=1)
# 	DomainVolumeFraction=0.1*i;
# 	therm_opt=optimiseTopologyThermal(DomainVolumeFraction,"thermal",false);
# 	solid_opt=optimiseTopologySolid(DomainVolumeFraction,"solid",false);
# 	pareto_opt=getParetoOptimalStructure((therm_opt[2],solid_opt[2]),DomainVolumeFraction);
# 	#multi_bcc=optimiseTopology(DomainVolumeFraction,pareto_opt[4],"multi_f2ccz_2_5mm_0$i",true);
# 	push!(optimal_weight,pareto_opt[4]);
# end