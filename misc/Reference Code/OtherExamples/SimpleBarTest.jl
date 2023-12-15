using Gridap
using GridapGmsh
using ForwardDiff
using AutoGrad

# Load Model
model = GmshDiscreteModel("meshes/SimpleBar.msh")
order = 1
reffe = ReferenceFE(lagrangian,Float64,order)

# Define Function Spaces
V = TestFESpace(model,reffe,dirichlet_tags=["Left","Right"])
U = TrialFESpace(V,[0,1])

# Define Domain and Boundary
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)

# Define PDE
function solvePDE(k)
	a(u,v) = ∫(k[1]* ∇(v)⋅∇(u) )dΩ
	l(v) = 0
	op = AffineFEOperator(a,l,U,V)
	ls = LUSolver()
	solver = LinearFESolver(ls)
	u = solve(solver,op)
end

g = (k) -> ForwardDiff.jacobian(solvePDE, [k])

x = AutoGrad.Param([2.0])
g = AutoGrad.@diff solvePDE(x)
