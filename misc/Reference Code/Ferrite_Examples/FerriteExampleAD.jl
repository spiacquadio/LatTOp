using Ferrite, SparseArrays
using ForwardDiff
using ReverseDiff
using Zygote

grid = generate_grid(Quadrilateral, (10, 10));

dim = 2
ip = Lagrange{dim, RefCube, 1}()
qr = QuadratureRule{dim, RefCube}(2)
cellvalues = CellScalarValues(qr, ip);

dh = DofHandler(grid)
add!(dh, :u, 1)
close!(dh);


ch = ConstraintHandler(dh);

∂Ω = union(
    getfaceset(grid, "left"),
    getfaceset(grid, "right"),
    getfaceset(grid, "top"),
    getfaceset(grid, "bottom"),
);

dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
add!(ch, dbc);

close!(ch)

function assemble_element!(k, Ke, fe, cellvalues::CellScalarValues)
    n_basefuncs = getnbasefunctions(cellvalues)
    for i = 1:n_basefuncs
        fe[i] = 0.0
        for j = 1:n_basefuncs
            Ke[i,j] = 0.0
        end
    end

    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            fe[i] += δu * dΩ
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end 
    return Ke, fe
end
function assemble_element(k, Ke, fe, cellvalues::CellScalarValues)
    KE = Zygote.Buffer(Ke)
    FE = Zygore.Buffer(fe)
    assemble_element!(k, KE, FE, cellvalues)
    return copy(KE), copy(FE)
end

function assemble_global(k, cellvalues::CellScalarValues, K::SparseMatrixCSC, dh::DofHandler)
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = Matrix{Float64}(undef, n_basefuncs,n_basefuncs)
	fe = Array{Float64}(undef, n_basefuncs)

    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        assemble_element!(k, Ke, fe, cellvalues)
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end

function get_matrices(k)
    K = create_sparsity_pattern(dh)
	K, f = assemble_global(k, cellvalues, K, dh)
    A = Zygote.Buffer(K)
    F = Zygote.Buffer(f)
    apply!(A, F, ch)
    return copy(A), copy(F)
end

function solveCompliance(k)
    K, f = get_matrices(k)
	u = K \ f;
	return u' * f
end

g = (k) -> ForwardDiff.derivative(solveCompliance, k)
g_ = (k) -> ReverseDiff.gradient(solveCompliance, [k])
gz = (k) -> Zygote.gradient(solveCompliance, k)

# vtk_grid("heat_equation_g", dh) do vtk
#     vtk_point_data(vtk, dh, u)
# end