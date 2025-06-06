# Below are several variants of the same function.
# TODO: refactor with multiple dispatch?

function assemble_matrices_comprehension(green_functions, mesh, wavenumber; direct=true)
    # Use comprehensions to build S and D matrices
    S = @inbounds [-1/2τ̅ * Complex(integral(green_functions, element(mesh, i), element(mesh, j), wavenumber)) for i in 1:mesh.nfaces, j in 1:mesh.nfaces]

    D = @inbounds [begin
            element_i = element(mesh, i)
            element_j = element(mesh, j)

            # Select the normal based on direct flag
            n = direct ? normal(element_j) : normal(element_i)

            -1/2τ̅ * Complex(n' * integral_gradient(green_functions, element_i, element_j, wavenumber; with_respect_to_first_variable=!direct))
        end for i in 1:mesh.nfaces, j in 1:mesh.nfaces]

    # Add diagonal elements to D
    @inbounds begin
        @views begin
            D1 = D .+ Diagonal(0.5 .* I(mesh.nfaces))
        end
    end

    return S, D1
end

function assemble_matrices_explicit_both(green_functions, mesh, wavenumber; direct=true)
    # Variant of the above using `both_integral_and_integral_gradient`
    # Probably not differentiable with Zygote
    S = Array{ComplexF64, 2}(undef, (mesh.nfaces, mesh.nfaces))
    D = Array{ComplexF64, 2}(undef, (mesh.nfaces, mesh.nfaces))
    for i in 1:mesh.nfaces
        for j in 1:mesh.nfaces
            element_i = element(mesh, i)
            element_j = element(mesh, j)
            Sij, Dij = both_integral_and_integral_gradient(green_functions, element_i, element_j, wavenumber; with_respect_to_first_variable=!direct)
            S[i, j] = -1/4π * Complex(Sij)
            n = direct ? normal(element_j) : normal(element_i)
            D[i, j] = -1/4π * Complex(n' * Dij)
        end
        D[i, i] = D[i, i] + 0.5
    end
    return S, D
end

function assemble_matrices_broadcasting(greens_functions, mesh, wavenumber; direct=true, arrtype=Array)
    # Another variant, that works with e.g. CuArray
    elements = arrtype([element(mesh, i) for i in 1:mesh.nfaces])
    co_elements = reshape(elements, (1, length(elements)))  # Same but as a (1, n) row vector
    S(e1, e2) = integral(greens_functions, e1, e2, wavenumber)
    S_matrix = -1/(4π) * S.(elements, co_elements)

    if direct
        D(e1, e2) = MarineHydro.normal(e2)' * integral_gradient(greens_functions, e1, e2, wavenumber, with_respect_to_first_variable=false)
        D_matrix = -1/(4π) * D.(elements, co_elements)
        return S_matrix, D_matrix + arrtype(0.5 .* I(mesh.nfaces))
    else
        K(e1, e2) = MarineHydro.normal(e1)' * integral_gradient(greens_functions, e1, e2, wavenumber, with_respect_to_first_variable=true)
        K_matrix = -1/(4π) * K.(elements, co_elements)
        return S_matrix, K_matrix + arrtype(0.5 .* I(mesh.nfaces))
    end
end


"""
    assemble_matrices(green_functions, mesh, wavenumber; direct=true)

Assembles the influence matrices based on the tuple of provided Green's functions, mesh, and wavenumber.

# Arguments
- `green_functions`: Iterable of `GreensFunction` objects.
- `mesh`: Floating body mesh with panel information such as vertices, faces, normals, areas etc.
- `wavenumber`: Incoming ocean wavenumber
- `direct=true`: A flag to specify whether to use direct BEM vs Indirect BEM.

# Returns
- A tuple of assembled matrices. S and (D or K) depending on the flag.
"""
# The default method:
const assemble_matrices = assemble_matrices_comprehension


function assemble_matrix_wu(mesh, wavenumber; direct=true)
    return assemble_matrices([Rankine(), RankineReflected(), GFWu()], mesh, wavenumber; direct)
end


function solve(D, S, bc; direct::Bool=true)
    if direct
        ϕ = implicit_linear(D,S*bc)
    else
        K = D
        sources = implicit_linear(K,bc)
        ϕ = S * sources
    end
    return ϕ
end

