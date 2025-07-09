using LinearAlgebra
using SparseArrays

function rk4_cpu(
    H0::Matrix{ComplexF64},
    mux::Matrix{ComplexF64},
    muy::Matrix{ComplexF64},
    Ex::Vector{Float64},
    Ey::Vector{Float64},
    psi0::Vector{ComplexF64},
    dt::Float64;
    return_traj::Bool = true,
    stride::Int = 1,
    renorm::Bool = false
)
    steps = (length(Ex) - 1) ÷ 2
    Ex3 = hcat(Ex[1:2:end-2], Ex[2:2:end-1], Ex[3:2:end])
    Ey3 = hcat(Ey[1:2:end-2], Ey[2:2:end-1], Ey[3:2:end])

    psi = copy(psi0)
    dim = length(psi)
    n_out = steps ÷ stride + 1
    out = Matrix{ComplexF64}(undef, dim, n_out)
    out[:, 1] .= psi
    idx = 2

    buf = similar(psi)
    k1 = similar(psi)
    k2 = similar(psi)
    k3 = similar(psi)
    k4 = similar(psi)

    for s in 1:steps
        ex1, ex2, ex4 = Ex3[s, :]
        ey1, ey2, ey4 = Ey3[s, :]

        H1 = H0 .+ mux .* ex1 .+ muy .* ey1
        H2 = H0 .+ mux .* ex2 .+ muy .* ey2
        H4 = H0 .+ mux .* ex4 .+ muy .* ey4

        k1 .= -im .* (H1 * psi)
        buf .= psi .+ 0.5 * dt .* k1
        k2 .= -im .* (H2 * buf)
        buf .= psi .+ 0.5 * dt .* k2
        k3 .= -im .* (H2 * buf)
        buf .= psi .+ dt .* k3
        k4 .= -im .* (H4 * buf)

        psi .+= (dt / 6.0) .* (k1 .+ 2k2 .+ 2k3 .+ k4)

        if renorm
            psi ./= sqrt(real(psi' * psi))
        end

        if return_traj && (s % stride == 0)
            out[:, idx] .= psi
            idx += 1
        end
    end

    return return_traj ? out : reshape(psi, :, 1)
end


function rk4_cpu_sparse(
    H0::SparseMatrixCSC{ComplexF64,Int},
    mux::SparseMatrixCSC{ComplexF64,Int},
    muy::SparseMatrixCSC{ComplexF64,Int},
    Ex::Vector{Float64},
    Ey::Vector{Float64},
    psi0::Vector{ComplexF64},
    dt::Float64;
    return_traj::Bool = true,
    stride::Int = 1,
    renorm::Bool = false
)
    steps = (length(Ex) - 1) ÷ 2
    Ex3 = hcat(Ex[1:2:end-2], Ex[2:2:end-1], Ex[3:2:end])
    Ey3 = hcat(Ey[1:2:end-2], Ey[2:2:end-1], Ey[3:2:end])

    psi = copy(psi0)
    dim = length(psi)
    n_out = steps ÷ stride + 1
    out = Matrix{ComplexF64}(undef, dim, n_out)
    out[:, 1] .= psi
    idx = 2

    # バッファの初期化
    buf = similar(psi)
    k1 = similar(psi)
    k2 = similar(psi)
    k3 = similar(psi)
    k4 = similar(psi)

    # スパース行列の準備
    rows = Vector{Int}()
    cols = Vector{Int}()
    vals = Vector{ComplexF64}()
    
    # 非ゼロ要素のパターンを収集
    for (i, j, v) in zip(findnz(H0)..., findnz(mux)..., findnz(muy)...)
        if !iszero(v)
            push!(rows, i)
            push!(cols, j)
        end
    end
    
    # 重複を除去
    unique_indices = unique(tuple.(rows, cols))
    rows = first.(unique_indices)
    cols = last.(unique_indices)
    
    # 各行列から値を抽出
    H0_data = [H0[i,j] for (i,j) in unique_indices]
    mux_data = [mux[i,j] for (i,j) in unique_indices]
    muy_data = [muy[i,j] for (i,j) in unique_indices]
    
    # 初期スパース行列の作成
    H = sparse(rows, cols, zeros(ComplexF64, length(rows)), dim, dim)

    for s in 1:steps
        ex1, ex2, ex4 = Ex3[s, :]
        ey1, ey2, ey4 = Ey3[s, :]

        # 行列要素の更新
        vals = H0_data .+ mux_data .* ex1 .+ muy_data .* ey1
        H = sparse(rows, cols, vals, dim, dim)
        k1 .= -im .* (H * psi)
        buf .= psi .+ 0.5 * dt .* k1

        vals = H0_data .+ mux_data .* ex2 .+ muy_data .* ey2
        H = sparse(rows, cols, vals, dim, dim)
        k2 .= -im .* (H * buf)
        buf .= psi .+ 0.5 * dt .* k2

        k3 .= -im .* (H * buf)
        buf .= psi .+ dt .* k3

        vals = H0_data .+ mux_data .* ex4 .+ muy_data .* ey4
        H = sparse(rows, cols, vals, dim, dim)
        k4 .= -im .* (H * buf)

        psi .+= (dt / 6.0) .* (k1 .+ 2k2 .+ 2k3 .+ k4)

        if renorm
            psi ./= sqrt(real(psi' * psi))
        end

        if return_traj && (s % stride == 0)
            out[:, idx] .= psi
            idx += 1
        end
    end

    return return_traj ? out : reshape(psi, :, 1)
end
