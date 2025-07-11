using LinearAlgebra
# パッケージのインストールと読み込み
using Pkg
Pkg.add("BenchmarkTools")
using BenchmarkTools

include("rk_schrodinger.jl") 
include("plots_config.jl")
using .PlotsConfig

# -------------------------------------------
# 二準位系の定義
# -------------------------------------------

hbar = 1.0  # プランク定数
omega = 1.0  # 遷移エネルギー差
mu = 0.1  # 双極子結合強度

dim = 2

# ハミルトニアン（対角成分のみ）
H0 = Diagonal(ComplexF64[0.0, hbar * omega])

# 双極子行列（ここでは x 成分だけを考える）
mux = ComplexF64[0 1; 1 0] * mu
muy = zeros(ComplexF64, dim, dim)  # y 成分なし

# -------------------------------------------
# 電場設定
# -------------------------------------------

T = 100.0  # 全体時間
dt = 0.01
t = 0:dt:T
steps = length(t)

E0 = 1.0  # 電場振幅
omega_L = omega  # 電場周波数（共鳴）

Ex = E0 .* sin.(omega_L .* t)
Ey = zeros(length(t))  # y 成分なし

# -------------------------------------------
# 初期状態
# -------------------------------------------

psi0 = ComplexF64[1.0, 0.0]  # 基底状態

# -------------------------------------------
# 計算
# -------------------------------------------

stride = 1  # サンプリング間隔

using SparseArrays  # sparse行列のサポート追加

# 行列をsparse形式に変換
H0_sparse = sparse(H0)
mux_sparse = sparse(mux)
muy_sparse = sparse(muy)

traj = rk4_cpu(
    Matrix(H0), Matrix(mux), Matrix(muy),
    Ex, Ey, psi0, dt*2;  # stride を使用
    return_traj=true, stride=stride, renorm=false
)

# traj = rk4_cpu_sparse(
#     H0_sparse, mux_sparse, muy_sparse, Ex, Ey, psi0, dt*2;
#     return_traj=true, stride=stride, renorm=false
# )


# -------------------------------------------
# ラビ振動の解析解
# -------------------------------------------

# ラビ周波数（共鳴条件下）
Omega_R = mu * E0

# 解析解の計算
t_analytical = t[1:stride*2:end]
pop1_analytical = cos.(Omega_R * t_analytical / 2).^2
pop2_analytical = sin.(Omega_R * t_analytical / 2).^2

# -------------------------------------------
# 結果のプロット
# -------------------------------------------

using Plots

pop1 = abs2.(traj[1, :])  # 状態 |0⟩ の確率
pop2 = abs2.(traj[2, :])  # 状態 |1⟩ の確率

p = plot(t[1:stride*2:end], pop1, label="Ground state (numerical)", lw=2)
plot!(p, t[1:stride*2:end], pop2, label="Excited state (numerical)", lw=2)
plot!(p, t_analytical, pop1_analytical, label="Ground state (analytical)", 
      ls=:dash, lw=2)
plot!(p, t_analytical, pop2_analytical, label="Excited state (analytical)", 
      ls=:dash, lw=2)
xlabel!("Time")
ylabel!("Population")
title!("Two-level system: Rabi oscillations")
display(p)  # プロットを表示

# -------------------------------------------
# 計算速度の測定
# -------------------------------------------

# 行列を事前に変換
H0_dense = Matrix(H0)
mux_dense = Matrix(mux)
muy_dense = Matrix(muy)

# 単一実行の時間計測
println("\n計算時間の測定:")
@time begin
    rk4_cpu(
        H0_dense, mux_dense, muy_dense,
        Ex, Ey, psi0, dt*2;
        return_traj=true, stride=stride, renorm=false
    )
end

# 詳細なベンチマーク
println("\n詳細なベンチマーク:")
bench = @benchmark rk4_cpu(
    $H0_dense, $mux_dense, $muy_dense,
    $Ex, $Ey, $psi0, $dt*2;
    return_traj=true, stride=$stride, renorm=false
)
display(bench)

# 統計情報の表示
println("\n統計情報:")
println("試行回数: ", length(bench.times), " 回")
println("最小実行時間: ", minimum(bench.times) / 1e6, " ms")
println("平均実行時間: ", mean(bench.times) / 1e6, " ms")
println("最大実行時間: ", maximum(bench.times) / 1e6, " ms")
println("標準偏差: ", std(bench.times) / 1e6, " ms")
println("メモリアロケーション: ", bench.memory / 1024, " KiB")
println("アロケーション回数: ", bench.allocs, " 回")

# スパース版の計算速度の測定
println("\nスパース版の計算時間の測定:")
@time begin
    rk4_cpu_sparse(
        H0_sparse, mux_sparse, muy_sparse,
        Ex, Ey, psi0, dt*2;
        return_traj=true, stride=stride, renorm=false
    )
end

# スパース版の詳細なベンチマーク
println("\nスパース版の詳細なベンチマーク:")
bench_sparse = @benchmark rk4_cpu_sparse(
    $H0_sparse, $mux_sparse, $muy_sparse,
    $Ex, $Ey, $psi0, $dt*2;
    return_traj=true, stride=$stride, renorm=false
)
display(bench_sparse)

# スパース版の統計情報の表示
println("\nスパース版の統計情報:")
println("試行回数: ", length(bench_sparse.times), " 回")
println("最小実行時間: ", minimum(bench_sparse.times) / 1e6, " ms")
println("平均実行時間: ", mean(bench_sparse.times) / 1e6, " ms")
println("最大実行時間: ", maximum(bench_sparse.times) / 1e6, " ms")
println("標準偏差: ", std(bench_sparse.times) / 1e6, " ms")
println("メモリアロケーション: ", bench_sparse.memory / 1024, " KiB")
println("アロケーション回数: ", bench_sparse.allocs, " 回")

# マルチスレッド版の計算速度の測定
println("\nマルチスレッド版の計算時間の測定:")
BLAS.set_num_threads(Threads.nthreads())  # BLASの並列化を有効化
@time begin
    rk4_cpu_sparse_mt(
        H0_sparse, mux_sparse, muy_sparse,
        Ex, Ey, psi0, dt*2;
        return_traj=true, stride=stride, renorm=false
    )
end

# マルチスレッド版の詳細なベンチマーク
println("\nマルチスレッド版の詳細なベンチマーク:")
bench_mt = @benchmark rk4_cpu_sparse_mt(
    $H0_sparse, $mux_sparse, $muy_sparse,
    $Ex, $Ey, $psi0, $dt*2;
    return_traj=true, stride=$stride, renorm=false
)
display(bench_mt)

# マルチスレッド版の統計情報の表示
println("\nマルチスレッド版の統計情報:")
println("試行回数: ", length(bench_mt.times), " 回")
println("最小実行時間: ", minimum(bench_mt.times) / 1e6, " ms")
println("平均実行時間: ", mean(bench_mt.times) / 1e6, " ms")
println("最大実行時間: ", maximum(bench_mt.times) / 1e6, " ms")
println("標準偏差: ", std(bench_mt.times) / 1e6, " ms")
println("メモリアロケーション: ", bench_mt.memory / 1024, " KiB")
println("アロケーション回数: ", bench_mt.allocs, " 回")

# SIMD最適化版の計算速度の測定
println("\nSIMD最適化版の計算時間の測定:")
@time begin
    rk4_cpu_sparse_simd(
        H0_sparse, mux_sparse, muy_sparse,
        Ex, Ey, psi0, dt*2;
        return_traj=true, stride=stride, renorm=false
    )
end

# SIMD最適化版の詳細なベンチマーク
println("\nSIMD最適化版の詳細なベンチマーク:")
bench_simd = @benchmark rk4_cpu_sparse_simd(
    $H0_sparse, $mux_sparse, $muy_sparse,
    $Ex, $Ey, $psi0, $dt*2;
    return_traj=true, stride=$stride, renorm=false
)
display(bench_simd)

# SIMD最適化版の統計情報の表示
println("\nSIMD最適化版の統計情報:")
println("試行回数: ", length(bench_simd.times), " 回")
println("最小実行時間: ", minimum(bench_simd.times) / 1e6, " ms")
println("平均実行時間: ", mean(bench_simd.times) / 1e6, " ms")
println("最大実行時間: ", maximum(bench_simd.times) / 1e6, " ms")
println("標準偏差: ", std(bench_simd.times) / 1e6, " ms")
println("メモリアロケーション: ", bench_simd.memory / 1024, " KiB")
println("アロケーション回数: ", bench_simd.allocs, " 回")

# システム情報の表示
println("\nシステム情報:")
println("利用可能なスレッド数: ", Threads.nthreads())
println("BLASスレッド数: ", BLAS.get_num_threads())
println("SIMD命令セット: ", LoopVectorization.VectorizationBase.register_size())
