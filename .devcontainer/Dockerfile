# ベースイメージとしてJulia 1.10.0を使用
FROM julia:1.10.0

# 必要なシステムパッケージをインストール
RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

# 作業ディレクトリを設定
WORKDIR /workspace

# 必要なパッケージを事前にインストール
RUN julia -e 'using Pkg; \
    Pkg.add([ \
        "BenchmarkTools", \
        "Plots", \
        "LoopVectorization", \
        "SparseArrays" \
    ]); \
    Pkg.precompile()'

# Project.tomlとManifest.tomlをコピー
COPY Project.toml .

# プロジェクトの依存関係を解決
RUN julia --project -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'

# ソースコードをコピー
COPY src/ src/

# コンテナ起動時のコマンド
CMD ["julia", "--project"]
