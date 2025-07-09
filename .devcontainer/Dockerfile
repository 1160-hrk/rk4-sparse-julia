# ベースイメージとしてJulia 1.10.0を使用
FROM julia:1.10.0

# 作業ディレクトリを設定
WORKDIR /workspace

# 必要なパッケージを事前にインストール
RUN julia -e 'using Pkg; \
    Pkg.add([ \
        "BenchmarkTools", \
        "Plots", \
        "LoopVectorization" \
    ]); \
    Pkg.precompile()'

# Project.tomlとManifest.tomlをコピー
COPY Project.toml .

# プロジェクトの依存関係を解決
RUN julia --project -e 'using Pkg; Pkg.resolve()'

# ソースコードをコピー
COPY src/ src/

# コンテナ起動時のコマンド
CMD ["julia", "--project"]
