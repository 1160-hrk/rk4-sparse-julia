# ベースイメージとしてJulia 1.10.0を使用
FROM julia:1.10.0

# 作業ディレクトリを設定
WORKDIR /workspace

# Project.tomlとManifest.tomlをコピー
COPY Project.toml .
COPY Manifest.toml* .

# パッケージのインストールと事前コンパイル
RUN julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# ソースコードをコピー
COPY src/ src/

# コンテナ起動時のコマンド
CMD ["julia", "--project"] 