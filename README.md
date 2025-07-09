# Two-Level System Simulation

This project implements and benchmarks various numerical methods for simulating a two-level quantum system.

## Requirements

- Docker
- VS Code with Remote - Containers extension

## Quick Start

1. Clone the repository:
```bash
git clone <repository-url>
cd <repository-name>
```

2. Open in VS Code with Dev Containers:
- Open VS Code
- Press F1 and select "Dev Containers: Open Folder in Container..."
- Select the project directory

3. Run the example:
```julia
julia --project
include("src/example_twolevelstate.jl")
```

## Project Structure

```
.
├── .devcontainer/      # Dev container configuration
├── src/                # Source code
│   ├── rk_schrodinger.jl          # Runge-Kutta implementations
│   └── example_twolevelstate.jl    # Example usage and benchmarks
├── Dockerfile          # Docker configuration
├── Project.toml        # Julia package dependencies
└── README.md          # This file
```

## Implementation Details

The project includes several implementations of the Runge-Kutta 4th order method:
1. Dense matrix implementation
2. Sparse matrix implementation
3. Multi-threaded sparse implementation
4. SIMD-optimized sparse implementation

## Benchmarking

Run the example script to see performance comparisons between different implementations:
```julia
julia --project src/example_twolevelstate.jl
```

## License

MIT 