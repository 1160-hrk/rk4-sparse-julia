# Benchmark Results for Two-Level System Simulation (v0.0.2)

_Version: v0.0.2_  
_Date: 2025-07-09_

## Overview
This document summarizes the performance comparison between dense, sparse, and multi-threaded sparse matrix implementations of the Runge-Kutta 4th order method for simulating a two-level quantum system.

## Dense Matrix Implementation

### Performance Metrics
- **Number of trials**: 2,893
- **Execution time**:
  - Minimum: 1.451 ms
  - Mean: 1.727 ms ± 0.246 ms
  - Maximum: 3.637 ms
- **Memory usage**: 6,094.80 KiB
- **Allocation count**: 55,024

### Benchmark Details
```
BenchmarkTools.Trial: 2893 samples with 1 evaluation per sample.
Range (min … max):  1.451 ms …   3.637 ms  ┊ GC (min … max): 0.00% … 0.00%
Time  (median):     1.634 ms               ┊ GC (median):    0.00%
Time  (mean ± σ):   1.727 ms ± 245.594 μs  ┊ GC (mean ± σ):  5.22% ± 9.85%
```

## Sparse Matrix Implementation

### Performance Metrics
- **Number of trials**: 4,901
- **Execution time**:
  - Minimum: 0.820 ms
  - Mean: 1.019 ms ± 0.261 ms
  - Maximum: 2.445 ms
- **Memory usage**: 4,222.52 KiB
- **Allocation count**: 40,061

### Benchmark Details
```
BenchmarkTools.Trial: 4901 samples with 1 evaluation per sample.
Range (min … max):  820.083 μs …   2.445 ms  ┊ GC (min … max): 0.00% … 55.27%
Time  (median):     933.167 μs               ┊ GC (median):    0.00%
Time  (mean ± σ):     1.019 ms ± 260.877 μs  ┊ GC (mean ± σ):  8.25% ± 14.37%
```

## Multi-threaded Sparse Implementation

### Performance Metrics
- **Number of trials**: 409
- **Execution time**:
  - Minimum: 11.071 ms
  - Mean: 12.222 ms ± 0.782 ms
  - Maximum: 17.114 ms
- **Memory usage**: 20,160.56 KiB
- **Allocation count**: 325,074

### Benchmark Details
```
BenchmarkTools.Trial: 409 samples with 1 evaluation per sample.
Range (min … max):  11.071 ms …  17.114 ms  ┊ GC (min … max): 0.00% … 10.65%
Time  (median):     11.807 ms               ┊ GC (median):    0.00%
Time  (mean ± σ):   12.222 ms ± 782.261 μs  ┊ GC (mean ± σ):  4.95% ±  5.17%
```

## Performance Comparison

### Execution Time (Minimum)
1. Sparse: 0.820 ms (Best)
2. Dense: 1.451 ms (+77% vs Sparse)
3. Multi-threaded Sparse: 11.071 ms (+1250% vs Sparse)

### Memory Usage
1. Sparse: 4,222.52 KiB (Best)
2. Dense: 6,094.80 KiB (+44% vs Sparse)
3. Multi-threaded Sparse: 20,160.56 KiB (+377% vs Sparse)

### Allocation Count
1. Sparse: 40,061 (Best)
2. Dense: 55,024 (+37% vs Sparse)
3. Multi-threaded Sparse: 325,074 (+711% vs Sparse)

## System Information
- Julia version: 1.10.0
- OS: Linux 6.10.14-linuxkit
- Available threads: 1
- BLAS threads: 1
- Implementation: Runge-Kutta 4th order method
- Test case: Two-level quantum system with Rabi oscillations

## Analysis

1. **Sparse vs Dense Implementation**:
   - Sparse implementation maintains its superior performance
   - Approximately 44% faster execution time
   - 31% less memory usage
   - 27% fewer allocations

2. **Multi-threaded Implementation**:
   - Currently shows degraded performance
   - Main factors:
     - Single thread environment (threads = 1)
     - Overhead from thread management
     - Increased memory allocations
   - Potential improvements:
     - Test in multi-core environment
     - Optimize thread synchronization
     - Reduce memory allocations

## Recommendations

1. **For Single-Thread Environment**:
   - Use the sparse implementation
   - Avoid multi-threaded version

2. **Future Optimizations**:
   - Test multi-threaded implementation in multi-core environment
   - Consider GPU implementation for larger systems
   - Optimize memory allocations in multi-threaded version
   - Implement SIMD optimizations 