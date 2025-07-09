# Benchmark Results for Two-Level System Simulation (v0.0.1)

_Version: v0.0.1_  
_Date: 2025-07-09_

## Overview
This document summarizes the performance comparison between dense and sparse matrix implementations of the Runge-Kutta 4th order method for simulating a two-level quantum system.

## Dense Matrix Implementation

### Performance Metrics
- **Number of trials**: 2,907
- **Execution time**:
  - Minimum: 1.436 ms
  - Mean: 1.717 ms ± 0.276 ms
  - Maximum: 9.571 ms
- **Memory usage**: 6,094.80 KiB
- **Allocation count**: 55,024

### Benchmark Details
```
BenchmarkTools.Trial: 2907 samples with 1 evaluation per sample.
Range (min … max):  1.436 ms …   9.571 ms  ┊ GC (min … max): 0.00% … 0.00%
Time  (median):     1.630 ms               ┊ GC (median):    0.00%
Time  (mean ± σ):   1.717 ms ± 276.149 μs  ┊ GC (mean ± σ):  5.15% ± 9.91%
```

## Sparse Matrix Implementation

### Performance Metrics
- **Number of trials**: 4,864
- **Execution time**:
  - Minimum: 0.824 ms
  - Mean: 1.026 ms ± 0.288 ms
  - Maximum: 6.965 ms
- **Memory usage**: 4,222.52 KiB
- **Allocation count**: 40,061

### Benchmark Details
```
BenchmarkTools.Trial: 4864 samples with 1 evaluation per sample.
Range (min … max):  824.333 μs …   6.965 ms  ┊ GC (min … max): 0.00% …  0.00%
Time  (median):     939.084 μs               ┊ GC (median):    0.00%
Time  (mean ± σ):     1.026 ms ± 287.725 μs  ┊ GC (mean ± σ):  7.97% ± 14.11%
```

## Performance Comparison

### Execution Time
- Sparse implementation is approximately **40% faster** than dense implementation (based on minimum execution time)
- Mean execution time improved by approximately **40%**

### Memory Usage
- Sparse implementation uses approximately **31% less memory**
- Allocation count reduced by approximately **27%**

## Conclusions
1. The sparse matrix implementation shows significant improvements in both execution time and memory usage
2. Most notable improvements:
   - Execution time reduced by ~40%
   - Memory usage reduced by ~31%
   - Fewer memory allocations
3. The sparse implementation shows more consistent performance with lower maximum execution times

## System Information
- Julia version: 1.10.0
- OS: Linux 6.10.14-linuxkit
- Implementation: Runge-Kutta 4th order method
- Test case: Two-level quantum system with Rabi oscillations 