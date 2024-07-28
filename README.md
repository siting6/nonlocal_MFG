# Nonlocal MFG

This repository contains MATLAB code for computing nonlocal mean field games (MFG).

## Getting Started

To test the moving obstacle examples described in [this paper](https://epubs.siam.org/doi/abs/10.1137/20M1334668), use the `fun_run_mfg_examples` function.

### Usage Instructions

1. **Windows Users:**
   - The provided `mex` functions should work out of the box.
   
2. **Non-Windows Users:**
   - If you encounter issues with the `mex` functions, use the alternative functions without the `mex` suffix.

### Optimization Details

- The program is optimized using MATLAB Coder, which generates the `mex` files to enhance performance.