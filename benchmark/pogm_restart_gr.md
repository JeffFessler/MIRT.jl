# POGM gradient restart benchmark

Run from the repository root (Julia standard libraries only):

```sh
julia --startup-file=no benchmark/pogm_restart_gr.jl
```

Optional arguments are sample count, row count, and column count:

```sh
julia --startup-file=no benchmark/pogm_restart_gr.jl 7 1000 80
```

The script loads solver source from this checkout and runs independently of the
CI benchmark suite. Synthetic fixtures are shared with the correctness tests.

## Method

- L1 logistic and half squared-hinge objectives: seed 17, lambda 0.015,
  dimensions 600 by 60, no intercept.
- A 10000-step original POGM run provides the reference; its mapping residual
  must be below 1e-10. Accuracy gates are residual <= 1e-7, absolute objective
  gap < 1e-8, and relative solution error < 1e-4.
- Compare original and isolated fixed-budget runs, offline calibrated fixed
  budgets, and early stopping checked every 1, 5, or 20 iterations. Offline
  calibration requires advance knowledge and is only a comparison baseline.
- Report medians of 5 warmed full calls with one BLAS thread. Initialization
  and history allocation are included; setup, reference computation,
  calibration, and diagnostics are excluded. GC runs before each sample.
- Allocated bytes measure full-call allocation volume, not peak memory or
  hot-loop-only allocations. A separate cheap quadratic probe reports
  incremental allocation per iteration, including gradient and history costs.

## Optional mutating callbacks

`f_grad!(dest, x)` and `g_prox!(dest, z, c)` can be supplied independently;
existing calls remain valid. Each callback fills its destination without
modifying inputs. The gradient buffer is reused, while proximal destinations
remain fresh to preserve callback history references. See the solver docstring
for the full interface and working-precision requirements.

The benchmark's mutating gradient uses `mul!` with preallocated scratch space.
Its gradient and prox callbacks each allocated zero bytes after warmup; the
complete solver still allocates iterate arrays and history.

## Results

Measured on Julia 1.13.3 on 2026-11-02 with the defaults above. All rows passed
the accuracy gates; ordinary and mutating callbacks stopped at the same iteration.

| Problem | Mode | Ordinary MiB | Mutating MiB | Ordinary ms | Mutating ms |
|---|---|---:|---:|---:|---:|
| Logistic | Fixed 2000 | 49.440 | 27.559 | 37.031 | 34.559 |
| Logistic | Check every 5 (45 iterations) | 1.282 | 0.681 | 1.037 | 0.991 |
| Squared hinge | Fixed 2000 | 50.116 | 28.235 | 20.123 | 16.627 |
| Squared hinge | Check every 5 (50 iterations) | 1.421 | 0.754 | 0.624 | 0.534 |

Early stopping avoids excess iterations but adds checking overhead. Mutating
callbacks reduce allocation when the caller supplies suitable implementations.
Timings depend on hardware and workload; these results do not establish a
universal speedup or an allocation-free solver.

## Validation

The complete `test/runtests.jl` passed on Julia 1.13.3, including Aqua (11/11),
original POGM tests (3/3), and the new solver tests (756/756). Tests cover original
trajectory agreement, composite stopping criteria, reference-solution accuracy,
optional mutating callbacks, and retained callback references. The documentation
site was not built.
