# Cookbook

Common tasks and examples for **SPARC-atomSFE** (`from atom import AtomicDFTSolver`).

`AtomicDFTSolver.solve()` returns a **dictionary** of arrays and metadata (see keys such as `energy`, `rho`, `converged`, `iterations`, `quadrature_nodes`, …).

## Basic tasks

### Calculate energy for a single atom

```python
from atom import AtomicDFTSolver

solver = AtomicDFTSolver(
    atomic_number=13,
    xc_functional="GGA_PBE",
)
result = solver.solve()
print(f"Total energy: {result['energy']:.6f} Ha")
```

### Compare different XC functionals

```python
from atom import AtomicDFTSolver

functionals = ["LDA_PW", "GGA_PBE", "SCAN"]
energies = {}

for func in functionals:
    solver = AtomicDFTSolver(atomic_number=1, xc_functional=func)
    result = solver.solve()
    energies[func] = result["energy"]
    print(f"{func}: {result['energy']:.6f} Ha")
```

### Electron density on the quadrature grid

```python
from atom import AtomicDFTSolver

solver = AtomicDFTSolver(atomic_number=1, xc_functional="LDA_PW")
result = solver.solve()

r = result["quadrature_nodes"]
rho = result["rho"]
print(f"Grid points: {len(r)}")
print(f"Density range: [{rho.min():.6e}, {rho.max():.6e}]")
```

## Visualization

### Plot electron density

```python
import matplotlib.pyplot as plt
from atom import AtomicDFTSolver

solver = AtomicDFTSolver(atomic_number=1, xc_functional="LDA_PW")
result = solver.solve()

r = result["quadrature_nodes"]
rho = result["rho"]

plt.figure(figsize=(8, 6))
plt.plot(r, rho, "b-", linewidth=2)
plt.xlabel("Radius (Bohr)")
plt.ylabel("Electron density")
plt.title("Hydrogen — electron density")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

## Advanced tasks

### Custom mesh parameters

```python
from atom import AtomicDFTSolver

solver = AtomicDFTSolver(
    atomic_number=1,
    xc_functional="LDA_PW",
    domain_size=30.0,
    finite_element_number=20,
    polynomial_order=25,
    mesh_type="polynomial",
    mesh_concentration=2.0,
)
result = solver.solve()
print(f"Energy with custom mesh: {result['energy']:.6f} Ha")
```

### Check convergence

```python
from atom import AtomicDFTSolver

solver = AtomicDFTSolver(
    atomic_number=1,
    xc_functional="LDA_PW",
    scf_tolerance=1e-10,
)
result = solver.solve()
print(f"Converged: {result['converged']}")
print(f"Iterations: {result['iterations']}")
print(f"rho residual: {result['rho_residual']:.2e}")
```

## Tips

1. Start with hydrogen (`atomic_number=1`) to validate the environment.
2. Check `result["converged"]` before trusting energies for production runs.
3. Heavier atoms often need a larger `domain_size` and/or more `finite_element_number`.
