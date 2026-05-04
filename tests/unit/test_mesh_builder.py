"""
Unit tests for mesh and interpolation building blocks.

This file contains the following test cases:
1) test_quadrature1d
   - Verifies Gauss-Legendre quadrature by integrating 1 / (1 + x^2) on [-1, 1].
   - Compares numerical integral with the analytical reference pi/2.

2) test_lobatto1d
   - Verifies Lobatto quadrature with the same benchmark integrand/reference.

3) test_mesh1d_grid
   - Builds a 1D exponential mesh and checks basic structural properties:
     node/element counts and monotonic increase of mesh nodes.

4) test_mesh1d_fe_nodes
   - Generates FE nodes from mesh nodes + Lobatto interpolation nodes.
   - Checks output dimensionality and expected flattened-node count.

5) test_fe_flat_to_block2d
   - Validates reshaping flat FE data into element blocks for:
     (a) endpoints_shared=True and (b) endpoints_shared=False.
   - Verifies invalid input lengths raise AssertionError for both layout modes.

6) test_lagrange_shape_functions_lagrange_basis_and_derivatives
   - Validates Lagrange basis/derivative evaluation:
     shape checks, partition of unity, derivative sum consistency,
     nodal interpolation identity, and finite-value checks.

How to run:
- Run only this file:
  python -m pytest tests/unit/test_mesh_builder.py -q

- Run only tests in this module by keyword:
  python -m pytest -k mesh_builder -q

- Run all unit tests:
  python -m pytest tests/unit -q
"""

from pathlib import Path
import sys

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from atom.mesh.builder import LagrangeShapeFunctions, Mesh1D, Quadrature1D


def _integrand(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + x**2)


def test_quadrature1d():
    r"""Validate Gauss-Legendre quadrature with an analytical benchmark.

    Reference integral:
        \int_{-1}^{1} 1/(1+x^2) dx = pi/2.
    """
    nodes, legendre_weights = Quadrature1D.gauss_legendre(95)
    integral = np.sum(_integrand(nodes) * legendre_weights)
    assert np.isclose(integral, np.pi / 2, atol=1e-6)


def test_lobatto1d():
    r"""Validate Lobatto quadrature with the same analytical benchmark.

    Uses the same target value/tolerance as test_quadrature1d for consistency.
    """
    nodes, lobatto_weights = Quadrature1D.lobatto(31)
    integral = np.sum(_integrand(nodes) * lobatto_weights)
    assert np.isclose(integral, np.pi / 2, atol=1e-6)


def test_mesh1d_grid():
    """Validate basic structural properties of a generated 1D mesh.

    Checks:
    - node count is finite_elements_number + 1
    - element-width count is finite_elements_number
    - mesh nodes are strictly increasing
    """
    mesh = Mesh1D(
        domain_size            = 10.0,
        finite_elements_number = 17,
        mesh_type              = "exponential",
        clustering_param       = 61.0,
        exp_shift              = 0.0,
    )
    mesh_nodes, mesh_width = mesh.generate_mesh_nodes_and_width()
    assert len(mesh_nodes) == 18
    assert len(mesh_width) == 17
    assert np.all(np.diff(mesh_nodes) > 0)


def test_mesh1d_fe_nodes():
    """Validate FE-node generation from mesh nodes and interpolation nodes.

    For a shared-endpoint FE layout, the expected total number of flattened
    FE nodes is:
        finite_elements_number * interpolation_nodes_number + 1

    Checks:
    - output is a 1D flattened array
    - FE node count is larger than coarse mesh node count
    - FE node count matches the reference formula above
    """

    FINITE_ELEMENTS_NUMBER = 17
    INTERPOLATION_NODES_NUMBER = 31
    FE_NODES_REFERENCE = FINITE_ELEMENTS_NUMBER * INTERPOLATION_NODES_NUMBER + 1
    mesh = Mesh1D(
        domain_size            = 10.0,
        finite_elements_number = FINITE_ELEMENTS_NUMBER,
        mesh_type              = "exponential",
        clustering_param       = 61.0,
        exp_shift              = 0.0,
    )
    mesh_nodes, _ = mesh.generate_mesh_nodes_and_width()
    interp_nodes, _ = Quadrature1D.lobatto(INTERPOLATION_NODES_NUMBER)
    fe_nodes = Mesh1D.generate_fe_nodes(mesh_nodes, interp_nodes)
    assert fe_nodes.ndim == 1
    assert len(fe_nodes) > len(mesh_nodes)
    assert len(fe_nodes) == FE_NODES_REFERENCE
    


def test_fe_flat_to_block2d():
    """Validate FE flat-array to block conversion under both layout modes.

    Cases covered:
    - endpoints_shared=True: element blocks share boundary points
    - endpoints_shared=False: direct block reshape without shared endpoints
    - invalid input length: function should raise AssertionError
    """
    flat1 = np.arange(10, dtype=float)
    out1 = Mesh1D.fe_flat_to_block2d(flat1, 3, endpoints_shared=True)
    expected1 = np.array([[0, 1, 2, 3], [3, 4, 5, 6], [6, 7, 8, 9]], dtype=float)
    assert out1.shape == (3, 4)
    np.testing.assert_array_equal(out1, expected1)

    flat2 = np.arange(12, dtype=float)
    out2 = Mesh1D.fe_flat_to_block2d(flat2, 3, endpoints_shared=False)
    expected2 = np.array([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]], dtype=float)
    assert out2.shape == (3, 4)
    np.testing.assert_array_equal(out2, expected2)

    bad_flat_shared = np.arange(9)
    with pytest.raises(AssertionError):
        Mesh1D.fe_flat_to_block2d(bad_flat_shared, 3, endpoints_shared=True)
    bad_flat_stacked = np.arange(10)
    with pytest.raises(AssertionError):
        Mesh1D.fe_flat_to_block2d(bad_flat_stacked, 3, endpoints_shared=False)


def test_lagrange_shape_functions_lagrange_basis_and_derivatives():
    """Validate core mathematical properties of Lagrange basis functions.

    This test checks:
    - output tensor shapes for basis and derivatives
    - partition of unity: sum_k L_k(x) = 1
    - derivative consistency: sum_k dL_k/dx = 0
    - nodal interpolation identity: L_k(x_j) = delta_{kj}
    - finite-value behavior (no NaN/Inf)
    """
    # Define a simple 1-element nodal setup and evaluation grid.
    nodes = np.array([0.0, 1.0, 2.0, 3.0])[None, :]
    x_eval = np.linspace(0.0, 3.0, 31)[None, :]

    # Evaluate Lagrange basis values and derivatives on x_eval.
    L, dLdx = LagrangeShapeFunctions.lagrange_basis_and_derivatives(nodes, x_eval)

    # Verify output tensor shapes: (batch, n_eval, n_nodes).
    assert L.shape == (1, x_eval.shape[1], nodes.shape[1])
    assert dLdx.shape == (1, x_eval.shape[1], nodes.shape[1])

    # Partition-of-unity check: sum of basis values equals 1.
    unity_error = np.max(np.abs(np.sum(L[0, :, :], axis=1) - 1.0))
    assert unity_error < 1e-12

    # Derivative-consistency check: sum of derivatives equals 0.
    deriv_sum_error = np.max(np.abs(np.sum(dLdx[0, :, :], axis=1)))
    assert deriv_sum_error < 1e-10

    # Nodal-identity check against the identity matrix.
    L_nodes, _ = LagrangeShapeFunctions.lagrange_basis_and_derivatives(nodes, nodes)
    identity_error = np.max(np.abs(L_nodes[0, :, :] - np.eye(nodes.shape[1])))
    assert identity_error < 1e-12

    # Numerical sanity check: all values should be finite.
    assert np.all(np.isfinite(L))
    assert np.all(np.isfinite(dLdx))
