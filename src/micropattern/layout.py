from typing import Literal

import numpy as np
from anndata import AnnData
from fa2 import ForceAtlas2

__all__ = ["force_layout"]


def _random_disk(r: float, n_points: int, noise: float, seed: int) -> np.ndarray:
    np.random.seed(seed)
    rho = r * np.sqrt(np.random.rand(n_points))
    theta = 2 * np.pi * np.random.rand(n_points)
    x = rho * np.cos(theta) + np.random.normal(0, noise, n_points)
    y = rho * np.sin(theta) + np.random.normal(0, noise, n_points)
    return np.column_stack([x, y])


def _random_ball(r: float, n_points: int, dim: int, noise: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    directions = rng.standard_normal((n_points, dim))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    radii = r * rng.uniform(0, 1, size=n_points) ** (1.0 / dim)
    pos = directions * radii[:, None]
    if noise > 0:
        pos += rng.normal(0, noise, size=pos.shape)
    return pos


def force_layout(
    adata: AnnData,
    use_rep: str,
    *,
    dim: int = 2,
    initial_position: Literal["random_disk", "random_ball"] | np.ndarray | None = None,
    distance_key: str = "connectivities",
    num_iterations: int = 10,
    scalingRatio: float = 2.0,
    gravity: float = 1.0,
    seed: int = 0,
    key_added: str = "X_fl",
) -> None:
    assert distance_key in adata.obsp, f"distance key '{distance_key}' not found in adata.obsp"
    knn_dist_mat = adata.obsp[distance_key]
    data = adata.X if use_rep == "X" else adata.obsm[use_rep]
    n_points = data.shape[0]

    if isinstance(initial_position, np.ndarray):
        pos = initial_position
    elif initial_position == "random_disk":
        assert dim == 2, "random_disk only supports dim=2, use random_ball for higher dimensions"
        pos = _random_disk(r=1, n_points=n_points, noise=0, seed=seed)
    elif initial_position == "random_ball":
        pos = _random_ball(r=1, n_points=n_points, dim=dim, noise=0, seed=seed)
    else:
        pos = None

    fa = ForceAtlas2(verbose=True, scalingRatio=scalingRatio, gravity=gravity, dim=dim)
    adata.obsm[key_added] = np.array(
        fa.forceatlas2(G=knn_dist_mat, pos=pos, iterations=num_iterations)
    )
