import math
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np


def relative_error(computed: float, exact: float) -> float:
    """Relative error defined as (x* - x) / x*, where x* is computed value."""
    if computed == 0.0:
        if exact == 0.0:
            return 0.0
        return math.copysign(math.inf, computed - exact)
    return (computed - exact) / computed


def exp_taylor_direct(x: float, tol: float = 1e-7, max_terms: int = 500, dtype=np.float32) -> Tuple[float, int]:
    """Compute e^x via direct Taylor series in the given floating precision."""
    x_val = dtype(x)
    term = dtype(1.0)
    s = dtype(1.0)
    n = 0

    while n < max_terms - 1:
        n += 1
        term = dtype(term * x_val / dtype(n))
        s = dtype(s + term)
        if abs(float(term)) < tol:
            break

    return float(s), n + 1


def exp_taylor_reciprocal_for_negative(x: float, tol: float = 1e-7, max_terms: int = 500, dtype=np.float32) -> Tuple[float, int]:
    """For x<0, compute e^x as 1 / e^{-x}."""
    if x >= 0.0:
        return exp_taylor_direct(x, tol=tol, max_terms=max_terms, dtype=dtype)

    pos_val, terms = exp_taylor_direct(-x, tol=tol, max_terms=max_terms, dtype=dtype)
    if pos_val == 0.0:
        return math.inf, terms
    return 1.0 / pos_val, terms


def evaluate_points(xs: List[float], method: str, tol: float, dtype=np.float32) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []

    for x in xs:
        if method == "direct":
            approx, terms = exp_taylor_direct(x, tol=tol, dtype=dtype)
        elif method == "reciprocal":
            approx, terms = exp_taylor_reciprocal_for_negative(x, tol=tol, dtype=dtype)
        else:
            raise ValueError("method must be 'direct' or 'reciprocal'")

        exact = math.exp(x)
        rel = relative_error(approx, exact)

        rows.append(
            {
                "x": x,
                "computed": approx,
                "exact": exact,
                "rel_error": rel,
                "terms": terms,
            }
        )

    return rows


def print_table(title: str, rows: List[Dict[str, float]]) -> None:
    print(title)
    print("-" * 88)
    print("x         computed(x*)       exact(x)          rel_err=(x*-x)/x*      terms")
    print("-" * 88)
    for r in rows:
        print(
            f"{r['x']:>6.1f}    {r['computed']:>13.20e}   {r['exact']:>13.20e}   "
            f"{r['rel_error']:>18.20e}   {int(r['terms']):>5d}"
        )
    print()


def visualize(rows_pos_direct, rows_neg_direct, rows_neg_reciprocal, output_path="visualization.png") -> None:
    x_pos = [r["x"] for r in rows_pos_direct]
    x_neg_abs = [abs(r["x"]) for r in rows_neg_direct]

    err_pos = [abs(r["rel_error"]) for r in rows_pos_direct]
    err_neg_direct = [abs(r["rel_error"]) for r in rows_neg_direct]
    err_neg_recip = [abs(r["rel_error"]) for r in rows_neg_reciprocal]

    tiny = np.finfo(np.float64).tiny
    err_pos = np.maximum(err_pos, tiny)
    err_neg_direct = np.maximum(err_neg_direct, tiny)
    err_neg_recip = np.maximum(err_neg_recip, tiny)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].semilogy(x_pos, err_pos, "o-", label="direct Taylor")
    axes[0].set_title("x > 0")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("abs((x*-x)/x*)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].semilogy(x_neg_abs, err_neg_direct, "o-", label="direct Taylor")
    axes[1].semilogy(x_neg_abs, err_neg_recip, "s-", label="reciprocal formula")
    axes[1].set_title("x < 0")
    axes[1].set_xlabel("|x|")
    axes[1].set_ylabel("abs((x*-x)/x*)")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    tol = 1e-7
    dtype = np.float32

    x_pos = [1.0, 5.0, 10.0, 15.0, 20.0]
    x_neg = [-1.0, -5.0, -10.0, -15.0, -20.0]

    rows_pos_direct = evaluate_points(x_pos, method="direct", tol=tol, dtype=dtype)
    rows_neg_direct = evaluate_points(x_neg, method="direct", tol=tol, dtype=dtype)
    rows_neg_recip = evaluate_points(x_neg, method="reciprocal", tol=tol, dtype=dtype)

    print("Problem 3: exp(x) by Taylor series (float32)")
    print("Relative error definition: (x* - x) / x*")
    print(f"Stopping criterion: |term| < {tol:.0e}")
    print()

    print_table("Table A: x = 1, 5, 10, 15, 20 (direct Taylor)", rows_pos_direct)
    print_table("Table B: x = -1, -5, -10, -15, -20 (direct Taylor)", rows_neg_direct)
    print_table("Table C: x = -1, -5, -10, -15, -20 (reciprocal formula)", rows_neg_recip)

    visualize(rows_pos_direct, rows_neg_direct, rows_neg_recip, output_path="visualization.png")
    print("Saved visualization: problem3/visualization.png")


if __name__ == "__main__":
    main()
