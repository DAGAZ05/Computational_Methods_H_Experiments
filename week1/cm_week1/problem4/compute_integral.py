import math
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np


def reference_integral(n: int, m: int = 200000) -> float:
    """Reference value from high-resolution trapezoidal integration."""
    x = np.linspace(0.0, 1.0, m + 1, dtype=np.float64)
    y = np.power(x, n) * np.exp(x - 1.0)
    return float(np.trapezoid(y, x))


def relative_error(computed: float, exact: float) -> float:
    """Relative error definition required by user: (x* - x) / x*."""
    if computed == 0.0:
        if exact == 0.0:
            return 0.0
        return math.copysign(math.inf, computed - exact)
    return (computed - exact) / computed


def forward_recurrence(n_max: int = 20, dtype=np.float32) -> np.ndarray:
    """y_n = 1 - n y_{n-1}, n=1..n_max, with y_0 = 1 - e^{-1}."""
    y = np.zeros(n_max + 1, dtype=dtype)
    y[0] = dtype(1.0 - math.exp(-1.0))
    for n in range(1, n_max + 1):
        y[n] = dtype(1.0 - n * y[n - 1])
    return y


def backward_recurrence(n_max: int = 20, y_nmax: float = 0.0, dtype=np.float32) -> np.ndarray:
    """y_{n-1} = (1 - y_n)/n, n=n_max..1, with given y_nmax."""
    y = np.zeros(n_max + 1, dtype=dtype)
    y[n_max] = dtype(y_nmax)
    for n in range(n_max, 0, -1):
        y[n - 1] = dtype((1.0 - y[n]) / n)
    return y


def build_results(n_max: int = 20, dtype=np.float32) -> Dict[str, np.ndarray]:
    exact = np.array([reference_integral(n) for n in range(n_max + 1)], dtype=np.float64)
    forward = forward_recurrence(n_max=n_max, dtype=dtype).astype(np.float64)
    backward0 = backward_recurrence(n_max=n_max, y_nmax=0.0, dtype=dtype).astype(np.float64)

    rel_f = np.array([relative_error(forward[n], exact[n]) for n in range(n_max + 1)], dtype=np.float64)
    rel_b = np.array([relative_error(backward0[n], exact[n]) for n in range(n_max + 1)], dtype=np.float64)

    return {
        "n": np.arange(n_max + 1, dtype=int),
        "exact": exact,
        "forward": forward,
        "backward0": backward0,
        "rel_forward": rel_f,
        "rel_backward0": rel_b,
    }


def print_key_table(results: Dict[str, np.ndarray]) -> None:
    print("Problem 4: y_n = integral_0^1 x^n * exp(x-1) dx")
    print("Relative error: (x* - x) / x*")
    print("Backward recurrence uses y_20 = 0")
    print()
    print("n   exact             forward           rel_err_fwd        backward(y20=0)    rel_err_bwd")
    print("-" * 96)

    for i in range(len(results["n"])):
        n = int(results["n"][i])
        ex = results["exact"][i]
        fw = results["forward"][i]
        bw = results["backward0"][i]
        erf = results["rel_forward"][i]
        erb = results["rel_backward0"][i]

        print(
            f"{n:>2d}  {ex: .8e}  {fw: .8e}  {erf: .8e}  {bw: .8e}  {erb: .8e}"
        )


def plot_results(results: Dict[str, np.ndarray], output_path: str = "visualization.png") -> None:
    n = results["n"]
    exact = results["exact"]
    forward = results["forward"]
    backward0 = results["backward0"]

    err_f = np.abs(results["rel_forward"])
    err_b = np.abs(results["rel_backward0"])
    floor = np.finfo(np.float64).tiny
    err_f = np.maximum(err_f, floor)
    err_b = np.maximum(err_b, floor)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    axes[0].plot(n, exact, "k-", label="reference", marker="o", markersize=3)
    axes[0].plot(n, forward, "r--", label="forward recurrence", marker="s", markersize=3)
    axes[0].plot(n, backward0, "b-.", label="backward recurrence (y20=0)", marker="^", markersize=3)
    axes[0].set_xlabel("n")
    axes[0].set_ylabel("y_n")
    axes[0].set_title("Computed y_n")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8)

    axes[1].semilogy(n, err_f, "r--", label="|rel_err_fwd|", marker="s", markersize=4)
    axes[1].semilogy(n, err_b, "b-.", label="|rel_err_bwd|", marker="^", markersize=4)
    axes[1].set_xlabel("n")
    axes[1].set_ylabel("abs((x*-x)/x*)")
    axes[1].set_title("Relative Error")
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    results = build_results(n_max=20, dtype=np.float32)
    print_key_table(results)
    plot_results(results, output_path="visualization.png")
    print("\nSaved visualization: problem4/visualization.png")


if __name__ == "__main__":
    main()
