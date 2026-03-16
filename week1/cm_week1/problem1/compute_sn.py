import argparse
import math
import time
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np

try:
    from numba import njit

    NUMBA_AVAILABLE = True
except Exception:
    NUMBA_AVAILABLE = False
    njit = None


def exact_sn(n: int) -> float:
    """Closed-form value in float64 (reference only)."""
    return 0.5 * (1.5 - 1.0 / n - 1.0 / (n + 1))



def relative_error_signed(x_star: float, x_exact: float) -> float:
    """Relative error defined as (x* - x) / x*."""
    if x_star == 0.0:
        if x_exact == 0.0:
            return 0.0
        return math.copysign(math.inf, x_star - x_exact)
    return (x_star - x_exact) / x_star
def count_significant_figures(approx: np.float32, exact: float) -> int:
    """Count significant digits based on |x*-x| <= 0.5 * 10^(m-n)."""
    if exact == 0.0:
        return 0

    abs_error = abs(float(approx) - exact)
    if abs_error == 0.0:
        return 15

    m = int(math.floor(math.log10(abs(exact)))) + 1
    n_digits = m - math.log10(2.0 * abs_error)
    return max(0, int(math.floor(n_digits)))


def _accumulate_terms_in_order(s: np.float32, terms: np.ndarray) -> np.float32:
    # Keep strict left-to-right float32 accumulation in C (np.cumsum).
    tmp = np.empty(terms.size + 1, dtype=np.float32)
    tmp[0] = s
    tmp[1:] = terms
    return np.cumsum(tmp, dtype=np.float32)[-1]


def _compute_forward_np(n: int) -> np.float32:
    s = np.float32(0.0)
    one = np.float32(1.0)
    chunk_size = 5_000_000

    start = 2
    while start <= n:
        end = min(n + 1, start + chunk_size)
        j = np.arange(start, end, dtype=np.int64)
        jf = j.astype(np.float32)
        terms = one / (jf * jf - one)
        s = _accumulate_terms_in_order(s, terms)
        start = end

    return np.float32(s)


def _compute_backward_np(n: int) -> np.float32:
    s = np.float32(0.0)
    one = np.float32(1.0)
    chunk_size = 5_000_000

    end = n
    while end >= 2:
        start = max(2, end - chunk_size + 1)
        j = np.arange(end, start - 1, -1, dtype=np.int64)
        jf = j.astype(np.float32)
        terms = one / (jf * jf - one)
        s = _accumulate_terms_in_order(s, terms)
        end = start - 1

    return np.float32(s)


if NUMBA_AVAILABLE:

    @njit(cache=True)
    def _compute_forward_numba(n: int) -> np.float32:
        s = np.float32(0.0)
        one = np.float32(1.0)
        for j in range(2, n + 1):
            jf = np.float32(j)
            s = np.float32(s + one / (np.float32(jf * jf) - one))
        return s

    @njit(cache=True)
    def _compute_backward_numba(n: int) -> np.float32:
        s = np.float32(0.0)
        one = np.float32(1.0)
        for j in range(n, 1, -1):
            jf = np.float32(j)
            s = np.float32(s + one / (np.float32(jf * jf) - one))
        return s


    def compute_forward(n: int) -> np.float32:
        return _compute_forward_numba(n)


    def compute_backward(n: int) -> np.float32:
        return _compute_backward_numba(n)

else:

    def compute_forward(n: int) -> np.float32:
        return _compute_forward_np(n)


    def compute_backward(n: int) -> np.float32:
        return _compute_backward_np(n)


def _warmup_jit() -> None:
    if NUMBA_AVAILABLE:
        _ = _compute_forward_numba(10)
        _ = _compute_backward_numba(10)


def _print_backend() -> None:
    print(f"Backend: {'numba-jit' if NUMBA_AVAILABLE else 'numpy-chunk-cumsum'}")


def run_forward_only(k_min: int = 1, k_max: int = 10) -> List[Dict[str, float]]:
    _warmup_jit()
    print("Problem 1 (1): forward order")
    _print_backend()

    results: List[Dict[str, float]] = []
    for k in range(k_min, k_max + 1):
        n = 10**k
        t0 = time.perf_counter()

        exact = exact_sn(n)
        forward = compute_forward(n)
        elapsed = time.perf_counter() - t0

        result = {
            "k": k,
            "n": n,
            "exact": exact,
            "forward": float(forward),
            "forward_digits": count_significant_figures(forward, exact),
            "time_s": elapsed,
        }
        results.append(result)

        print(
            f"k={k:2d} n={n:>11d} time={elapsed:8.3f}s "
            f"exact={exact:.9e} forward={result['forward']:.9e} "
            f"digits={result['forward_digits']}"
        )

    return results


def run_backward_only(k_min: int = 1, k_max: int = 10) -> List[Dict[str, float]]:
    _warmup_jit()
    print("Problem 1 (2): backward order")
    _print_backend()

    results: List[Dict[str, float]] = []
    for k in range(k_min, k_max + 1):
        n = 10**k
        t0 = time.perf_counter()

        exact = exact_sn(n)
        backward = compute_backward(n)
        elapsed = time.perf_counter() - t0

        result = {
            "k": k,
            "n": n,
            "exact": exact,
            "backward": float(backward),
            "backward_digits": count_significant_figures(backward, exact),
            "time_s": elapsed,
        }
        results.append(result)

        print(
            f"k={k:2d} n={n:>11d} time={elapsed:8.3f}s "
            f"exact={exact:.9e} backward={result['backward']:.9e} "
            f"digits={result['backward_digits']}"
        )

    return results


def run_benchmark_both(k_min: int = 1, k_max: int = 10) -> List[Dict[str, float]]:
    _warmup_jit()

    print("Compute S_n in float32 for n = 10^k")
    _print_backend()

    results: List[Dict[str, float]] = []
    for k in range(k_min, k_max + 1):
        n = 10**k
        t0 = time.perf_counter()

        exact = exact_sn(n)
        forward = compute_forward(n)
        backward = compute_backward(n)

        elapsed = time.perf_counter() - t0

        result = {
            "k": k,
            "n": n,
            "exact": exact,
            "forward": float(forward),
            "backward": float(backward),
            "forward_digits": count_significant_figures(forward, exact),
            "backward_digits": count_significant_figures(backward, exact),
            "time_s": elapsed,
        }
        results.append(result)

        print(
            f"k={k:2d} n={n:>11d} time={elapsed:8.3f}s "
            f"forward={result['forward']:.9e} ({result['forward_digits']:2d} digits) "
            f"backward={result['backward']:.9e} ({result['backward_digits']:2d} digits)"
        )

    return results


def visualize_results(results: List[Dict[str, float]], output_path: str = "visualization.png") -> None:
    k_values = [r["k"] for r in results]
    forward_errors = [abs(relative_error_signed(r["forward"], r["exact"])) for r in results]
    backward_errors = [abs(relative_error_signed(r["backward"], r["exact"])) for r in results]
    forward_digits = [r["forward_digits"] for r in results]
    backward_digits = [r["backward_digits"] for r in results]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].semilogy(k_values, forward_errors, "o-", label="forward")
    axes[0].semilogy(k_values, backward_errors, "s-", label="backward")
    axes[0].set_xlabel("k (n=10^k)")
    axes[0].set_ylabel("abs((x* - x)/x*)")
    axes[0].set_title("Relative Error")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(k_values, forward_digits, "o-", label="forward")
    axes[1].plot(k_values, backward_digits, "s-", label="backward")
    axes[1].set_xlabel("k (n=10^k)")
    axes[1].set_ylabel("significant digits")
    axes[1].set_title("Significant Digits")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved visualization: {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute S_n with float32 in different summation orders.")
    parser.add_argument(
        "--mode",
        choices=["forward", "backward", "both"],
        default="both",
        help="forward/backward/both",
    )
    parser.add_argument("--k-min", type=int, default=1, help="minimum k, n=10^k")
    parser.add_argument("--k-max", type=int, default=10, help="maximum k, n=10^k")
    parser.add_argument("--plot", default="problem1/visualization.png", help="plot output path (only for mode=both)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.k_min < 1 or args.k_max < args.k_min:
        raise ValueError("Require 1 <= k_min <= k_max")

    if args.mode == "forward":
        run_forward_only(args.k_min, args.k_max)
    elif args.mode == "backward":
        run_backward_only(args.k_min, args.k_max)
    else:
        results = run_benchmark_both(args.k_min, args.k_max)
        visualize_results(results, args.plot)


if __name__ == "__main__":
    main()

