#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational Methods Week 2: Solve equation 4cos(x) = e^x
Implement four iterative methods and compare their efficiency
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import time
from typing import Tuple, List, Dict


class EquationSolver:
    """Solver for equation 4cos(x) = e^x"""

    def __init__(self, epsilon: float = 1e-4, verbose: bool = False):
        """
        Initialize solver

        Args:
            epsilon: Tolerance for convergence
            verbose: Print iteration details
        """
        self.epsilon = epsilon
        self.max_iter = 1000
        self.verbose = verbose

    def f(self, x: float) -> float:
        """Original equation: f(x) = e^x - 4cos(x) = 0"""
        return np.exp(x) - 4 * np.cos(x)

    def df(self, x: float) -> float:
        """Derivative of f(x)"""
        return np.exp(x) + 4 * np.sin(x)

    def phi(self, x: float) -> float:
        """Iteration function for simple iteration: x = phi(x)"""
        # From e^x = 4cos(x), we get x = arccos(e^x / 4)
        return np.arccos(np.exp(x) / 4)

    def simple_iteration(self, x0: float) -> Tuple[float, int, List[float]]:
        """
        Simple iteration method

        Args:
            x0: Initial value

        Returns:
            (root, iterations, history)
        """
        x = x0
        history = [x0]

        if self.verbose:
            print(f"{'k':<5} {'x_k':<15} {'|x_k - x_(k-1)|':<20}")
            print("-" * 40)
            print(f"{0:<5} {x0:<15.10f} {'-':<20}")

        for i in range(self.max_iter):
            x_new = self.phi(x)
            history.append(x_new)
            diff = abs(x_new - x)

            if self.verbose:
                print(f"{i+1:<5} {x_new:<15.10f} {diff:<20.10e}")

            if diff < self.epsilon:
                return x_new, i + 1, history
            x = x_new
        raise ValueError("Simple iteration method did not converge")

    def steffensen(self, x0: float) -> Tuple[float, int, List[float]]:
        """
        Steffensen's iteration method (accelerated convergence)

        Args:
            x0: Initial value

        Returns:
            (root, iterations, history)
        """
        x = x0
        history = [x0]

        if self.verbose:
            print(f"{'k':<5} {'x_k':<15} {'|x_k - x_(k-1)|':<20}")
            print("-" * 40)
            print(f"{0:<5} {x0:<15.10f} {'-':<20}")

        for i in range(self.max_iter):
            phi_x = self.phi(x)
            phi_phi_x = self.phi(phi_x)
            denominator = phi_phi_x - 2 * phi_x + x
            if abs(denominator) < 1e-12:
                raise ValueError("Steffensen method: denominator near zero")
            x_new = x - (phi_x - x) ** 2 / denominator
            history.append(x_new)
            diff = abs(x_new - x)

            if self.verbose:
                print(f"{i+1:<5} {x_new:<15.10f} {diff:<20.10e}")

            if diff < self.epsilon:
                return x_new, i + 1, history
            x = x_new
        raise ValueError("Steffensen method did not converge")

    def newton(self, x0: float) -> Tuple[float, int, List[float]]:
        """
        Newton's iteration method

        Args:
            x0: Initial value

        Returns:
            (root, iterations, history)
        """
        x = x0
        history = [x0]

        if self.verbose:
            print(f"{'k':<5} {'x_k':<15} {'|x_k - x_(k-1)|':<20}")
            print("-" * 40)
            print(f"{0:<5} {x0:<15.10f} {'-':<20}")

        for i in range(self.max_iter):
            fx = self.f(x)
            dfx = self.df(x)
            if abs(dfx) < 1e-12:
                raise ValueError("Newton method: derivative near zero")
            x_new = x - fx / dfx
            history.append(x_new)
            diff = abs(x_new - x)

            if self.verbose:
                print(f"{i+1:<5} {x_new:<15.10f} {diff:<20.10e}")

            if diff < self.epsilon:
                return x_new, i + 1, history
            x = x_new
        raise ValueError("Newton method did not converge")

    def secant(self, x0: float, x1: float) -> Tuple[float, int, List[float]]:
        """
        Secant method (two-point chord method)

        Args:
            x0: First initial value
            x1: Second initial value

        Returns:
            (root, iterations, history)
        """
        history = [x0, x1]

        if self.verbose:
            print(f"{'k':<5} {'x_k':<15} {'|x_k - x_(k-1)|':<20}")
            print("-" * 40)
            print(f"{0:<5} {x0:<15.10f} {'-':<20}")
            print(f"{1:<5} {x1:<15.10f} {abs(x1-x0):<20.10e}")

        for i in range(self.max_iter):
            f0 = self.f(x0)
            f1 = self.f(x1)
            if abs(f1 - f0) < 1e-12:
                raise ValueError("Secant method: denominator near zero")
            x_new = x1 - f1 * (x1 - x0) / (f1 - f0)
            history.append(x_new)
            diff = abs(x_new - x1)

            if self.verbose:
                print(f"{i+2:<5} {x_new:<15.10f} {diff:<20.10e}")

            if diff < self.epsilon:
                return x_new, i + 1, history
            x0, x1 = x1, x_new
        raise ValueError("Secant method did not converge")


def run_method(solver: EquationSolver, method: str) -> Dict:
    """Run specified method and return results"""
    x0 = np.pi / 4
    x1 = np.pi / 2

    start_time = time.perf_counter()

    try:
        if method == "simple":
            root, iterations, history = solver.simple_iteration(x0)
        elif method == "steffensen":
            root, iterations, history = solver.steffensen(x0)
        elif method == "newton":
            root, iterations, history = solver.newton(x0)
        elif method == "secant":
            root, iterations, history = solver.secant(x0, x1)
        else:
            raise ValueError(f"Unknown method: {method}")

        elapsed_time = time.perf_counter() - start_time

        return {
            "method": method,
            "root": root,
            "iterations": iterations,
            "time": elapsed_time,
            "history": history,
            "error": abs(solver.f(root)),
        }
    except Exception as e:
        return {"method": method if "method_name" in locals() else method, "error": str(e)}


def plot_results(results: List[Dict], solver: EquationSolver):
    """Plot scientific-style visualization"""
    # Set scientific publication style
    plt.style.use("seaborn-v0_8-paper")

    # Configure font and fix minus sign issue
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial"],
            "axes.unicode_minus": False,  # Fix minus sign display
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.figsize": (12, 8),
            "figure.dpi": 100,
        }
    )

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # 1. Function plot and root locations
    ax1 = fig.add_subplot(gs[0, :])
    x_range = np.linspace(-1, 2, 500)
    y1 = np.exp(x_range)
    y2 = 4 * np.cos(x_range)
    ax1.plot(x_range, y1, "b-", linewidth=1.5, label=r"$e^x$")
    ax1.plot(x_range, y2, "r-", linewidth=1.5, label=r"$4\cos(x)$")
    ax1.axhline(y=0, color="k", linestyle="--", linewidth=0.5, alpha=0.3)
    ax1.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)
    for result in results:
        if "root" in result:
            ax1.plot(result["root"], np.exp(result["root"]), "o", markersize=8, label=f'{result["method"]}')
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_title(r"Equation: $4\cos(x) = e^x$")
    ax1.legend(loc="best", frameon=True, shadow=False)

    # 2. Iteration count comparison
    ax2 = fig.add_subplot(gs[1, 0])
    methods = [r["method"] for r in results if "iterations" in r]
    iterations = [r["iterations"] for r in results if "iterations" in r]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    bars = ax2.bar(range(len(methods)), iterations, color=colors[: len(methods)], alpha=0.8, edgecolor="black", linewidth=0.8)
    ax2.set_xticks(range(len(methods)))
    ax2.set_xticklabels(methods, rotation=15, ha="right")
    ax2.set_ylabel("Iterations")
    ax2.set_title("Iteration Count Comparison")
    ax2.grid(True, axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
    for i, (bar, val) in enumerate(zip(bars, iterations)):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5, str(val), ha="center", va="bottom", fontsize=9)

    # 3. Computation time comparison
    ax3 = fig.add_subplot(gs[1, 1])
    times = [r["time"] * 1000 for r in results if "time" in r]
    bars = ax3.bar(range(len(methods)), times, color=colors[: len(methods)], alpha=0.8, edgecolor="black", linewidth=0.8)
    ax3.set_xticks(range(len(methods)))
    ax3.set_xticklabels(methods, rotation=15, ha="right")
    ax3.set_ylabel("Time (ms)")
    ax3.set_title("Computation Time Comparison")
    ax3.grid(True, axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
    for i, (bar, val) in enumerate(zip(bars, times)):
        ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f"{val:.4f}", ha="center", va="bottom", fontsize=8)

    # 4. Convergence process (error vs iteration)
    ax4 = fig.add_subplot(gs[2, 0])
    for i, result in enumerate(results):
        if "history" in result:
            errors = [abs(solver.f(x)) for x in result["history"]]
            ax4.semilogy(range(len(errors)), errors, "o-", linewidth=1.5, markersize=4, label=result["method"], color=colors[i])
    ax4.set_xlabel("Iteration")
    ax4.set_ylabel(r"$|f(x)|$ (log scale)")
    ax4.set_title("Convergence Process")
    ax4.legend(loc="best", frameon=True, shadow=False)
    ax4.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)

    # 5. Iteration point evolution
    ax5 = fig.add_subplot(gs[2, 1])
    for i, result in enumerate(results):
        if "history" in result:
            ax5.plot(range(len(result["history"])), result["history"], "o-", linewidth=1.5, markersize=4, label=result["method"], color=colors[i])
    ax5.set_xlabel("Iteration")
    ax5.set_ylabel("x")
    ax5.set_title("Iteration Point Evolution")
    ax5.legend(loc="best", frameon=True, shadow=False)
    ax5.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.savefig("results.png", dpi=300, bbox_inches="tight")
    print("\nVisualization saved to results.png")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Solve equation 4cos(x) = e^x")
    parser.add_argument(
        "method",
        choices=["simple", "steffensen", "newton", "secant", "all"],
        help="Choose method: simple, steffensen, newton, secant, or all",
    )
    parser.add_argument("--epsilon", type=float, default=1e-4, help="Tolerance (default: 1e-4)")
    parser.add_argument("--verbose", action="store_true", help="Print iteration details")

    args = parser.parse_args()

    solver = EquationSolver(epsilon=args.epsilon, verbose=args.verbose)

    print(f"\n{'=' * 60}")
    print(f"Solving equation: 4cos(x) = e^x")
    print(f"Tolerance: epsilon = {args.epsilon}")
    print(f"{'=' * 60}\n")

    if args.method == "all":
        methods = ["simple", "steffensen", "newton", "secant"]
        results = []
        for method in methods:
            print(f"\n--- {method.upper()} METHOD ---")
            result = run_method(solver, method)
            results.append(result)
            if "root" in result:
                print(f"\n{result['method']}:")
                print(f"  Root: x = {result['root']:.10f}")
                print(f"  Iterations: {result['iterations']}")
                print(f"  Time: {result['time']*1000:.4f} ms")
                print(f"  |f(x)|: {result['error']:.2e}")
            else:
                print(f"{result['method']}: Failed - {result['error']}")

        print(f"\n{'=' * 60}")
        print("EFFICIENCY COMPARISON:")
        print(f"{'=' * 60}")
        valid_results = [r for r in results if "iterations" in r]
        if valid_results:
            min_iter = min(r["iterations"] for r in valid_results)
            min_time = min(r["time"] for r in valid_results)
            for r in valid_results:
                iter_ratio = r["iterations"] / min_iter
                time_ratio = r["time"] / min_time
                print(f"{r['method']:20s}: Iter ratio = {iter_ratio:.2f}x, Time ratio = {time_ratio:.2f}x")

        plot_results(results, solver)
    else:
        print(f"\n--- {args.method.upper()} METHOD ---")
        result = run_method(solver, args.method)
        if "root" in result:
            print(f"\n{result['method']}:")
            print(f"  Root: x = {result['root']:.10f}")
            print(f"  Iterations: {result['iterations']}")
            print(f"  Time: {result['time']*1000:.4f} ms")
            print(f"  |f(x)|: {result['error']:.2e}")
        else:
            print(f"{result['method']}: Failed - {result['error']}")


if __name__ == "__main__":
    main()

