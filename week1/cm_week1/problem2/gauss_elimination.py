import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def gauss_elimination(A, b, dtype=np.float32, pivoting=False, singular_tol=1e-30):
    """Solve Ax=b by Gaussian elimination."""
    A = np.array(A, dtype=dtype, copy=True)
    b = np.array(b, dtype=dtype, copy=True)
    n = b.size

    for i in range(n):
        if pivoting:
            max_row = i + np.argmax(np.abs(A[i:, i]))
            if max_row != i:
                A[[i, max_row]] = A[[max_row, i]]
                b[[i, max_row]] = b[[max_row, i]]

        pivot = A[i, i]
        if np.abs(pivot) <= singular_tol:
            raise ValueError(f"Zero or near-zero pivot at row {i}: {pivot}")

        for j in range(i + 1, n):
            factor = A[j, i] / pivot
            A[j, i:] = A[j, i:] - factor * A[i, i:]
            b[j] = b[j] - factor * b[i]

    x = np.zeros(n, dtype=dtype)
    for i in range(n - 1, -1, -1):
        if np.abs(A[i, i]) <= singular_tol:
            raise ValueError(f"Zero or near-zero pivot during back substitution at row {i}")
        x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]

    return x


def system1(epsilon, dtype=np.float32):
    # [eps, 1; 1, 1] [x1, x2]^T = [1+eps, 2]^T
    A = np.array([[epsilon, 1.0], [1.0, 1.0]], dtype=dtype)
    b = np.array([1.0 + epsilon, 2.0], dtype=dtype)
    return A, b


def system2(epsilon, dtype=np.float32):
    # [1, 1; eps, 1] [x2, x1]^T = [2, 1+eps]^T
    A = np.array([[1.0, 1.0], [epsilon, 1.0]], dtype=dtype)
    b = np.array([2.0, 1.0 + epsilon], dtype=dtype)
    return A, b


def component_relative_error_signed(x_num, x_exact):
    """Compute signed componentwise relative error: (x* - x)/x*.

    Here x* is the computed (machine) value x_num, and x is the exact value.
    """
    x_num = x_num.astype(np.float64)
    x_exact = x_exact.astype(np.float64)
    rel = np.full_like(x_num, np.nan, dtype=np.float64)
    nonzero = x_num != 0.0
    rel[nonzero] = (x_num[nonzero] - x_exact[nonzero]) / x_num[nonzero]

    # If x*=0, the ratio is unbounded unless x is also 0.
    zero_denom = ~nonzero
    numer = x_num - x_exact
    rel[zero_denom] = np.sign(numer[zero_denom]) * np.inf
    rel[zero_denom & (numer == 0.0)] = 0.0
    return rel

def effective_digits_from_rel_error(rel_error_value, dtype=np.float32):
    """If |e_r| <= 0.5*10^{-n}, n significant digits are guaranteed."""
    if np.isnan(rel_error_value):
        return np.nan
    if np.isinf(rel_error_value):
        return 0.0
    rel_abs = np.abs(rel_error_value)
    if rel_abs <= 0.0:
        return float(np.finfo(dtype).precision)
    n = np.floor(-np.log10(2.0 * rel_abs))
    return float(max(0.0, n))

def evaluate_components(x_num, x_exact, dtype=np.float32):
    rel = component_relative_error_signed(x_num, x_exact)
    d1 = effective_digits_from_rel_error(rel[0], dtype=dtype)
    d2 = effective_digits_from_rel_error(rel[1], dtype=dtype)
    return rel, np.array([d1, d2], dtype=np.float64)


def print_table(title, rows):
    print(title)
    print('-' * 114)
    print(
        'k   eps          x1               x2               '
        'rel_err_x1        rel_err_x2        digits_x1  digits_x2'
    )
    print('-' * 114)

    for row in rows:
        if row['status'] != 'ok':
            print(f"{row['k']:<3d} {row['eps']:>10.1e}   failed: {row['message']}")
            continue

        print(
            f"{row['k']:<3d} {row['eps']:>10.1e}   "
            f"{row['x'][0]: .8e}   {row['x'][1]: .8e}   "
            f"{row['rel'][0]: .8e}   {row['rel'][1]: .8e}   "
            f"{row['digits'][0]:>5.1f}      {row['digits'][1]:>5.1f}"
        )

    print()


def plot_results(rows_sys1, rows_sys2, output_path='visualization.png'):
    ok1 = [r for r in rows_sys1 if r['status'] == 'ok']
    ok2 = [r for r in rows_sys2 if r['status'] == 'ok']

    k1 = np.array([r['k'] for r in ok1], dtype=int)
    k2 = np.array([r['k'] for r in ok2], dtype=int)

    rel1_x1 = np.abs(np.array([r['rel'][0] for r in ok1], dtype=np.float64))
    rel1_x2 = np.abs(np.array([r['rel'][1] for r in ok1], dtype=np.float64))
    rel2_x1 = np.abs(np.array([r['rel'][0] for r in ok2], dtype=np.float64))
    rel2_x2 = np.abs(np.array([r['rel'][1] for r in ok2], dtype=np.float64))

    dig1_x1 = np.array([r['digits'][0] for r in ok1], dtype=np.float64)
    dig1_x2 = np.array([r['digits'][1] for r in ok1], dtype=np.float64)
    dig2_x1 = np.array([r['digits'][0] for r in ok2], dtype=np.float64)
    dig2_x2 = np.array([r['digits'][1] for r in ok2], dtype=np.float64)

    floor = np.finfo(np.float64).tiny
    rel1_x1 = np.maximum(rel1_x1, floor)
    rel1_x2 = np.maximum(rel1_x2, floor)
    rel2_x1 = np.maximum(rel2_x1, floor)
    rel2_x2 = np.maximum(rel2_x2, floor)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex='col')

    axes[0, 0].semilogy(k1, rel1_x1, 'o-', label='|rel_err_x1|')
    axes[0, 0].semilogy(k1, rel1_x2, 's-', label='|rel_err_x2|')
    axes[0, 0].set_title('Problem 2(1): Relative Error')
    axes[0, 0].set_ylabel('abs((x*-x)/|x*|)')
    axes[0, 0].grid(True, which='both', ls='--', alpha=0.4)
    axes[0, 0].legend()

    axes[1, 0].plot(k1, dig1_x1, 'o-', label='digits_x1')
    axes[1, 0].plot(k1, dig1_x2, 's-', label='digits_x2')
    axes[1, 0].set_xlabel('k (epsilon = 10^{-2k})')
    axes[1, 0].set_ylabel('significant digits')
    axes[1, 0].set_xticks(k1)
    axes[1, 0].grid(True, ls='--', alpha=0.4)
    axes[1, 0].legend()

    axes[0, 1].semilogy(k2, rel2_x1, 'o-', label='|rel_err_x1|')
    axes[0, 1].semilogy(k2, rel2_x2, 's-', label='|rel_err_x2|')
    axes[0, 1].set_title('Problem 2(2): Relative Error')
    axes[0, 1].set_ylabel('abs((x*-x)/|x*|)')
    axes[0, 1].grid(True, which='both', ls='--', alpha=0.4)
    axes[0, 1].legend()

    axes[1, 1].plot(k2, dig2_x1, 'o-', label='digits_x1')
    axes[1, 1].plot(k2, dig2_x2, 's-', label='digits_x2')
    axes[1, 1].set_xlabel('k (epsilon = 10^{-2k})')
    axes[1, 1].set_ylabel('significant digits')
    axes[1, 1].set_xticks(k2)
    axes[1, 1].grid(True, ls='--', alpha=0.4)
    axes[1, 1].legend()

    fig.suptitle('Problem 2: Componentwise Relative Error and Significant Digits')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def run_experiment(dtype=np.float32, pivoting=False):
    x_exact = np.array([1.0, 1.0], dtype=np.float64)

    rows_sys1 = []
    rows_sys2 = []

    for k in range(0, 11):
        eps = 10.0 ** (-2 * k)

        A1, b1 = system1(eps, dtype=dtype)
        A2, b2 = system2(eps, dtype=dtype)

        row1 = {'k': k, 'eps': eps, 'status': 'ok'}
        row2 = {'k': k, 'eps': eps, 'status': 'ok'}

        try:
            x1 = gauss_elimination(A1, b1, dtype=dtype, pivoting=pivoting)
            rel1, dig1 = evaluate_components(x1, x_exact, dtype=dtype)
            row1.update({'x': x1.astype(np.float64), 'rel': rel1, 'digits': dig1})
        except ValueError as exc:
            row1.update({'status': 'failed', 'message': str(exc)})

        try:
            x2_raw = gauss_elimination(A2, b2, dtype=dtype, pivoting=pivoting)
            x2 = np.array([x2_raw[1], x2_raw[0]], dtype=dtype)
            rel2, dig2 = evaluate_components(x2, x_exact, dtype=dtype)
            row2.update({'x': x2.astype(np.float64), 'rel': rel2, 'digits': dig2})
        except ValueError as exc:
            row2.update({'status': 'failed', 'message': str(exc)})

        rows_sys1.append(row1)
        rows_sys2.append(row2)

    print('Gaussian Elimination Experiment for Problem 2')
    print(f'dtype={np.dtype(dtype).name}, pivoting={pivoting}')
    print('epsilon = 10^(-2k), k=0..10')
    print()

    print_table('Table A: Problem 2(1)', rows_sys1)
    print_table('Table B: Problem 2(2)', rows_sys2)

    plot_results(rows_sys1, rows_sys2, output_path='visualization.png')
    print('Visualization generated: visualization.png')

    return rows_sys1, rows_sys2


if __name__ == '__main__':
    # single precision + no pivoting highlights instability
    run_experiment(dtype=np.float32, pivoting=False)





