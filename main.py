import numpy as np
import matplotlib.pyplot as plt

def progonka(rhs):
	"""
	Решение методом прогонки (быстрее чем Гаусик)
	"""
	rhs = np.asarray(rhs, dtype=float)
	n = len(rhs)

	lower = np.full(n - 1, -1.0)
	main = np.full(n, 2.0)
	upper = np.full(n - 1, -1.0)

	d = rhs.copy()

	for k in range(1, n):
		coef = lower[k - 1] / main[k - 1]
		main[k] -= coef * upper[k - 1]
		d[k] -= coef * d[k - 1]

	sol = np.zeros(n)
	sol[-1] = d[-1] / main[-1]

	for k in range(n - 2, -1, -1):
		sol[k] = (d[k] - upper[k] * sol[k + 1]) / main[k]

	return sol


def exact_u(x):
	return np.sin(4 * x) * np.cos(3 * x)


def second_derivative(x):
	return -0.5 * (49 * np.sin(7 * x) + np.sin(x))


def build_rhs(N):
	step = 1.0 / N
	vec = np.zeros(N - 1)

	for j in range(1, N):
		point = j * step
		vec[j - 1] = -second_derivative(point)

	return vec


def solve_problem(N, left_val, right_val):
	rhs = build_rhs(N)
	step = 1.0 / N

	rhs *= step ** 2

	rhs[0] += left_val
	rhs[-1] += right_val

	inner = progonka(rhs)

	result = np.zeros(N + 1)
	result[0] = left_val
	result[1:-1] = inner
	result[-1] = right_val

	return result


def test_convergence(grid_sizes, left_val, right_val):
	errors_L2 = []
	errors_C = []
	steps = []

	for N in grid_sizes:
		step = 1.0 / N

		approx = solve_problem(N, left_val, right_val)

		grid = np.linspace(0, 1, N + 1)
		exact = exact_u(grid)

		diff = approx - exact

		err_L2 = np.sqrt(np.sum(diff ** 2) * step)
		err_C = np.max(np.abs(diff))

		errors_L2.append(err_L2)
		errors_C.append(err_C)
		steps.append(step)

	steps = np.array(steps)

	plt.figure(figsize=(8, 6))

	plt.loglog(steps, errors_L2, marker="o", label="L2")
	plt.loglog(steps, errors_C, marker="s", label="C")
	plt.loglog(steps, steps ** 2, linestyle="--", linewidth=2, label="O(h²)")

	plt.xlabel("шаг сетки h")
	plt.ylabel("величина ошибки")
	plt.title("Сходимость разностной схемы")

	plt.grid(True, which="both")
	plt.legend()

	plt.savefig("graph_convergence.png", dpi=100)
	plt.show()

	return errors_L2, errors_C, steps


if __name__ == "__main__":
	grids = [2 ** k for k in range(1, 16)]

	left = 0.0
	right = np.sin(4) * np.cos(3)

	test_convergence(grids, left, right)
