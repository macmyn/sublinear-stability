import numpy as np
# import pytest
import threading
from scipy import integrate

def derivative(t, N, r, a, K, exponent_intra, exponent_inter):
    N = np.array(N)  # Ensure N is a NumPy array
    sum_aN_off_diag = np.dot(a - np.diag(np.diag(a)), N ** exponent_inter)
    sum_aN_diag = np.diag(a) * np.where(N > 1e-10, N ** exponent_intra, 0)
    sum_aN = np.real(sum_aN_off_diag + np.sign(exponent_intra) * sum_aN_diag)
    dotN = N * (np.sign(exponent_intra) * r - sum_aN) / K
    dotN = np.where(N < 0, 0, dotN)
    return dotN

# Generate interaction matrix
def generate_interaction_matrix(num_species, connectivity, mean, sd):
    a = np.zeros((num_species, num_species))
    np.fill_diagonal(a, 1)
    for i in range(num_species):
        for j in range(num_species):
            if i != j and np.random.random() < connectivity:
                a[i, j] = np.random.normal(loc=mean, scale=sd)
    return a

# Get eigenvalues of the Jacobian matrix
def get_eigenvalues(N, r, K, a, exponent_intra, exponent_inter):
    num_species = len(N)
    jacobian_matrix = np.zeros((num_species, num_species))
    for i in range(num_species):
        for j in range(num_species):
            if i == j:
                sum_inter = sum(a[i, k] * (N[k]) ** exponent_inter / K[i] for k in range(num_species)) - a[i, i] * (N[i]) ** exponent_inter / K[i]
                jacobian_matrix[i, j] = np.sign(exponent_intra) * (r[i] - (1 + exponent_intra) * a[i, i] * np.real(N[i] ** exponent_intra) / K[i] - sum_inter)
            else:
                jacobian_matrix[i, j] = -a[i, j] * exponent_inter * N[i] * np.real(N[j] ** (exponent_inter - 1)) / K[i]
    if not np.any(np.isnan(jacobian_matrix)) and not np.any(np.isinf(jacobian_matrix)):
        return np.linalg.eigvals(jacobian_matrix)
    else:
        return np.nan

# Timeout handler using threading
class ODESolverWithTimeout:
    def __init__(self, fun, t_span, y0, args=(), max_time=5, t_eval=None, precision=1e-10, max_step=0.1):
        self.fun = fun
        self.t_span = t_span
        self.y0 = y0
        self.args = args
        self.max_time = max_time
        self.t_eval = t_eval
        self.precision = precision
        self.max_step = max_step
        self.result = None
        self.thread = None

    def solve(self):
        self.thread = threading.Thread(target=self._solve_ode)
        self.thread.start()
        self.thread.join(timeout=self.max_time)
        if self.thread.is_alive():
            print("Solver timed out")
            self.result = None
        return self.result

    def _solve_ode(self):
        try:
            sol = integrate.solve_ivp(self.fun, self.t_span, self.y0, args=self.args, t_eval=self.t_eval, method='BDF', atol=self.precision, rtol=self.precision, max_step=self.max_step)
            self.result = sol
        except Exception as e:
            print(f"Solver failed: {e}")
            self.result = None