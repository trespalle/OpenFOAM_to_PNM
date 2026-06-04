import sys 
import numpy as np
from mpmath import mp
import scipy.stats as st
from scipy.optimize import curve_fit
from scipy.special import gamma as gamma_func

# Set the desired precision for calculations 
mp.dps = 10 # dps = decimal places


N_BINS_FOR_FIT = 50 
PERCENTILE_TO_IGNORE = 30.0 


if len(sys.argv) != 3:
    print("Usage: python calculate_fit.py <input_data_filepath> <output_filepath_base>")
    sys.exit(1)

# Get parameters from command line 
data_path = sys.argv[1]
base_path = sys.argv[2]

print(f"Loading raw data from {data_path}...")
try:
    # Load data
    data = np.loadtxt(data_path)
    data = data[data > 1e-9] 
    if len(data) == 0:
        raise ValueError("No valid data points found after filtering.")

    # Build histogram
    y_data, bin_edges = np.histogram(data, bins=N_BINS_FOR_FIT, density=True)
    x_data = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    # Filter data
    x_min_threshold = np.percentile(data, PERCENTILE_TO_IGNORE)
    
    # Create mask
    fit_mask = (x_data > x_min_threshold)
    
    x_fit = x_data[fit_mask]
    y_fit = y_data[fit_mask]
    
    if len(x_fit) < 2:
        raise ValueError(f"Not enough data points left to fit after filtering at {PERCENTILE_TO_IGNORE}th percentile.")
        
    print(f"Fitting on histogram bins where q > {x_min_threshold:.2e} (ignoring bottom {PERCENTILE_TO_IGNORE}%)")

    # Define fit function
    def gamma_pdf_for_fit(x, k, theta):
        if k <= 0 or theta <= 0:
            return np.full_like(x, np.inf)
        x_safe = np.clip(x, 1e-30, None)
        return (1.0 / (gamma_func(k) * (theta**k))) * (x_safe**(k-1)) * np.exp(-x_safe/theta)

    # Inital guess of parameters
    mu_guess = np.mean(data)
    var_guess = np.var(data)
    k_guess = mu_guess**2 / var_guess
    theta_guess = var_guess / mu_guess

    print("Fitting Gamma PDF using Least-Squares on filtered histogram...")
    # Execute LS fit
    popt, pcov = curve_fit(gamma_pdf_for_fit, x_fit, y_fit, p0=[k_guess, theta_guess])
    
    k, theta = popt
    mu = k * theta
    sigma2 = k * (theta**2)
    print(f"Least-Squares Fit: k = {k:.4f}, theta = {theta:.4f}")

except Exception as e:
    print(f"An error occurred: {e}. Exiting.")
    sys.exit(1)


# Define functions g(q) and f(q) 
def g(q, k, theta):
    return (1.0 / float(mp.gamma(k)) / theta**k) * q**(k-1) * np.exp(-q/theta)
    
def f(q, k, theta):
    if k <= 1: return 0.0
    return float(mp.gammainc(k - 1, a=q/theta)) / (float(mp.gamma(k)) * theta)

# Generate the points for the curve 
q_max = np.mean(data) + 3 * np.std(data)
q_values = np.linspace(q_max / 500.0, q_max, 500) 


p_values = np.zeros_like(q_values)
g_values = np.zeros_like(q_values)
f_values = np.zeros_like(q_values)

for i, q in enumerate(q_values):
    g_val = g(q, k, theta)
    f_val = f(q, k, theta)
    g_values[i] = g_val
    f_values[i] = f_val
    p_values[i] = (1.0/3.0) * g_val + (2.0/3.0) * f_val

np.savetxt(f"{base_path}_P_fit.dat", np.column_stack((q_values, p_values)), header="q P(q)")
np.savetxt(f"{base_path}_g_fit.dat", np.column_stack((q_values, g_values)), header="q g(q)")
np.savetxt(f"{base_path}_f_fit.dat", np.column_stack((q_values, f_values)), header="q f(q)")

print(f"Analytical curve data saved to {base_path}_*.dat files")