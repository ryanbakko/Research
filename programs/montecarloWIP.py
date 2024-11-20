# Import necessary libraries
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c

# Import any other specific models from custom modules
import CurrentHIICode


# Define the forward model
def forward_model(freqs, physical_params, nu_0):
    """
       Predicts an observed radio recombination line and radio continuum
       spectra based on input frequencies and physical parameters.

    Inputs:
           freqs :: 1-D array of scalars (GHz)
               Frequencies at which to evaluate the spectrum
           physical_params :: 3-length array of scalars
               physical_params[0] = Electron temperature (K)
               physical_params[1] = Emission measure (pc cm-6)
               physical_params[2] = Non-thermal FWHM line width in velocity units (km s-1)
           nu_0 :: scalar (GHz)
               RRL rest frequency

      Returns:
           spectrum :: 1-D array of scalars (K)
               Brightness temperature spectrum
    """
    # EM, Te, ne = physical_params
    T_e, EM, non_thermal_fwhm = physical_params

    # add units
    freqs = freqs * u.GHz
    T_e = T_e * u.K
    nu_0 = nu_0 * u.GHz
    EM = EM * u.pc / u.cm**6
    non_thermal_fwhm = non_thermal_fwhm * u.km / u.s

    # evaluate physics
    non_thermal_fwhm_freq = nu_0 * non_thermal_fwhm / c.c
    thermal_fwhm_freq = CurrentHIICode.fwhm(T_e, nu_0)
    fwhm_freq = np.sqrt(thermal_fwhm_freq**2.0 + non_thermal_fwhm_freq**2.0)
    tau_c = CurrentHIICode.continuum_optical_depth(EM, freqs, T_e)
    tau_l = CurrentHIICode.line_center_opacity(EM, fwhm_freq, T_e)
    tau_l_profile = CurrentHIICode.line_profile(freqs, nu_0, tau_l, fwhm_freq)
    brightness_temp = CurrentHIICode.brightness_temperature(T_e, tau_c + tau_l_profile)

    return brightness_temp.to("K").value


# Define function to add noise to the spectrum
def add_noise(spectrum, loc=0.0, scale=1.0):
    noise = np.random.normal(loc=loc, scale=scale, size=spectrum.shape)
    return spectrum + noise


# Simulate observations based on the forward model
def simulate_observation(freqs, true_physical_params, nu_0, noise):
    """
    Simulate a noisy spectrum based on the true physical parameters.
    """
    true_spectrum = forward_model(freqs, true_physical_params, nu_0)
    observed_spectrum = add_noise(true_spectrum, loc=0.0, scale=noise)
    return observed_spectrum


# Perform Monte Carlo least squares fitting
def monte_carlo_least_squares(
    freqs, true_physical_params, nu_0, noise, num_iterations=1000
):
    """
    Perform multiple least squares fits using simulated data to generate distributions.
    """
    fitted_params = []

    for _ in range(num_iterations):
        # Simulate observed data
        observed_spectrum = simulate_observation(
            freqs, true_physical_params, nu_0, noise
        )
        # Set up initial guess for fitting
        initial_guess = [7500.0, 950.0, 20.0]
        # Perform fitting using least squares
        params, _ = fit_spectrum(freqs, observed_spectrum, nu_0, initial_guess)
        fitted_params.append(params)

    return np.array(fitted_params)


# Define the fitting function
def fit_spectrum(freqs, data, nu_0, initial_guess):
    """
    Fits the spectrum model to the observed data using non-linear least squares.
    """

    def residuals(phys_params, freqs, data, nu_0):
        """
        Residuals function for least squares fitting.
        """
        model = forward_model(freqs, phys_params, nu_0)  # Generate the model spectrum
        return model - data  # Residuals between the model and observed data

    # Perform the least squares fitting
    result = least_squares(residuals, initial_guess, args=(freqs, data, nu_0))

    # Calculate uncertainties from the covariance matrix
    jacobian = result.jac  # Jacobian matrix from the optimization result
    covariance = np.linalg.inv(
        jacobian.T @ jacobian
    )  # Inverse of the Hessian approximation
    uncertainties = np.sqrt(
        np.diag(covariance)
    )  # Standard deviations of the fitted parameters

    return result.x, uncertainties  # Return fitted parameters and uncertainties


# Define a function to plot the Monte Carlo results
def plot_monte_carlo_results(fitted_params, param_names):
    """
    Plots the Monte Carlo results, including posterior distributions and uncertainty intervals.
    """
    num_params = fitted_params.shape[1]
    fig, axes = plt.subplots(num_params, 2, figsize=(10, 5 * num_params))

    for i in range(num_params):
        param_vals = fitted_params[:, i]

        # Plot histogram for posterior distribution
        axes[i, 0].hist(param_vals, bins=30, density=True, alpha=0.75)
        axes[i, 0].set_title(f"Posterior Distribution for {param_names[i]}")

        # Plot uncertainty interval (mean and standard deviation)
        mean = np.mean(param_vals)
        std_dev = np.std(param_vals)
        axes[i, 1].errorbar(1, mean, yerr=std_dev, fmt="o")
        axes[i, 1].set_title(f"Uncertainty Interval for {param_names[i]}")

    plt.tight_layout()
    plt.show()


# Add a main section to test the workflow
if __name__ == "__main__":
    # Example frequency range and true parameters
    freqs = np.linspace(6.99, 7.01, 100)  # Frequency range in GHz
    true_physical_params = [8000.0, 1000.0, 25.0]  # Example physical parameters
    nu_0 = 7.0  # Reference frequency
    noise = 0.05  # Noise level for the observations

    # Run Monte Carlo simulation to get fitted parameters
    fitted_params = monte_carlo_least_squares(
        freqs, true_physical_params, nu_0, noise, num_iterations=500
    )

    # Define parameter names for plotting
    param_names = ["T_e", "EM", "v_turb"]

    # Plot the Monte Carlo results
    plot_monte_carlo_results(fitted_params, param_names)
