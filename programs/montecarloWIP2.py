"""
monte_carlo.py
Experiments with forward modeling (despite the name!)
Ryan Bakko & Trey Wenger - August, September 2024

"""

import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c
import CurrentHIICode


# Simulate Observation
def simulate_observation(freqs, physical_params, nu_0, noise):
    true_spectrum = forward_model(freqs, physical_params, nu_0)
    observed_spectrum = add_noise(true_spectrum, loc=0.0, scale=noise)
    return observed_spectrum


# Adding noise function
def add_noise(spectrum, loc=0.0, scale=1.0):
    noise = np.random.normal(loc=loc, scale=scale, size=spectrum.shape)
    return spectrum + noise


# Define Forward Model
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
    T_e, EM, non_thermal_fwhm = physical_params
    freqs = freqs * u.GHz
    T_e = T_e * u.K
    nu_0 = nu_0 * u.GHz
    EM = EM * u.pc / u.cm**6
    non_thermal_fwhm = non_thermal_fwhm * u.km / u.s

    non_thermal_fwhm_freq = nu_0 * non_thermal_fwhm / c.c
    thermal_fwhm_freq = CurrentHIICode.fwhm(T_e, nu_0)
    fwhm_freq = np.sqrt(thermal_fwhm_freq**2.0 + non_thermal_fwhm_freq**2.0)
    tau_c = CurrentHIICode.continuum_optical_depth(EM, freqs, T_e)
    tau_l = CurrentHIICode.line_center_opacity(EM, fwhm_freq, T_e)
    tau_l_profile = CurrentHIICode.line_profile(freqs, nu_0, tau_l, fwhm_freq)
    brightness_temp = CurrentHIICode.brightness_temperature(T_e, tau_c + tau_l_profile)

    return brightness_temp.to("K").value


# Monte Carlo Least Squares Fitting
def monte_carlo_least_squares(
    freqs, true_physical_params, nu_0, noise, num_iterations=1000
):
    """
    Performs Monte Carlo least squares fitting multiple times on shuffled datasets.

    Parameters:
    num_iterations: int
        Number of Monte Carlo iterations.

    Returns:
    np.array
        Array of fitted parameters from each iteration.
    """
    fitted_params = []
    for i in range(num_iterations):
        print(i, end="\r")
        observed_spectrum = simulate_observation(
            freqs, true_physical_params, nu_0, noise
        )
        initial_guess = [7500.0, 950.0, 20.0]  # Initial guess
        result = least_squares(
            loss_function, initial_guess, args=(freqs, nu_0, observed_spectrum)
        )
        fitted_params.append(result.x)
    return np.array(fitted_params)


# Loss Function
def loss_function(physical_params, freqs, nu_0, observed_spectrum):
    predicted_spectrum = forward_model(freqs, physical_params, nu_0)
    return observed_spectrum - predicted_spectrum


# Posterior Predictive Samples
def posterior_predictive_samples(fitted_params, freqs, nu_0, noise, num_samples=100):
    """
    Generates posterior predictive samples based on fitted parameters.

    Parameters:
    fitted_params: np.array
        Array of fitted parameters from Monte Carlo simulation.
    freqs: np.array
        Frequencies for the model.
    nu_0: float
        Reference frequency.
    noise: float
        Noise level to add to the simulated spectra.
    num_samples: int
        Number of posterior predictive samples to generate.

    Returns:
    np.array
        Array of posterior predictive spectra.
    """
    posterior_samples = []
    for i in range(num_samples):
        random_params = fitted_params[np.random.randint(0, len(fitted_params))]
        simulated_spectrum = simulate_observation(freqs, random_params, nu_0, noise)
        posterior_samples.append(simulated_spectrum)
    return np.array(posterior_samples)


# Main Workflow
if __name__ == "__main__":
    freqs = np.linspace(1.4, 1.6, 100)  # Frequency range in GHz
    true_physical_params = [7500.0, 950.0, 20.0]  # Example true parameters
    nu_0 = 1.5  # Reference frequency
    noise = 0.1  # Example noise level

    # Run Monte Carlo simulation to get fitted parameters
    fitted_params = monte_carlo_least_squares(
        freqs, true_physical_params, nu_0, noise, num_iterations=500
    )

    # Generate posterior predictive samples
    posterior_samples = posterior_predictive_samples(
        fitted_params, freqs, nu_0, noise, num_samples=100
    )

    # Plot the posterior predictive samples
    plt.figure(figsize=(10, 5))
    for sample in posterior_samples:
        plt.plot(freqs, sample, color="gray", alpha=0.3)  # Plot posterior samples
    plt.title("Posterior Predictive Samples")
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Flux Density")
    plt.show()
