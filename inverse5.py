"""
inverse5.py
Experiments with forward modeling (despite the name!)
Ryan Bakko & Trey Wenger - August 2024

Worked on calculating uncertainties in parameters and drawing the uncertainty interval on my plot
TODO: Debug and check if I'm on the right track
"""
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c
import CurrentHIICode.py


# Step 1 of Monte Carlo: Shuffle Data
def shuffle_data(data):
    """
    Shuffles the input data.
    
    Parameters:
    data: np.array
        Array of data points to shuffle.
    
    Returns:
    np.array
        Shuffled data.
    """
    shuffled_data = np.copy(data)
    np.random.shuffle(shuffled_data)
    return shuffled_data

# Step 2: Monte Carlo Least Squares Fitting
def monte_carlo_least_squares(data_x, data_y, model_func, p0, num_iterations=1000):
    """
    Performs Monte Carlo least squares fitting multiple times on shuffled datasets.
    
    Parameters:
    data_x: np.array
        Independent variable data points (e.g., frequencies).
    data_y: np.array
        Dependent variable data points (e.g., brightness temperature).
    model_func: callable
        The forward model function to fit.
    p0: np.array
        Initial guess for the parameters.
    num_iterations: int
        Number of Monte Carlo iterations.
        
    Returns:
    np.array
        Array of fitted parameters from each iteration.
    """
    fitted_params = []

    for i in range(num_iterations):
        # Shuffle the data (here I shuffle the y-data, but it could be x-data as well)
        shuffled_y = shuffle_data(data_y)
        
        # Perform least squares fitting
        result = least_squares(lambda params: model_func(data_x, params) - shuffled_y, p0)
        
        # Save the fitted parameters
        fitted_params.append(result.x)

    return np.array(fitted_params)

# Step 3: Plotting the results
def plot_monte_carlo_results(fitted_params, param_names):
    """
    Plots the Monte Carlo results, including posterior distributions and uncertainty intervals.
    
    Parameters:
    fitted_params: np.array
        Array of fitted parameters from Monte Carlo simulation.
    param_names: list
        List of parameter names corresponding to the columns in fitted_params.
    """
    num_params = fitted_params.shape[1]
    
    fig, axes = plt.subplots(num_params, 2, figsize=(10, 5 * num_params))
    
    for i in range(num_params):
        param_vals = fitted_params[:, i]
        
        # Plot histogram for posterior distribution
        axes[i, 0].hist(param_vals, bins=30, density=True, alpha=0.75)
        axes[i, 0].set_title(f'Posterior Distribution for {param_names[i]}')
        
        # Plot uncertainty interval (mean and standard deviation)
        mean = np.mean(param_vals)
        std_dev = np.std(param_vals)
        axes[i, 1].errorbar(1, mean, yerr=std_dev, fmt='o')
        axes[i, 1].set_title(f'Uncertainty Interval for {param_names[i]}')

    plt.tight_layout()
    plt.show()

# End of New Monte Carlo additions

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
    #EM, Te, ne = physical_params
    T_e, EM, non_thermal_fwhm = physical_params

    # add units
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

    return brightness_temp.to('K').value


#Adding noise function
def add_noise(spectrum, loc=0.0, scale=1.0):
    noise = np.random.normal(loc=loc, scale=scale, size=spectrum.shape)
    return spectrum + noise

# Loss function
def loss_function(physical_params, freqs, nu_0, observed_spectrum):
    predicted_spectrum = forward_model(freqs, physical_params, nu_0)
    residuals = observed_spectrum - predicted_spectrum
    return residuals

def residuals(physical_params, freqs, data, nu_0):
    """
    Computes the residuals between the observed data and the model prediction.
    
    Parameters:
    physical_params (list): List of parameters to fit.
    freqs (np.ndarray): Array of frequency values.
    data (np.ndarray): Observed data to fit the model against.
    nu_0 (float): Reference frequency.
    
    Returns:
    np.ndarray: Residuals between the model and observed data.
    """
    model = forward_model(freqs, physical_params, nu_0)
    return model - data

def fit_spectrum(freqs, data, nu_0, initial_guess):
    """
    Fits the spectrum model to the observed data using non-linear least squares.
    
    Parameters:
    freqs (np.ndarray): Array of frequency values.
    data (np.ndarray): Observed data to fit the model against.
    nu_0 (float): Reference frequency.
    initial_guess (list): Initial guess for the physical parameters.
    
    Returns:
    tuple: Fitted parameters and their uncertainties.
    """
    # Prepare your data (generate synthetic observations)
    freqs = np.linspace(6.99, 7.01, 1000) * u.GHz  # Frequency range for the spectrum - New
    true_physical_params = [8000.0, 1000.0, 25.0]  # Replace with true physical parameters
    nu_0 = 7.0 # RRL rest frequency (GHz)
    true_spectrum = forward_model(freqs, true_physical_params, nu_0)
    noise = 0.001
    observed_spectrum = add_noise(true_spectrum, loc=0.0, scale=noise)  # Adjust scale

    # Set up initial parameters
    initial_guess = [7500.0, 950.0, 20.0]  # Initial guess for the physical parameters
        
    #Perform the least squares fitting
    # #result = least_squares(residuals, initial_guess, args=(freqs, data, nu_0))
    # Run least squares fitting
    result = least_squares(loss_function, initial_guess, args=(freqs, nu_0, observed_spectrum,))

    # Calculate uncertainties from the covariance matrix
    jacobian = result.jac
    covariance = np.linalg.inv(jacobian.T @ jacobian)
    uncertainties = np.sqrt(np.diag(covariance))
    return result.x, uncertainties


def plot_fit_with_uncertainty(freqs, data, nu_0, fitted_params, uncertainties):
    """
    Plots the fitted model along with the uncertainty intervals.
    
    Parameters:
    freqs (np.ndarray): Array of frequency values.
    data (np.ndarray): Observed data to fit the model against.
    nu_0 (float): Reference frequency.
    fitted_params (list): Fitted physical parameters.
    uncertainties (np.ndarray): Uncertainties in the fitted parameters.
    
    Returns:
    None
    """
    # Generate the model data using the fitted parameters
    fitted_model = forward_model(freqs, fitted_params, nu_0)
    
    # Plot the data and the fitted model
    plt.figure(figsize=(10, 6))
    plt.plot(freqs, data, 'o', label='Observed Data')
    plt.plot(freqs, fitted_model, '-', label='Fitted Model', color='red')

    # Plot the uncertainty interval
    upper_bound = forward_model(freqs, fitted_params + uncertainties, nu_0)
    lower_bound = forward_model(freqs, fitted_params - uncertainties, nu_0)
    plt.fill_between(freqs, lower_bound, upper_bound, color='red', alpha=0.3, label='Uncertainty Interval')

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Observed Spectrum')
    plt.title('Fitted Spectrum with Uncertainty Interval')
    plt.legend()
    plt.show()

# Example usage (this should be in a separate script or a testing function)
if __name__ == "__main__":
    # Example data
    freqs = np.linspace(1e9, 1.5e9, 100)  # Example frequency range
    true_params = [1e4, 1e3, 1e2]  # True EM, Te, ne
    nu_0 = 1e9  # Reference frequency
    data = forward_model(freqs, true_params, nu_0) + np.random.normal(0, 0.05, freqs.shape)  # Add some noise
    
    # Initial guess for fitting
    initial_guess = [1e3, 1e2, 1e1]
    
    # Perform fitting
    fitted_params, uncertainties = fit_spectrum(freqs, data, nu_0, initial_guess)
    
    # Print the results
    print("True EM, Te, ne:", true_params)
    print("Fitted EM, Te, ne:", fitted_params)
    print("Uncertainties in fitted parameters:", uncertainties)
    
    # Plot the results with uncertainties
    plot_fit_with_uncertainty(freqs, data, nu_0, fitted_params, uncertainties)
    
