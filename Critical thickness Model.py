#-----------------------------------------------------------#
#                                                           #
# - By Lawrence Hopper                                      #
# - 22/02/2025                                              #
# - Calculates the critical thickness of the user selected  #
#   ternary epilayer on a binary substrate, plots the data  #
#   and saves it to a text file for further formatting      #
#                                                           #
#-----------------------------------------------------------#


import pandas as pd # For table analysis
import math #Math operators
import matplotlib.pyplot as plt #Plots graphs

# Load the Excel file
file_path = "Data.xlsx"  # The dataset containing properties of binary compounds
df = pd.read_excel(file_path, sheet_name="Data")  # Reads the "Data" sheet from the Excel file

# User Inputs
substrate_material = input("Enter the substrate material: ")  # User provides substrate material
A = input("Enter the first element in the ternary compound (A): ")  # First element of ternary compound
B = input("Enter the second element in the ternary compound (B): ")  # Second element of ternary compound
C = input("Enter the common element in the ternary compound (C): ")  # Common element in the ternary compound

# Form ternary compound names by combining elements
AC_compound = A + C  # Example: InAs
BC_compound = B + C  # Example: GaAs

# Function to fetch material properties from the dataset
def get_property(compound, property_name):
    result = df.loc[df["Compound"] == compound, property_name]
    if not result.empty:
        return float(result.values[0])  # Converts result to float for calculations
    raise ValueError(f"{property_name} for {compound} not found in the dataset!")  # Error handling if compound is missing

#Fetches lattice constants for the substrate, AC, and BC compounds
lattice_constant_substrate = get_property(substrate_material, "a (Å)")
lattice_constant_AC = get_property(AC_compound, "a (Å)")
lattice_constant_BC = get_property(BC_compound, "a (Å)")

#Fetches elastic constants (c11 and c12) for Poisson's ratio calculation
c11_AC = get_property(AC_compound, "c₁₁ (10¹¹ N/m²)")
c12_AC = get_property(AC_compound, "c₁₂ (10¹¹ N/m²)")
c11_BC = get_property(BC_compound, "c₁₁ (10¹¹ N/m²)")
c12_BC = get_property(BC_compound, "c₁₂ (10¹¹ N/m²)")

#Generates x values from 0 to 1 in steps of 0.01 (for ternary compound compositions)
x_values = [round(i * 0.01, 2) for i in range(101)]

# Compute lattice constants of the epilayer using linear interpolation
lattice_constants = [
    x * lattice_constant_AC + (1 - x) * lattice_constant_BC for x in x_values
]

# Compute Poisson's ratio (v) using a weighted sum
poisson_ratios = [
    (x * c12_AC + (1 - x) * c12_BC) / (x * c11_AC + (1 - x) * c11_BC + x * c12_AC + (1 - x) * c12_BC)
    for x in x_values
]

# Function to compute critical layer thickness
def compute_critical_thickness(lattice_epilayer, lattice_substrate, poisson_ratio):
    
    fvalue = (lattice_epilayer - lattice_substrate) / lattice_substrate  # Strain calculation
    # Handle zero or extremely small strain values to avoid division errors
    if abs(fvalue) < 1e-6:
        return float('nan')  # Return NaN if strain is too small
    b = 4  # Burgers vector approximation (unit Å)
    thickness = 50  # Initial guess for thickness (Å)

    # Iterative calculation for convergence
    for _ in range(100):  
        if thickness <= b:  # Prevent division errors when thickness is too small
            return float('nan')
        
        log_term = math.log(max(thickness / b, 1e-6))  # Ensure valid log input

        # People and bean equation for critical thickness
        new_thickness = ((1 - poisson_ratio) / (1 + poisson_ratio)) * (1 / (16 * math.pi * math.sqrt(2))) \
                        * ((b ** 2) / lattice_epilayer) * (1 / fvalue ** 2) * log_term
        
        if abs(new_thickness - thickness) < 1e-5:  # Convergence condition
            return new_thickness
        
        thickness = new_thickness  # Update thickness for next iteration

    return thickness  # Return final calculated thickness

# Compute critical layer thicknesses for different composition fractions
critical_layer_thicknesses = [
    compute_critical_thickness(lattice_constants[i], lattice_constant_substrate, poisson_ratios[i])
    for i in range(len(x_values))
]

# Create a DataFrame to store results
results_df = pd.DataFrame({
    "Composition Fraction (x)": x_values,
    "Lattice Constant of Epilayer (Å)": lattice_constants,
    "Critical Layer Thickness (Å)": critical_layer_thicknesses
})

# Save results to a text file
txt_filename = "lattice_critical_thickness_results.txt"
results_df.to_csv(txt_filename, sep=",", index=False, float_format="%.5f")
print(f"Results saved to {txt_filename}")

# Function to plot results
def plot_graphs():
    """
    Plots the critical layer thickness as a function of composition fraction.
    """
    plt.figure(figsize=(8, 5))  # Set figure size
    plt.plot(x_values, critical_layer_thicknesses, marker='s', linestyle='-', color='r', label="Critical Layer Thickness (Å)")
    plt.xlabel("Composition Fraction (x)")  # X-axis label
    plt.ylabel("Critical Layer Thickness (Å)")  # Y-axis label
    plt.title("Critical Layer Thickness vs. Composition Fraction")  # Graph title
    plt.yscale("log")  # Set logarithmic scale for better visualization
    plt.grid(True)  # Add grid for readability
    plt.legend()  # Show legend
    plt.show()  # Display plot

# Call function to display the plot
plot_graphs()








