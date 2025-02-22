#-----------------------------------------------------------#
#                                                           #
# - By Lawrence Hopper                                      #
# - 21/02/2025                                              #
# - This code uses the Krijn model to simulate a binary     #
#   epilayer on binary substrate and ternary epilayer       #
#   on a binary substrate and calculates band shift         #
#                                                           #
#-----------------------------------------------------------#


import pandas as pd # For table analysis
from tabulate import tabulate #Formats tables neatly
import numpy as np #Math operators
import matplotlib.pyplot as plt #Plots graphs


#reads the excel file and displays it in the terminal
df = pd.read_excel("Data.xlsx", sheet_name="Data")    #Reads Excel file containing the binary compound data
pd.set_option("display.max_rows", None)               #Show all rows
pd.set_option("display.max_columns", None)            #Show all columns
print(tabulate(df, headers="keys", tablefmt="grid"))  #Print the full Table


#asks the user for choice of simulation either b-b or t-b
choice = int(input("Which would you like to simulate?\n1: Binary Compound\n2: Ternary Compound\n "))

if choice == 1:
    print("Simulating the binary compound")
    substrate = input("Please enter the compound to be used as the Substrate: ")# Ask the user for the binary substrate compound
    epilayer = input("Please enter the compound to be used as the Epilayer: ")# Ask the user for the binary epilayer compound
    print(f"The selected Substrate is: {substrate}")
    print(f"The selected Epilayer is: {epilayer}")


    try:
        substrate_row = df[df["Compound"] == substrate].iloc[0]  # First matching row for Substrate
        substrate_data = substrate_row.to_dict()                 # Convert to a dictionary
        epilayer_row = df[df["Compound"] == epilayer].iloc[0]    # First matching row for Epilayer
        epilayer_data = epilayer_row.to_dict()                   # Convert to a dictionary
        print(f"Substrate data for {substrate}:")
        print(tabulate([substrate_data.values()], headers=substrate_data.keys(), tablefmt="grid"))
        print(f"Epilayer data for {epilayer}:")
        print(tabulate([epilayer_data.values()], headers=epilayer_data.keys(), tablefmt="grid"))

        #The data for the Epilayer and Substrate is printed to the terminal


    except IndexError:
        print("Error: One of the compounds could not be found. Please check your inputs.")
        #If the user inputs cannot be found, an error is printed
        

                ####################KRIJN EQUATIONS####################

        
    substrate_a_parallel = substrate_row["a (Å)"] #Retrieves the parallel lattice constants
    epilayer_a_parallel = substrate_row["a (Å)"]  
    
    substrate_D001 = 2 * (substrate_row["c₁₂ (10¹¹ N/m²)"] / substrate_row["c₁₁ (10¹¹ N/m²)"]) #Calculates D001 for the sub
    epilayer_D001 = 2 * (epilayer_row["c₁₂ (10¹¹ N/m²)"] / epilayer_row["c₁₁ (10¹¹ N/m²)"])    
    
    substrate_a_perp = epilayer_row["a (Å)"] * (1 - (substrate_D001 * ((substrate_a_parallel/epilayer_row["a (Å)"] ) - 1)))
    epilayer_a_perp = epilayer_row["a (Å)"] * (1 - (epilayer_D001 * ((epilayer_a_parallel/epilayer_row["a (Å)"] ) - 1)))
    #Calculates perpendicular lattice constants
    
    substrate_epsilon_perp = (substrate_a_perp/epilayer_row["a (Å)"]) - 1 #Perpendicular strain calculation
    epilayer_epsilon_perp = (epilayer_a_perp/epilayer_row["a (Å)"]) - 1   
    
    substrate_epsilon_parallel = (substrate_a_parallel/epilayer_row["a (Å)"]) - 1 #Parallel strain calculation
    epilayer_epsilon_parallel = (epilayer_a_parallel/epilayer_row["a (Å)"]) - 1   
    
    substrate_hydro_average = substrate_row["aᵥ (eV)"] * (2*substrate_epsilon_parallel + substrate_epsilon_perp)#Hydrostatic average
    epilayer_hydro_average = epilayer_row["aᵥ (eV)"] * (2*epilayer_epsilon_parallel + epilayer_epsilon_perp)
    
    substrate_hydro_conduction = substrate_row["aₑ(Γ) (eV)"] * (2*substrate_epsilon_parallel + substrate_epsilon_perp)#Hydrostatic Conduction
    epilayer_hydro_conduction = epilayer_row["aₑ(Γ) (eV)"] * (2*epilayer_epsilon_parallel + epilayer_epsilon_perp)
    
    substrate_strainshift = 2*substrate_row["b (eV)"]*(substrate_epsilon_perp-substrate_epsilon_parallel)#Strain shift calculations
    epilayer_strainshift = 2*epilayer_row["b (eV)"]*(epilayer_epsilon_perp-epilayer_epsilon_parallel)
    
    substrate_hh = -0.5 * substrate_strainshift #Heavy hole calculations
    epilayer_hh = -0.5 * epilayer_strainshift
    
    substrate_lh = -0.5 * substrate_row["Δ₀ (eV)"] + 0.25 * substrate_strainshift + \
                   0.5 * ((substrate_row["Δ₀ (eV)"])**2 + substrate_row["Δ₀ (eV)"] * substrate_strainshift + 2.25 * (substrate_strainshift)**2) ** 0.5

    epilayer_lh = -0.5 * epilayer_row["Δ₀ (eV)"] + 0.25 * epilayer_strainshift + \
                  0.5 * ((epilayer_row["Δ₀ (eV)"])**2 + epilayer_row["Δ₀ (eV)"] * epilayer_strainshift + 2.25 * (epilayer_strainshift)**2) ** 0.5

    substrate_so = -0.5 * substrate_row["Δ₀ (eV)"] + 0.25 * substrate_strainshift - \
                   0.5 * ((substrate_row["Δ₀ (eV)"])**2 + substrate_row["Δ₀ (eV)"] * substrate_strainshift + 2.25 * (substrate_strainshift)**2) ** 0.5

    epilayer_so = -0.5 * epilayer_row["Δ₀ (eV)"] + 0.25 * epilayer_strainshift - \
                  0.5 * ((epilayer_row["Δ₀ (eV)"])**2 + epilayer_row["Δ₀ (eV)"] * epilayer_strainshift + 2.25 * (epilayer_strainshift)**2) ** 0.5


    substrate_valance = substrate_row["Eᵥ,ₐᵥ (eV)"] + (substrate_row["Δ₀ (eV)"]/3) + substrate_hydro_average +max(substrate_hh,substrate_lh)
    epilayer_valance = epilayer_row["Eᵥ,ₐᵥ (eV)"] + (epilayer_row["Δ₀ (eV)"]/3) + epilayer_hydro_average +max(epilayer_hh,epilayer_lh)
    #Valence band calculations
    
    substrate_conduction = substrate_row["Eᵥ,ₐᵥ (eV)"] + (substrate_row["Δ₀ (eV)"]/3) + substrate_row["Eg(Γ) (eV)"] + substrate_hydro_conduction
    epilayer_conduction = epilayer_row["Eᵥ,ₐᵥ (eV)"] + (epilayer_row["Δ₀ (eV)"]/3) + epilayer_row["Eg(Γ) (eV)"] + epilayer_hydro_conduction
    #Conduction band calculations
    
    substrate_bandgap = substrate_conduction - substrate_valance #Bandgap calculations
    epilayer_bandgap = epilayer_conduction - epilayer_valance
    
    print(f"The substrate bandgap is: {substrate_bandgap}") #Outputs the bandgaps to the terminal
    print(f"The epilayer bandgap is: {epilayer_bandgap}")

    
if choice == 2:
    print("Simulating the ternary compound")
    df1 = pd.read_excel("Data.xlsx", sheet_name="Ternary Data") #Reads Excel file containing the bowing parameters
    pd.set_option("display.max_rows", None)                     #Show all rows
    pd.set_option("display.max_columns", None)                  #Show all columns
    print(tabulate(df1, headers="keys", tablefmt="grid"))       #Print the full bowing parameter table
    substrate = input("Please enter the compound to be used as the Substrate: ") # Ask the user for the binary Substrate compound
    epilayer = int(input("Please enter the ternary compound to be used as the Epilayer, 1-14: "))#Asks user for the ternary epilayer compound
    #The ternary compounds are listed 1 to 14 so just the number is needed when selecting the epilaye
    print(f"The selected Substrate is: {substrate}") #Displays the substrate and epilayer selected to the terminal
    print(f"The selected Epilayer is: {epilayer}")

    
    try:
        substrate_row = df[df["Compound"] == substrate].iloc[0]  # First matching row for Substrate is stored 
        substrate_data = substrate_row.to_dict()  # Convert to a dictionary
        t_epilayer_row = df1[df1["Compound Number"] == epilayer].iloc[0]    # First matching row for Epilayer is stored 
        t_epilayer_data = t_epilayer_row.to_dict()  # Convert to a dictionary

        print(f"Substrate data for {substrate}:")        # Prints rows to terminal
        print(tabulate([substrate_data.values()], headers=substrate_data.keys(), tablefmt="grid"))
        print(f"Epilayer data for {epilayer}:")
        print(tabulate([t_epilayer_data.values()], headers=t_epilayer_data.keys(), tablefmt="grid"))

    except IndexError: # If one of the inputs is wrong or not found, an error message is shown
        print("Error: One of the compounds could not be found. Please check your inputs.")

    x_values = np.arange(0.99, -0.01, -0.01) #x values are defined in an array, going from 0.99 to 0.01 in increments of 0.01
    x_values = np.round(x_values, decimals=2)#Rounds the values to 2 dp
    binary_AB = input("What is the Epilayer B(AB): ")
    binary_AB_row = df[df["Compound"] == binary_AB].iloc[0]  # First matching row for AC
    binary_AB_data = binary_AB_row.to_dict()  # Convert to a dictionary
    binary_AC = input("What is the Epilayer B(AC): ")
    binary_AC_row = df[df["Compound"] == binary_AC].iloc[0]  # First matching row for BC
    binary_AC_data = binary_AC_row.to_dict()  # Convert to a dictionary
    #The binary compounds making up the ternary alloy need to be defined so they can be extracted from the excel file containing the data
    print(f"Substrate data for {binary_AB}:") #Prints the data for AB and AC for reference
    print(tabulate([binary_AB_data.values()], headers=binary_AB_data.keys(), tablefmt="grid"))
    print(f"Epilayer data for {binary_AC}:")
    print(tabulate([binary_AC_data.values()], headers=binary_AC_data.keys(), tablefmt="grid"))


    #This small section calculates the average of the ternary alloy parameters dependent on the composition
    a_ternary = x_values * binary_AB_row["a (Å)"] + (1-x_values) * binary_AC_row["a (Å)"] 
    c12_ternary = binary_AB_row["c₁₂ (10¹¹ N/m²)"] * x_values + (1-x_values) * binary_AC_row["c₁₂ (10¹¹ N/m²)"]
    c11_ternary = binary_AB_row["c₁₁ (10¹¹ N/m²)"] * x_values + (1-x_values) * binary_AC_row["c₁₁ (10¹¹ N/m²)"]
    Cvav_ternary = 3 * (binary_AB_row["aᵥ (eV)"]-binary_AC_row["aᵥ (eV)"])*((binary_AB_row["a (Å)"]-binary_AC_row["a (Å)"])/substrate_row["a (Å)"])
    Evav_ternary = binary_AB_row["Eᵥ,ₐᵥ (eV)"] * x_values + (1-x_values) * binary_AC_row["Eᵥ,ₐᵥ (eV)"] + x_values*(x_values-1)* Cvav_ternary
    delta0_ternary = binary_AB_row["Δ₀ (eV)"] * x_values + (1-x_values) * binary_AC_row["Δ₀ (eV)"] + x_values*(x_values-1)*t_epilayer_row["C(0)"]
    Eg_ternary = binary_AB_row["Eg(Γ) (eV)"] * x_values + (1-x_values) * binary_AC_row["Eg(Γ) (eV)"] + x_values*(x_values-1)*t_epilayer_row["C(Eg)"]
    av_ternary = binary_AB_row["aᵥ (eV)"] * x_values + (1-x_values) * binary_AC_row["aᵥ (eV)"]
    ac_ternary = binary_AB_row["aₑ(Γ) (eV)"] * x_values + (1-x_values) * binary_AC_row["aₑ(Γ) (eV)"]
    b_ternary = binary_AB_row["b (eV)"] * x_values + (1-x_values) * binary_AC_row["b (eV)"]

                #########################KRIJN EQUATIONS#########################

    substrate_a_parallel = substrate_row["a (Å)"]
    epilayer_a_parallel = substrate_row["a (Å)"]
    substrate_D001 = 2 * (substrate_row["c₁₂ (10¹¹ N/m²)"] / substrate_row["c₁₁ (10¹¹ N/m²)"])
    epilayer_D001 = 2 * (c12_ternary / c11_ternary)
    substrate_a_perp = a_ternary * (1 - (substrate_D001 * ((substrate_a_parallel/a_ternary ) - 1)))
    epilayer_a_perp = a_ternary * (1 - (epilayer_D001 * ((epilayer_a_parallel/a_ternary ) - 1)))
    substrate_epsilon_perp = (substrate_a_perp/a_ternary) - 1
    epilayer_epsilon_perp = (epilayer_a_perp/a_ternary) - 1
    substrate_epsilon_parallel = (substrate_a_parallel/a_ternary) - 1
    epilayer_epsilon_parallel = (epilayer_a_parallel/a_ternary) - 1
    substrate_hydro_average = substrate_row["aᵥ (eV)"] * (2*substrate_epsilon_parallel + substrate_epsilon_perp)
    epilayer_hydro_average = av_ternary * (2*epilayer_epsilon_parallel + epilayer_epsilon_perp)
    substrate_hydro_conduction = substrate_row["aₑ(Γ) (eV)"] * (2*substrate_epsilon_parallel + substrate_epsilon_perp)
    epilayer_hydro_conduction = ac_ternary * (2*epilayer_epsilon_parallel + epilayer_epsilon_perp)
    substrate_strainshift = 2*substrate_row["b (eV)"]*(substrate_epsilon_perp-substrate_epsilon_parallel)
    epilayer_strainshift = 2*b_ternary*(epilayer_epsilon_perp-epilayer_epsilon_parallel)
    substrate_hh = -0.5 * substrate_strainshift
    epilayer_hh = -0.5 * epilayer_strainshift
    substrate_lh = -0.5 * substrate_row["Δ₀ (eV)"] + 0.25 * substrate_strainshift + 0.5\
                   * ((substrate_row["Δ₀ (eV)"])**2+substrate_row["Δ₀ (eV)"]*substrate_strainshift+2.25*(substrate_strainshift)**2)**0.5
    epilayer_lh = -0.5 * delta0_ternary + 0.25 * epilayer_strainshift + 0.5 * ((delta0_ternary)**2+delta0_ternary*epilayer_strainshift+2.25*(epilayer_strainshift)**2)**0.5
    substrate_so = -0.5 * substrate_row["Δ₀ (eV)"] + 0.25 * substrate_strainshift - 0.5\
                   * ((substrate_row["Δ₀ (eV)"])**2+substrate_row["Δ₀ (eV)"]*substrate_strainshift+2.25*(substrate_strainshift)**2)**0.5
    epilayer_so = -0.5 * delta0_ternary + 0.25 * epilayer_strainshift - 0.5 * ((delta0_ternary)**2+delta0_ternary*epilayer_strainshift+2.25*(epilayer_strainshift)**2)**0.5

    max_value = [max(a, b) for a, b in zip(substrate_hh, substrate_lh)] # Calculates the maximum between the heavy hole and light hole
    substrate_valance = substrate_row["Eᵥ,ₐᵥ (eV)"] + (substrate_row["Δ₀ (eV)"]/3) + substrate_hydro_average + max_value
    max_value1 = [max(a, b) for a, b in zip(epilayer_hh, epilayer_lh)]# Calculates the maximum between the heavy hole and light hole
    epilayer_valance = Evav_ternary + (delta0_ternary/3) + epilayer_hydro_average + max_value1
    
    substrate_conduction = substrate_row["Eᵥ,ₐᵥ (eV)"] + (substrate_row["Δ₀ (eV)"]/3) + substrate_row["Eg(Γ) (eV)"] + substrate_hydro_conduction
    epilayer_conduction = Evav_ternary + (delta0_ternary/3) + Eg_ternary + epilayer_hydro_conduction #Conduction calculations

    substrate_bandgap = substrate_conduction - substrate_valance #Bandgap calculations
    epilayer_bandgap = epilayer_conduction - epilayer_valance

    plt.scatter(x_values, epilayer_bandgap, color='blue', label="Bandgap vs x", marker='o')#Plot info
    plt.xlabel("x")
    plt.ylabel("Bandgap (eV)")
    plt.title("Bandgap vs x")
    plt.legend()
    plt.show()    #Shows the plot

    

    
    
    


    








    
