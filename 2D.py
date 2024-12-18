#https://chatgpt.com/share/675eea89-a650-8001-b6b7-bb59d3d7e2f0
import numpy as np
import matplotlib.pyplot as plt
# File path
file_path = 'raw_2D_code.txt'

Re = 1.8e5
mu = 18.23e-6 #Pa*s at T = 22.2 C
c = 0.16

# Constant probe locations (chord positions, % chord length)
probe_positions_u = np.array([
    0, 0.35626, 1.33331, 3.66108, 7.2922, 11.35604, 15.59135, 19.91328, 
    24.28443, 28.68627, 33.10518, 37.53128, 41.95991, 46.38793, 50.8156, 
    55.2486, 59.69223, 64.13685, 68.579, 73.02401, 77.47357, 81.93114, 
    86.38589, 90.8108, 100
])

probe_positions_l = np.array([
    0, 0.43123, 1.47147, 3.92479, 7.79506, 
    12.0143, 16.32276, 20.67013, 25.03792, 29.41554, 33.79772, 38.18675, 
    42.57527, 46.96278, 51.35062, 55.73662, 60.12075, 64.50502, 68.8901, 
    73.28011, 77.67783, 82.07965, 86.47978, 100
])

x_boundaries_u = [0.0]
x_integration_len_u = []
x_boundaries_l = [0.0]
x_integration_len_l = []

for i in range(len(probe_positions_u)):
    if(i == len(probe_positions_u) - 1):
        x_boundaries_u.append(100.0)
    else:
        x_boundaries_u.append((probe_positions_u[i] + probe_positions_u[i+1])/2)
    x_integration_len_u.append(x_boundaries_u[i+1] - x_boundaries_u[i])
    
for i in range(len(probe_positions_l)):
    if(i == len(probe_positions_l) - 1):
        x_boundaries_l.append(100.0)
    else:
        x_boundaries_l.append((probe_positions_l[i] + probe_positions_l[i+1])/2)
    x_integration_len_l.append(x_boundaries_l[i+1] - x_boundaries_l[i])

print(x_boundaries_l)

# Read the file and extract headers and data
with open(file_path, 'r') as f:
    lines = f.readlines()

# Extract headers (first line contains the headers)
headers = lines[0].strip().split()
units = lines[1].strip().split()
data = []

# Skip the first two lines (headers and units) and read the data
for line in lines[2:]:
    if line.strip():  # Skip empty lines
        # Convert each value to float if possible, else leave it as a string
        row = []
        for value in line.strip().split():
            try:
                row.append(float(value))
            except ValueError:
                row.append(value)
        data.append(row)

# Convert data to a NumPy array
data_array = np.array(data, dtype=object)

# Create a dictionary of NumPy arrays for each column
columns = {header: data_array[:, idx] for idx, header in enumerate(headers)}

# Access confirmation
print("Columns available:", list(columns.keys()))

# Example: Accessing specific column and specific value
column_name = 'P001'  # Replace with desired column name
row_index = 2  # Replace with desired row index (0-based)

if column_name in columns:
    column = columns[column_name]
    value = column[row_index]
    print(column)
    print(f"Value at column '{column_name}', row {row_index}: {value}")
else:
    print(f"Column '{column_name}' not found!")

#get the normal force coefficient without
def getCn(AOA):
    """
    Extracts the row corresponding to the specified AOA and organizes probe data into an array.
    
    :param AOA: Angle of Attack in degrees (numeric)
    """
    # Ensure the AOA column exists
    if 'Alpha' not in columns:
        print("Column 'AOA' not found in the data.")
        return None
    
    # Convert AOA column to numeric values
    aoa_column = np.array(columns['Alpha'], dtype=float)
    rho = 0
    # Find the row corresponding to the specified AOA
    for row_index in range(31):  # Only search within rows 0 to 31
        if np.isclose(aoa_column[row_index], AOA, atol=1e-6):  # Compare with tolerance for floats
            # Extract probe data for the given AOA
            probe_data = np.array([
                float(columns[f'P{str(i).zfill(3)}'][row_index]) for i in range(1, 50)
            ])  # Adjust range based on probe columns
            print(f"Probe data for AOA = {AOA}:\n{probe_data}")
            rho = float(columns[f'rho'][row_index])
    
    V_inf = mu * Re / (rho * c)
    C_pl = 0
    C_pu = 0
    print(len(x_integration_len_l))
    for i in range(0,25):
        C_pu += probe_data[i] * x_integration_len_u[i]
        print(i)
    for i in range(0,24):
        C_pl += probe_data[i+25] * x_integration_len_l[i]
        print(i)
    print(C_pu)
    print(C_pl)        
    C_pl = C_pl * (1/(0.5 * rho * V_inf**2))
    C_pu = C_pu * (1/(0.5 * rho * V_inf**2))
    #calculate the normal force coefficient

# Example usage:
probe_data = getCn(5.0)  # Replace with the desired AOA
"""
def plot_pressure_distribution(AOA):

    #Plots the pressure distribution (upper and lower surfaces) for a given AOA.
    #:param AOA: Angle of Attack in degrees
    
    if AOA not in pressure_data:
        print(f"No data available for AOA = {AOA}")
        return

    # Extract data for given AOA
    pressures = np.array(pressure_data[AOA])

    # Split into upper and lower surfaces
    upper_surface_x = probe_positions[:25]
    upper_surface_y = pressures[:25]

    lower_surface_x = probe_positions[25:]
    lower_surface_y = pressures[25:]

    # Plot the pressure distribution
    plt.figure(figsize=(10, 6))
    plt.plot(upper_surface_x, upper_surface_y, 'o-', label="Upper Surface")
    plt.plot(lower_surface_x, lower_surface_y, 'o-', label="Lower Surface")
    plt.gca().invert_yaxis()  # Pressure decreases downward

    plt.title(f"Pressure Distribution at AOA = {AOA}°")
    plt.xlabel("Chord Position (%)")
    plt.ylabel("Static Pressure Difference")
    plt.legend()
    plt.grid()
    plt.show()

# Example usage: Plot for AOA = 5°
plot_pressure_distribution(5)
"""