#https://chatgpt.com/share/675eea89-a650-8001-b6b7-bb59d3d7e2f0
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
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
#boundary value array for upper and lower part of the airfoil
x_boundaries_u = [0.0]
x_boundaries_l = [0.0]

x_integration_len_u = []
x_integration_len_l = []

for i in range(len(probe_positions_u)):
    if(i == len(probe_positions_u) - 1):
        x_boundaries_u.append(100.0)
    else:
        x_boundaries_u.append((probe_positions_u[i] + probe_positions_u[i+1])/2) #boundary x coordinates for using each probe (for discrete sum integration)
    x_integration_len_u.append(x_boundaries_u[i+1] - x_boundaries_u[i]) #lengths over which the integration happens
    
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
#print("Columns available:", list(columns.keys()))

#get the normal force coefficient
def getCn(AOA):
    # Ensure the AOA column exists
    if 'Alpha' not in columns:
        print("Column 'Alpha' not found in the data.")
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
            probe_data_u = probe_data[:25]
            probe_data_l = probe_data[25:]
    
    V_inf = mu * Re / (rho * c)
            
    C_pl = probe_data_l * (1/(0.5 * rho * V_inf**2))
    C_pu = probe_data_u * (1/(0.5 * rho * V_inf**2))

    plt.figure()
    plt.plot(probe_positions_u, C_pu, label='C_p upper part', marker='o')
    plt.plot(probe_positions_l, C_pl, label='C_p lower part', marker='o')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    
    print(C_pu)
    print(C_pl)
    #calculate the normal force coefficient

# Example usage:
probe_data = getCn(5.0)  # Replace with the desired AOA