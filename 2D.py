#https://chatgpt.com/share/675eea89-a650-8001-b6b7-bb59d3d7e2f0
#https://chatgpt.com/share/676be685-ff4c-8001-958e-3e83bf74113e
import numpy as np
import math
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

probe_positions_u /= 100
probe_positions_l /= 100

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
def getCp(AOA):
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
            #print(f"Probe data for AOA = {AOA}:\n{probe_data}")
            rho = float(columns[f'rho'][row_index])
            probe_data_u = probe_data[:25]
            probe_data_l = probe_data[25:]
    
    V_inf = mu * Re / (rho * c)
            
    C_pl = probe_data_l * (1/(0.5 * rho * V_inf**2))
    C_pu = probe_data_u * (1/(0.5 * rho * V_inf**2))


    return C_pu, C_pl

def getForceCoeffs(AOA, C_pl, C_pu):
    integral_upper = np.trapz(C_pu, probe_positions_u)
    integral_lower = np.trapz(C_pl, probe_positions_l)
    C_n = (integral_lower - integral_upper)
    print(f'Integral Upper: {integral_upper}, Integral lower: {integral_lower}, C_n: {C_n}, AOA: {AOA}')
    integral_upper = np.trapz(C_pu * probe_positions_u, probe_positions_u)
    integral_lower = np.trapz(C_pl * probe_positions_l, probe_positions_l)
    C_m = integral_lower - integral_upper
    print(f'Integral Upper: {integral_upper}, Integral lower: {integral_lower}, C_m: {C_m}, AOA: {AOA}')
    

    # get the cl and cd using aoa & trig

# Define the AOAs to process
aoas = [-6.0, -3.0, 0.0, 3.0, 6.0, 14.0]  # Example with more AOAs

# Initialize arrays to store results
results = [getCp(aoa) for aoa in aoas]

# Calculate the grid size for subplots (3 columns per row)
columns = 3
rows = math.ceil(len(aoas) / columns)  # Number of rows needed

# Create subplots dynamically
plt.figure(figsize=(columns * 5, rows * 4))  # Adjust figure size based on rows and columns

for i, (aoa, (C_pu, C_pl)) in enumerate(zip(aoas, results), start=1):
    plt.subplot(rows, columns, i)  # Create subplot with dynamic rows and columns
    plt.plot(probe_positions_u, C_pu, label='C_p upper part', marker='o', color='#187795')
    plt.plot(probe_positions_l, C_pl, label='C_p lower part', marker='o', color='#F76F8E')
    plt.xlabel('x/c')
    plt.ylabel(f'C_p with AOA = {aoa} deg')
    plt.gca().invert_yaxis()
    if i == 1:  # Add legend to the first subplot only
        plt.legend()
    getForceCoeffs(aoa, C_pl, C_pu)

plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to avoid overlap with title
plt.show()