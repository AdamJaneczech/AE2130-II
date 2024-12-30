#https://chatgpt.com/share/676e7c21-9ca8-8001-afc3-92def3482e39
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
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

probe_positions_total = np.array([
    0, 12, 21, 27, 33, 39, 45, 51, 57, 63, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 156, 162, 168, 174, 180, 186, 195, 207,  219
])

probe_positions_static = np.array([
    43.5, 55.5, 67.5, 79.5, 91.5, 103.5, 115.5, 127.5, 139.5, 151.5, 163.5, 175.5
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
def getCp(AOA, plot = False):
    # Ensure the AOA column exists
    if 'Alpha' not in columns:
        print("Column 'Alpha' not found in the data.")
        return None
    
    aoa_column = np.array(columns['Alpha'], dtype=float)
    for row_index in range(31):
        if np.isclose(aoa_column[row_index], AOA, atol=1e-6):
            probe_data = np.array([
                float(columns[f'P{str(i).zfill(3)}'][row_index]) for i in range(1, 50)
            ])
            rho = float(columns['rho'][row_index])
            probe_data_u = probe_data[:25]
            probe_data_l = probe_data[25:]
            break
    else:
        print(f"No data found for AOA = {AOA} degrees.")
        return None

    V_inf = mu * Re / (rho * c)
    C_pu = probe_data_u / (0.5 * rho * V_inf**2)
    C_pl = probe_data_l / (0.5 * rho * V_inf**2)

    # Plot Cp
    if(plot):
        plt.figure(figsize=(10, 6))
        plt.plot(probe_positions_u, C_pu, label='C_p (upper)', marker='o', color='#3DA5D9')
        plt.plot(probe_positions_l, C_pl, label='C_p (lower)', marker='o', color='#D7263D')
        plt.xlabel('x/c')
        plt.ylabel('C_p')
        plt.gca().invert_yaxis()
        plt.title(f'Pressure Coefficient at AOA = {AOA}°')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return C_pu, C_pl

def getCn(AOA):
    if(getCp(AOA) == None):
        return None
    else:
        C_pu, C_pl = getCp(AOA)
    #integrate pressure coefficients on both sides of the airfoil
    integral_upper = np.trapz(C_pu, probe_positions_u)
    integral_lower = np.trapz(C_pl, probe_positions_l)
    #compute the normal force coefficient
    C_n = (integral_lower - integral_upper)
    print(f'Integral Upper: {integral_upper}, Integral lower: {integral_lower}, C_n: {C_n}, AOA: {AOA}')

    return C_n

def getCm(AOA):
    if(getCp(AOA) == None):
        return None
    else:
        C_pu, C_pl = getCp(AOA)
    #integrate pressure coefficients multiplied by the chord-wise distance on both sides of the airfoil
    integral_upper = np.trapz(C_pu * probe_positions_u, probe_positions_u)
    integral_lower = np.trapz(C_pl * probe_positions_l, probe_positions_l)
    #compute C_m
    C_m = integral_lower - integral_upper
    print(f'Integral Upper: {integral_upper}, Integral lower: {integral_lower}, C_m: {C_m}, AOA: {AOA}')

    return C_m

def getVelocityProfile(AOA, plot = False):
    """
    Computes the velocity profile at a given angle of attack (AOA), interpolating
    static pressures only within the boundaries of the static pressure probes.

    Parameters:
        AOA (float): Angle of attack in degrees.
    """
    if 'Alpha' not in columns:
        print("Column 'Alpha' not found in the data.")
        return None

    aoa_column = np.array(columns['Alpha'], dtype=float)
    for row_index in range(31):
        if np.isclose(aoa_column[row_index], AOA, atol=1e-6):
            total_pressures = np.array([
                float(columns[f'P{str(i).zfill(3)}'][row_index]) for i in range(50, 97)
            ])
            static_pressures = np.array([
                float(columns[f'P{str(i).zfill(3)}'][row_index]) for i in range(98, 110)
            ])
            rho = float(columns['rho'][row_index])
            break
    else:
        print(f"No data found for AOA = {AOA} degrees.")
        return None

    # Interpolate static pressures strictly within static probe boundaries
    static_pressure_interp = interp1d(
        probe_positions_static, static_pressures, kind='linear',
        bounds_error=False, fill_value="extrapolate"
    )

    valid_indices = (probe_positions_total >= probe_positions_static[0]) & (probe_positions_total <= probe_positions_static[-1])
    restricted_total_pressures = total_pressures[valid_indices]
    restricted_positions = probe_positions_total[valid_indices]
    restricted_interpolated_static_pressures = static_pressure_interp(restricted_positions)

    # Compute velocity using the restricted range
    velocities = np.sqrt(2 * (restricted_total_pressures - restricted_interpolated_static_pressures) / rho)
    V_inf = mu * Re / (rho * c)
    velocity_deficit = V_inf - velocities

    # Plot velocity profile
    if(plot):
        plt.figure(figsize=(10, 6))
        plt.plot(restricted_positions, velocities, label='Velocity (m/s)', marker='o', color='#3DA5D9')
        plt.axhline(V_inf, color='#50514F', linestyle='--', label='Free-stream Velocity')
        plt.xlabel('Transverse Axis (y)')
        plt.ylabel('Velocity (m/s)')
        plt.title(f'Velocity Profile at AOA = {AOA}° (Within Static Probe Range)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return velocities, velocity_deficit

def getCt(AOA):
    """
    Compute the tangential force coefficient (C_t) for a given angle of attack (AOA).
    """
    if getCp(AOA) is None:
        return None
    else:
        C_pu, C_pl = getCp(AOA)
    
    # Load airfoil coordinates
    airfoil_coords = np.loadtxt('SD6060-104-88_180.dat')

    # Upper surface: Read first 90 lines and reverse
    upper_surface = airfoil_coords[:90][::-1]

    # Lower surface: Read from line 91 till the end
    lower_surface = airfoil_coords[90:]

    # Separate x and y coordinates for each surface
    x_upper, y_upper = upper_surface[:, 0], upper_surface[:, 1]
    x_lower, y_lower = lower_surface[:, 0], lower_surface[:, 1]

    # Compute slopes (dy/dx) for upper and lower surfaces
    dy_dx_upper = np.gradient(y_upper, x_upper)
    dy_dx_lower = np.gradient(y_lower, x_lower)

    # Compute local angles (theta) for upper and lower surfaces
    theta_upper = np.arctan(dy_dx_upper)
    theta_lower = np.arctan(dy_dx_lower)

    # Interpolate theta arrays to match probe positions
    interp_theta_upper = np.interp(probe_positions_u, x_upper, theta_upper)
    interp_theta_lower = np.interp(probe_positions_l, x_lower, theta_lower)

    # Compute tangential force coefficient contributions
    delta_x_upper = np.diff(probe_positions_u)
    delta_x_lower = np.diff(probe_positions_l)

    Ct_upper = np.sum(C_pu[:-1] * np.cos(interp_theta_upper[:-1]) * delta_x_upper)
    Ct_lower = np.sum(C_pl[:-1] * np.cos(interp_theta_lower[:-1]) * delta_x_lower)

    # Total tangential force coefficient
    C_t = Ct_upper + Ct_lower

    print(f"Tangential Force Coefficient (C_t): {C_t:.4f}, AOA: {AOA}")
    return C_t

# Initialize arrays to store AOA, C_n, and C_m values
aoa_range = range(-15, 15)  # Define AOA range
cn_values = []
cm_values = []
ct_values = []
valid_aoa = []

for aoa in aoa_range:
    # Calculate C_n and C_m only if data is available
    try:
        # Get C_n and C_m using respective functions
        C_n = getCn(aoa)
        C_m = getCm(aoa)
        C_t = getCt(aoa)
        cn_values.append(C_n)
        cm_values.append(C_m)
        ct_values.append(C_t)

        # Store valid AOA
        valid_aoa.append(aoa)
    except:
        # Skip angles of attack without data
        continue

# Filter out None values from cn_values and cm_values along with their corresponding valid_aoa
filtered_cn_values = [cn for cn in cn_values if cn is not None]
filtered_cm_values = [cm for cm in cm_values if cm is not None]
filtered_ct_values = [ct for ct in ct_values if ct is not None]
filtered_valid_aoa = [aoa for idx, aoa in enumerate(valid_aoa) if cn_values[idx] is not None]

# Plot C_n vs AOA
plt.figure(figsize=(10, 6))
plt.plot(filtered_valid_aoa, filtered_cn_values, label="C_n", marker="o", color='#3DA5D9')
plt.xlabel("Angle of Attack (°)")
plt.ylabel("C_n")
plt.title("Normal Force Coefficient (C_n) vs Angle of Attack")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Plot C_m vs AOA
plt.figure(figsize=(10, 6))
plt.plot(filtered_valid_aoa, filtered_cm_values, label="C_m", marker="o", color='#3DA5D9')
plt.xlabel("Angle of Attack (°)")
plt.ylabel("C_m")
plt.title("Moment Coefficient (C_m) vs Angle of Attack")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Plot C_t vs AOA
plt.figure(figsize=(10, 6))
plt.plot(filtered_valid_aoa, filtered_ct_values, label="C_t", marker="o", color='#3DA5D9')
plt.xlabel("Angle of Attack (°)")
plt.ylabel("C_t")
plt.title("Tangential Force Coefficient (C_t) vs Angle of Attack")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()