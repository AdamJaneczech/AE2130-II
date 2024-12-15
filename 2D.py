#https://chatgpt.com/share/675eea89-a650-8001-b6b7-bb59d3d7e2f0
import numpy as np

# File path
file_path = 'raw_2D_code.txt'

# Read the file and extract headers and data
with open(file_path, 'r') as f:
    lines = f.readlines()

# Extract headers (second line contains the headers)
headers = lines[1].strip().split()
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
column_name = 'degrees'  # Replace with desired column name
row_index = 2  # Replace with desired row index (0-based)

if column_name in columns:
    column = columns[column_name]
    value = column[row_index]
    print(column)
    print(f"Value at column '{column_name}', row {row_index}: {value}")
else:
    print(f"Column '{column_name}' not found!")