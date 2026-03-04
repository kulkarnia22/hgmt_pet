import numpy as np
import os

"""# Load binary file of doubles
binary_file = "data/det_angle_output.data"
angles = np.fromfile(binary_file, dtype=np.float64)

# Save to text file (one value per line)
text_file = "plots/det_angule_output.txt"
np.savetxt(text_file, angles, fmt="%.6f")"""

"""photon_data = np.fromfile("data/Incoming_Photons.data", dtype=np.float64)

#print(f"Converted {len(angles)} values from binary to {text_file}")
print(photon_data[:100])"""

import numpy as np

filename = "data/HGMTPointVac.lor"

# Each LOR has 9 doubles
lor_size = 9

# Read file as double array
data = np.fromfile(filename, dtype=np.float64)

# Reshape into rows of 9 doubles
lors = data.reshape(-1, lor_size)

print("First 20 LORs:\n")

for i in range(min(20, len(lors))):
    x, y, z, c_xx, c_xy, c_xz, c_yy, c_yz, c_zz = lors[i]
    
    print(f"LOR {i}:")
    print(f"  center = ({x:.4f}, {y:.4f}, {z:.4f})")
    print(f"  covariance = [{c_xx:.4f}, {c_xy:.4f}, {c_xz:.4f}, {c_yy:.4f}, {c_yz:.4f}, {c_zz:.4f}]")
    print()
