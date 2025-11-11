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

#data check
lor_file = "data/HGMTPointVac.lor"
#photon_file = "data/Incoming_Photons.data"

# file sizes
lor_bytes = os.path.getsize(lor_file)
#photon_bytes = os.path.getsize(photon_file)

# number of doubles
lor_doubles = lor_bytes // 8
#photon_doubles = photon_bytes // 8

# number of LORs
num_lors = lor_doubles // 9
# number of photons
#num_photons = photon_doubles

#print("Num photons:", num_photons)
print("Num LORs:", num_lors)
print("Sensitivity:", num_lors/1000000)
#print("Photon/2:", num_photons // 2)
#print("Match?", num_photons // 2 == num_lors)
