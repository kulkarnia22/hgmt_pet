import numpy as np

lor11_file = "data/HGMTPointVac11.lor"
lor12_file = "data/HGMTPointVac12.lor"
lor21_file = "data/HGMTPointVac21.lor"
lor22_file = "data/HGMTPointVac22.lor"

data11 = np.fromfile(lor11_file, dtype = np.float64)
data12 = np.fromfile(lor12_file, dtype = np.float64)
data21 = np.fromfile(lor21_file, dtype = np.float64)
data22 = np.fromfile(lor22_file, dtype = np.float64)

num_lors11 = data11.size/9
num_lors12 = data12.size/9
num_lors21 = data21.size/9
num_lors22 = data22.size/9

print("(1,1) number = " + str(num_lors11/1000000))
print("(1,2) number = " + str(num_lors12/1000000))
print("(2,1) number = " + str(num_lors21/1000000))
print("(2,2) number = " + str(num_lors22/1000000))