import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

input_file_path = 'clean_2d_trajectory.dat'

# Open the input file in read mode
with open(input_file_path, 'r') as input_file:
    # Read all lines from the file
    lines = input_file.readlines()

# Extract numerical values and create a list of lists
data_list = [list(map(float, line.split())) for line in lines]

# Convert the list of lists to a 2D NumPy array
data1 = np.array(data_list)
data2 = np.copy(data1)
data3 = np.copy(data1)
data4 = np.copy(data1)

data2[:, 1] = -data1[:, 1]
data3[:, 0] = -data3[:, 0]
data3[:, 1] = -data3[:, 1]
data4[:, 0] = -data1[:, 0]

data2 = data2[::-1, :]
data4 = data4[::-1, :]

data = np.concatenate((data1, data2, data3, data4), 0)
#data = data[::74]
data_n = np.array([
    [0.0,    0.0492],
    #[0.02432,0.04377],
    [0.03363,0.0365],
    [0.0504, 0.0],
    [0.03363,-0.0365],
    #[0.02432,-0.04377],
    [0.0,-0.0492],
    #[-0.02432,-0.04377],
    [-0.03363,-0.0365],
    [-0.0504,0.0],
    [-0.03363,0.0365],
    #[-0.02432,0.04377],
    [0.0, 0.0492]
        ])

x = data_n[:, 0]
y = data_n[:, 1]
ax = plt.figure()
# Create a spline interpolation
param = np.linspace(0, 1, x.size)
spl = make_interp_spline(param, np.c_[x,y], k=2)
xnew, y_smooth = spl(np.linspace(0, 1, x.size * 100)).T
#plt.scatter(data[:,0], data[:,1])
plt.plot(xnew, y_smooth, color='black')
plt.savefig('trajectory.png', dpi=300)
plt.xlim(-0.05,0.05)
plt.ylim(-0.05,0.05)
plt.legend()
plt.show()
