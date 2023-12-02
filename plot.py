import matplotlib.pyplot as plt
import numpy as np


file = open('3denergy.csv', 'r')

lines = file.readlines()

data = []
for line in lines:
	line_formatted = line.split(",")
	line_complete = []
	for i in line_formatted:
		line_complete.append(float(i))
	data.append(line_complete)

data = np.array(data)
fig = plt.figure()

ax = fig.add_subplot(projection='3d')
ax.scatter(data[:,0],data[:,1],data[:,2])
plt.xlim(-0.1,0.1)
plt.ylim(-0.1,0.1)
#plt.zlim(-0.2,0.2)
plt.savefig("btp.jpg", dpi=500)

plt.show()
