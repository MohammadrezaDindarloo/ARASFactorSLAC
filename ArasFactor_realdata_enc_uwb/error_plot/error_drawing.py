import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import csv
import os

direc = os.getcwd()


rows = []
with open(direc + "/error_plot/pos_i_cpp_test.csv", 'r') as file:
    csvreader = csv.reader(file)
    header = next(csvreader)
    rows.append(header)
    for row in csvreader:
        rows.append(row)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Optimizer
ax.scatter(-2.00908 , -8.29336, 8.52047, marker='o', c='r')
ax.scatter(2.52119, -8.39907, 8.4358, marker='o', c='r')
ax.scatter(2.7661, 4.78913, 8.31409, marker='o', c='r')
ax.scatter(-1.71295, 4.90232, 8.30433, marker='o', c='r')

# GT
ax.scatter(-1.9874742 , -8.31965637, 8.47184658, marker='^', c='b')
ax.scatter(2.52022147, -8.38887501, 8.46931362, marker='^', c='b')
ax.scatter(2.71799795, 4.77520639, 8.36416322, marker='^', c='b')
ax.scatter(-1.79662371, 4.83333111, 8.37001991, marker='^', c='b')

pose_x = [float(row[0]) for row in rows]
pose_y = [float(row[1]) for row in rows]
pose_z = [float(row[2]) for row in rows]

ax.scatter(pose_x, pose_y, pose_z, label='Points', c='black',  s=10)

ax.plot(pose_x, pose_y, pose_z, label='Connected Lines', linestyle='-' , c='r')

# Change the color of the first point to green
ax.scatter(pose_x[0], pose_y[0], pose_z[0], label='First Point', c='green', s=30)

# Draw a sphere around the specified point
center = np.array([-1.9874742 , -8.31965637, 8.47184658])
radius = 0.7  # Adjust the radius as needed
phi, theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
x_sphere = center[0] + radius * np.sin(phi) * np.cos(theta)
y_sphere = center[1] + radius * np.sin(phi) * np.sin(theta)
z_sphere = center[2] + radius * np.cos(phi)
ax.plot_surface(x_sphere, y_sphere, z_sphere, color='cyan', alpha=0.2, label='Sphere')

center = np.array([2.52022147, -8.38887501, 8.46931362])
radius = 0.7  # Adjust the radius as needed
phi, theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
x_sphere = center[0] + radius * np.sin(phi) * np.cos(theta)
y_sphere = center[1] + radius * np.sin(phi) * np.sin(theta)
z_sphere = center[2] + radius * np.cos(phi)
ax.plot_surface(x_sphere, y_sphere, z_sphere, color='cyan', alpha=0.2, label='Sphere')

center = np.array([2.71799795, 4.77520639, 8.36416322])
radius = 0.7  # Adjust the radius as needed
phi, theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
x_sphere = center[0] + radius * np.sin(phi) * np.cos(theta)
y_sphere = center[1] + radius * np.sin(phi) * np.sin(theta)
z_sphere = center[2] + radius * np.cos(phi)
ax.plot_surface(x_sphere, y_sphere, z_sphere, color='cyan', alpha=0.2, label='Sphere')

center = np.array([-1.79662371, 4.83333111, 8.37001991])
radius = 0.7  # Adjust the radius as needed
phi, theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
x_sphere = center[0] + radius * np.sin(phi) * np.cos(theta)
y_sphere = center[1] + radius * np.sin(phi) * np.sin(theta)
z_sphere = center[2] + radius * np.cos(phi)
ax.plot_surface(x_sphere, y_sphere, z_sphere, color='cyan', alpha=0.2, label='Sphere')

ax.quiver(0, 0, 0, 1, 0, 0, color='blue', label='X-axis')
ax.quiver(0, 0, 0, 0, 1, 0, color='green', label='Y-axis')
ax.quiver(0, 0, 0, 0, 0, 1, color='purple', label='Z-axis')


# Define the vertices of the rectangular box
x_min, y_min, z_min = np.min(pose_x) -0.2, np.min(pose_y) -0.5, np.min(pose_z) -0.1
x_max, y_max, z_max = np.max(pose_x) +0.2, np.max(pose_y) +0.5, np.max(pose_z) +0.5

vertices = [
    (x_min, y_min, z_min),
    (x_min, y_max, z_min),
    (x_max, y_max, z_min),
    (x_max, y_min, z_min),
    (x_min, y_min, z_max),
    (x_min, y_max, z_max),
    (x_max, y_max, z_max),
    (x_max, y_min, z_max)
]

# Define the rectangular box faces
faces = [
    [vertices[0], vertices[1], vertices[2], vertices[3]],
    [vertices[4], vertices[5], vertices[6], vertices[7]],
    [vertices[0], vertices[3], vertices[7], vertices[4]],
    [vertices[1], vertices[2], vertices[6], vertices[5]],
    [vertices[0], vertices[1], vertices[5], vertices[4]],
    [vertices[2], vertices[3], vertices[7], vertices[6]]
]

# Plot the rectangular box
ax.add_collection3d(Poly3DCollection(faces, facecolors=[(135/255, 206/255, 235/255)]*6, linewidths=0.7, edgecolors='black', alpha=0.1))




ax.set(xlim=(-2.8, 3.5), ylim=(-8.5, 5), zlim=(0, 9.5),
       xlabel='X', ylabel='Y', zlabel='Z')

ax.view_init(elev=23, azim=-21)
plt.show()


