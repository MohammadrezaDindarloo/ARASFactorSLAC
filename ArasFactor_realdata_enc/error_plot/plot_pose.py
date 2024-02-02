import matplotlib.pyplot as plt
import numpy as np
import csv
import os

direc = os.getcwd()

GTPOSE = []
with open(direc + "/error_plot/All_GT_Position.csv", 'r') as file:
    csvreader = csv.reader(file)
    header = next(csvreader)
    GTPOSE.append(header)
    for row in csvreader:
        GTPOSE.append(row)

OPTIMIZATIONPOSE = []
with open(direc + "/error_plot/OptimizedPose.csv", 'r') as file:
    csvreader = csv.reader(file)
    header = next(csvreader)
    OPTIMIZATIONPOSE.append(header)
    for row in csvreader:
        OPTIMIZATIONPOSE.append(row)

pose_x_GT = [float(row[0]) for row in GTPOSE]
pose_y_GT = [float(row[1]) for row in GTPOSE]
pose_z_GT = [float(row[2]) for row in GTPOSE]

pose_x_OP = [float(row[0]) for row in OPTIMIZATIONPOSE]
pose_y_OP = [float(row[1]) for row in OPTIMIZATIONPOSE]
pose_z_OP = [float(row[2]) for row in OPTIMIZATIONPOSE]

tr_errors_x =[]
tr_errors_y =[]
tr_errors_z =[]
for i in range(len(GTPOSE)):
    tr_errors_x.append((float(GTPOSE[i][0]) - float(OPTIMIZATIONPOSE[i][0]))*1000)
    tr_errors_y.append((float(GTPOSE[i][1]) - float(OPTIMIZATIONPOSE[i][1]))*1000)
    tr_errors_z.append((float(GTPOSE[i][2]) - float(OPTIMIZATIONPOSE[i][2]))*1000)

fig = plt.figure(figsize=(8,10))
ax = fig.subplots(6,1, sharex = True)
ax[0].plot(pose_x_OP,'r*')
ax[1].plot(pose_y_OP,'r*')
ax[2].plot(pose_z_OP,'r*')
ax[0].plot(pose_x_GT,'b')
ax[1].plot(pose_y_GT,'b')
ax[2].plot(pose_z_GT,'b')
ax[0].set_title('FG vs GT')
ax[0].legend(['FG', 'GT'])

ax[3].plot(tr_errors_x,'r*')
ax[4].plot(tr_errors_y,'r*')
ax[5].plot(tr_errors_z,'r*')
ax[3].axhline(np.mean(tr_errors_x), color='black', linestyle='--', label='Mean')
ax[4].axhline(np.mean(tr_errors_y), color='black', linestyle='--', label='Mean')
ax[5].axhline(np.mean(tr_errors_z), color='black', linestyle='--', label='Mean')
ax[3].set_title('FG vs GT (mm)')
[ax[i].grid(True) for i in range(6)]
plt.tight_layout()

plt.show()