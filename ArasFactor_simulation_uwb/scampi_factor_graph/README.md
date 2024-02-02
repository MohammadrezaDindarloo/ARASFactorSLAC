## scampi_factor_graph_optimization
This file contain the Cpp code for solver in factor graph optimizer. 
There are two factor here developed:
- Inverse kinematic factor which gives the position and initial rotation of robot and find the optimize value of cabel forces, cables length, catenary variables and rotation of robot.
- Forward kinematic factor which gives the cales length and initial rotation of robot and find the optimize value of cabel forces, robot position, catenary variables and rotation of robot.
### include files
The necessaray header files for optimization should add to `lib/include` file. This FK and IK are the residuals and costs for factors. These includes will be generaing in `cpp/symforce/sym` in this directory and are needed to paste to `lib/include` without changing in any name. 
### Run
Running this code just needed to modify your initial pose in `master.c` and compile it.

