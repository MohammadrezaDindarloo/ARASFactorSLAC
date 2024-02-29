## scampi_factor_graph_optimization
For generating the results and error a jupyter notebook script is provided. 
1. At first for saving the FK, IK and L_measure outputs we need to make a senario and feed some initial pos for master code to save the results. 
2. The real data contain 1000 frame and all the data is saved in `matlab_data`, from these data the initial poses for all these frame saved as csv file and are located in `raw_python_data_fedding_cpp_optimization`. 
3. Next with compiling and running `master.cpp`, using these initial poses, the inverse and forward results for fator graph optimization will be saved in `data_from_cpp_to_python`.
4. At the last by running the notebook the results will be generated in `Results/real_data`.
