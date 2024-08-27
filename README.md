# ARASFactorSLAC

[![IROS 2023](https://img.shields.io/badge/Conference-IROS%202023-orange.svg)](https://www.iros2023.org/)
[![C++](https://img.shields.io/badge/Language-C++-blue.svg)](https://isocpp.org/)
[![Python](https://img.shields.io/badge/Language-Python-green.svg)](https://www.python.org/)
[![GTSAM](https://img.shields.io/badge/Library-GTSAM-brightgreen.svg)](https://gtsam.org/)
[![Symforce](https://img.shields.io/badge/Library-Symforce-yellow.svg)](https://symforce.org/)
[![RecurDyn](https://img.shields.io/badge/Tool-RecurDyn-red.svg)](https://www.functionbay.org/recurdyn/)


Welcome to the **ARASFactorSLAC** repository! This project contains the implementation of a self-calibration and localization framework for Cable-Driven Parallel Robots (CDPRs) with sagging cables. The framework is designed to iteratively refine kinematic parameters, including anchor point locations and initial cable lengths, while simultaneously localizing the end-effector using a graph-based approach that considers the sagging effect of cables.

---

## 📚 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Usage](#usage)
- [Repository Structure](#repository-structure)
- [Results](#results)
- [Citation](#citation)

---

## Overview

This repository provides the code and resources needed to implement and test the simultaneous localization and calibration (SLAC) procedure for large-scale and small-scale CDPRs, taking into account cable sagging effects. The framework is validated against Finite Element (FE) simulations to ensure accuracy and applicability across different robot configurations.

The project is based on the paper titled *"A Graph-Based Self-Calibration Technique for Cable-Driven Robots with Sagging Cable."* The proposed method uses a factor graph that integrates onboard sensor data and a catenary cable model to improve both the localization and calibration processes for CDPRs.

---

## 🚀 Features

- **Graph-Based SLAC:** Utilizes factor graphs for efficient simultaneous localization and calibration of kinematic parameters.
- **📏 Cable Sag Modeling:** Incorporates a catenary cable model to address sagging effects.
- **🔄 Scalability:** Applicable to both small-scale and large-scale CDPRs.
- **✅ Validation:** Validated through FE simulations in RecurDyn software.

---


## Usage

To run the calibration and localization procedure using ARASFactorSLAC, follow these steps:

1. **Setting Your Robot Specific Details:**
   - Before running the main calibration and localization procedure, you need to configure the system model to match your specific robot.
   - **Provide Robot Characteristics:**
     - Input your robot's details into the system model, including:
       - Robot mass
       - Robot center of mass (CoM)
       - Robot size
       - Cable attachment points on the robot
   - **Generate the Model:**
     - Use the `model_generator.py` script to create the model for your robot. This script will generate all the necessary model files, including the Jacobians, and save them in the `model_system_output` directory.
   - **Replace Model Files:**
     - After generating the model, replace the existing robot model files in the project with the newly generated ones from the `model_system_output` directory. These files should be placed in the appropriate directories as required by your scenario.

2. **Prepare the Data:**
   - Collect the required sensor data, including:
     - End-effector poses
     - Relative cable-length measurements
     - Tension values for the reference cable at the end-effector attachment point
   - Ensure your data is formatted as required by the example scripts provided.

3. **Run Calibration and Localization:**
   - Use the provided scripts to execute the ARASFactorSLAC procedure. Compile and run the `main.cpp` corresponding to your scenario.

4. **Analyze Results:**
   - The results will be saved in the specified output directory. These will include:
     - Calibrated kinematic parameters
     - Localization data for the end-effector
     - Evaluation metrics such as error rates and accuracy
   - Review the output files and consider visualizing the calibration and localization results for further analysis.

---


## Our Results from Provided Dataset

The ARASFactorSLAC framework has been rigorously validated through Finite Element (FE) simulations using the RecurDyn software. The results demonstrate the accuracy and effectiveness of the proposed simultaneous localization and calibration (SLAC) procedure for both small and large-scale Cable-Driven Parallel Robots (CDPRs).

### Model Verification

| **Metric**                                   | **Large-Scale Robot**           | **Small-Scale Robot**          |
|----------------------------------------------|---------------------------------|--------------------------------|
| **Mean Percentage Error in Cable Length (MPE-L)** | 0.0089%                          | 0.0074%                         |
| **Mean Percentage Error in Force (MPE-F)**   | 0.9154%                          | 0.9046%                         |
| **Cable Length Range**                       | 134.7 to 200.7 meters            | 10.8 to 25.9 meters             |
| **Cable Force Range**                        | 489.1 to 954.3 N                 | 12.4 to 26.2 N                  |

### Calibration Results

| **Scenario**                | **Pulley Average Error (m)** | **Offset Average Error (m)** |
|-----------------------------|-----------------------------|------------------------------|
| **Large-Scale Robot (9 Poses)**  | 0.387                       | 0.380                        |
| **Large-Scale Robot (18 Poses)** | 0.226                       | 0.217                        |
| **Large-Scale Robot (35 Poses)** | 0.191                       | 0.173                        |
| **Small-Scale Robot (8 Poses)**  | 0.149                       | 0.115                        |

### Importance of Cable Sag Consideration

The significance of considering cable sag in the calibration process was highlighted by comparing results obtained using a simplified rigid cable model versus the catenary model employed in ARASFactorSLAC. The omission of cable sag led to a significant degradation in calibration accuracy, with the mean absolute error increasing from **0.19 meters** to **2.34 meters** for the large-scale robot, underscoring the critical importance of accurate cable modeling.

---

These results validate the robustness and precision of the ARASFactorSLAC framework, making it a powerful tool for the simultaneous localization and calibration of CDPRs across varying scales and conditions.

## 📄 Citation

If you use the ARASFactorSLAC framework in your research or project, please cite the following papers:

```bibtex
@inproceedings{khorrambakht2023graph,
  title={Graph-Based Visual-Kinematic Fusion and Monte Carlo Initialization for Fast-Deployable Cable-Driven Robots},
  author={Khorrambakht, Rooholla and Damirchi, Hamed and Dindarloo, MR and Saki, A and Khalilpour, SA and Taghirad, Hamid D and Weiss, Stephan},
  booktitle={2023 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)},
  pages={1832--1839},
  year={2023},
  organization={IEEE}
}

@inproceedings{Coming IROS 2024 ...,
  title={A Graph-Based Self-Calibration Technique for Cable-Driven Robots with Sagging Cable},
  author={Dindarloo, MR and Mirjalili, AS and Khalilpour, SA and Khorrambakht, Rooholla and Weiss, Stephan and Taghirad, Hamid D},
  booktitle={2024 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)},
  year={2024},
  organization={IEEE}
}
