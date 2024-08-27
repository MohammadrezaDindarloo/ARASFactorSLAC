# SCAMPI Calibration Modified

Welcome to the **SCAMPI Calibration Modified** repository! This project contains the implementation of a self-calibration framework for Cable-Driven Parallel Robots (CDPRs) with sagging cables. The framework is designed to iteratively refine kinematic parameters, including anchor point locations and initial cable lengths, using a graph-based approach that considers the sagging effect of cables.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Repository Structure](#repository-structure)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Overview

This repository provides the code and resources needed to implement and test the self-calibration procedure for large-scale and small-scale CDPRs, taking into account cable sagging effects. The framework is validated against Finite Element (FE) simulations to ensure accuracy and applicability across different robot configurations.

The project is based on the paper titled *"A Graph-Based Self-Calibration Technique for Cable-Driven Robots with Sagging Cable."* The proposed method uses a factor graph that integrates onboard sensor data and a catenary cable model to improve the calibration process for CDPRs.

## Features

- **Graph-Based Calibration:** Utilizes factor graphs for efficient calibration of kinematic parameters.
- **Cable Sag Modeling:** Incorporates a catenary cable model to address sagging effects.
- **Scalability:** Applicable to both small-scale and large-scale CDPRs.
- **Validation:** Validated through FE simulations in RecurDyn software.

## Installation

To get started with the SCAMPI Calibration Modified project, follow these steps:

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/MohammadrezaDindarloo/scampi_calibration_modified.git
   cd scampi_calibration_modified
