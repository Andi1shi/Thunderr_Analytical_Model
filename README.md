# Thunderstorm Simulation and Optimization

## Description
This MATLAB script is designed for simulating horizontal slowly-varying mean wind velocity caused by traveling downbursts. It combines three velocity fields: radial impinging jet, downburst translational, and ambient ABL (Atmospheric Boundary Layer) wind velocities. The model employs a parametric-analytical approach for spatiotemporal reconstruction of downburst events and optimizes downburst geometric and kinematics parameters against real measurements using the Teaching-Learning-Based Optimization (TLBO) algorithm.

## Features
- **Parametric-Analytical Model**: Uses a well-defined analytical model to simulate downburst wind fields.
- **Optimization**: Employs TLBO to optimize the geometric and kinematic parameters of downbursts based on recorded data.
- **Visualization**: Includes functionality to plot and analyze the simulation results and optimization process.

## Installation
No specific installation is required other than having MATLAB installed on your system. This script was developed and tested in MATLAB R2021a.

## Usage
- **Load Data**: Ensure that the thunderstorm data file `Thunderstorm_Romania_MA.mat` is in the MATLAB path or in the same directory as the script.
- **Run the Script**: Open the script in MATLAB and press Run. The script clears the current MATLAB environment, so ensure you save your work before running the script.
- **Adjust Parameters**: Users can modify simulation time settings and station location parameters directly in the script under the "Simulation Parameters" section.

## Configuration
Before running the simulation, ensure that you configure the simulation parameters according to your specific requirements. The parameters include the simulation time settings, station locations, and bounds for various simulation parameters like the downburst radius and velocity.

## Dependencies
- MATLAB (tested on R2021a)
- No additional MATLAB toolboxes are required.

## Citation
If you use this code in your research or wish to refer to it in your academic publication, please cite it as follows:


Replace `https://github.com/username/repository` with the actual URL of your GitHub repository. Be sure to update this citation with a DOI if you archive the repository on Zenodo or a similar service.

## License
MIT License

Copyright (c) 2024, University of Genoa, Department of Civil, Chemical, and Environmental Engineering.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this MATLAB script, including without limitation the rights
to use, copy, modify, merge and publish. 


## Acknowledgments
- Xhelaj, A., Burlando, M., Solari, G. (2020). A general-purpose analytical model for reconstructing the thunderstorm outflows of travelling downbursts immersed in ABL flows. Journal of Wind Engineering and Industrial Aerodynamics, 207, 104373.
- Xhelaj, A., Burlando, M.: Application of metaheuristic optimization algorithms to evaluate the geometric and kinematic parameters of downbursts. Advances in Engineering Software. Volume 173, November 2022, 103203.

