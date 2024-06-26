# Thunderstorm Simulation and Optimization

## Description
This MATLAB script is designed for simulating horizontal slowly-varying mean wind velocity caused by traveling downbursts. It combines three velocity fields: radial impinging jet, downburst translational, and ambient ABL (Atmospheric Boundary Layer) wind velocities. The model employs a parametric-analytical approach for spatiotemporal reconstruction of downburst events and optimizes downburst geometric and kinematics parameters against real measurements using the Teaching-Learning-Based Optimization (TLBO) algorithm.

## Features
- **Parametric-Analytical Model**: Uses an analytical model to simulate downburst wind fields.
- **Optimization**: Employs TLBO to optimize the geometric and kinematic parameters of downbursts based on recorded data.
- **Visualization**: Includes functionality to plot and analyze the simulation results and optimization process.

## Installation
No specific installation is required other than having MATLAB installed on your system. This script was developed and tested in MATLAB R2021a.

## Usage
- **Load Data**: Ensure that the thunderstorm data file `Thunderstorm_Romania_MA.mat` or a similar dataset is accessible within the MATLAB path or in the same directory as the script.
- **Run the Script**: Open the script `MultiRunV01.m` in MATLAB and execute it. Note that the script utilizes parallel computing (parfor), enhancing computational efficiency.
- **Adjust Parameters**: Users can change simulation settings, refine downburst parameters by adjusting the lower and upper bounds, and specify the number of runs for metaheuristic optimization. Conducting multiple runs ensures robust optimization, leading to the identification of the best solution.

## Configuration
Before running the simulation, ensure that you configure the simulation parameters according to your specific requirements. The parameters include the simulation time settings, station locations, and bounds for various simulation parameters like the downburst radius and velocity.

## Dependencies
- MATLAB (tested on R2021a)
- No additional MATLAB toolboxes are required.


## Citation
If you use this code in your research or wish to refer to it in your academic publication, please cite it as follows:

```bibtex
@software{thunderstorm_simulation_optimization,
  author = {Andi Xhelaj},
  title = {Thunderstorm Simulation and Optimization Script},
  year = 2024,
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/Andi1shi/Thunderr_Analytical_Model}}
}
```
## License
MIT License

Copyright (c) 2024, University of Genoa, Department of Civil, Chemical, and Environmental Engineering.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this MATLAB script, including without limitation the rights
to use, copy, modify, merge and publish. 


## References
- Xhelaj, A., Burlando, M., Solari, G. (2020). A general-purpose analytical model for reconstructing the thunderstorm outflows of travelling downbursts immersed in ABL flows. Journal of Wind Engineering and Industrial Aerodynamics, 207, 104373.
- Xhelaj, A., Burlando, M. (2022): Application of metaheuristic optimization algorithms to evaluate the geometric and kinematic parameters of downbursts. Advances in Engineering Software. Volume 173, 103203.

