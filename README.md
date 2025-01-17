# Solution Space Optimization (SSO) Toolbox

## Overview
The SSO Toolbox is a MATLAB-based collection of functions and scripts dedicated to the computation and optimization of solution spaces.

## Features
- Optimization of box-shaped solution spaces  
- Optimization of component-specific solution spaces  
- Interfaces for both MATLAB- and Python-based system functions (bottom-up mappings)  
- Multiple auxiliary functions for visualizing solution spaces  
- Specialized solvers for certain systems (e.g., trusses)

## Project Structure
- *BottomUpMappings*: System functions written in MATLAB  
- *BottomUpMappingsPython*: System functions written in Python  
- *Data*: General folder for extra inputs  
- *SolutionSpacesToolbox*: Functions specifically related to computing solution spaces  
- *SupportFunctions*: Wrappers for generic optimization functions and general system solvers  
- *TestScripts*: Scripts demonstrating how to use the toolbox functions for desired outcomes  
- *UtilityFunctions*: General-purpose functions not intrinsically tied to solution spaces

## Full Documentation
Each function includes a header that can be accessed in MATLAB using the `help` command. These headers contain all relevant information for each function.

## Requirements
- MATLAB 2019b or newer

## Installation
Clone or download the repository and place all files in your active MATLAB workspace.

## Usage
To add these functions to the MATLAB path, run `setup_sso_toolbox`.  
For usage examples and demonstrations, refer to the scripts in the **TestScripts** folder.

## Authors and Acknowledgment
**Eduardo Rodrigues Della Noce**  
From 10.01.2022 to 09.01.2025, the development of this code was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – Project Number 454149634.

## License
Licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).
