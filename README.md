# Solution Space Optimization (SSO) Toolbox

## Overview
This is a MATLAB-based toolbox dedicated to the computation and optimization of solution spaces.

## Features
- Optimization of box-shaped solution spaces.
- Optimization of component solution spaces.
- Interfaces for both MATLAB- and Python-based system functions (bottom-up mappings).
- Many smaller functions for the visualization of said solution spaces.
- Some system-specific solvers (like for trusses).

## Project Structure
- 'BottomUpMappings': system functions written in MATLAB.
- 'BottomUpMappingsPython': system functions written in Python.
- 'Data': general space for extra inputs.
- 'SolutionSpacesToolbox': functions specifically related to the computation of solution spaces.
- 'SupportFunctions': specific wrappers for generic optimization functions and general system solvers.
- 'TestScripts': scripts where the functions are used to find the desired results.
- 'UtilityFunctions': functions with general usefulness (not intrinsically tied to solution spaces).

## Full Documentation
Each function has an implemented header, which can be seen in matalb using 'help'. These headers contain all relevant information for each function. 

## Requirements
MATLAB 2019b or newer.

## Installation
Clone or download the repository and put the files ont he active workspace being used.

## Usage
To add this functions to the MATLAB path, run 'setup_sso_toolbox'.
For examples on how to use each functionality, check the tests in 'TestScripts'.

## Authors and acknowledgment
Eduardo Rodrigues Della Noce
From 10.01.2022 to 09.01.2025, the development of this code was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Project Number 454149634.

## License
Licensed under Apache License, Version 2.0.
https://www.apache.org/licenses/LICENSE-2.0