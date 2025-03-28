# OrcaFunctionalHub

## Overview

**OrcaFunctionalHub** is a Julia script designed to facilitate molecular modeling and computational chemistry tasks using the Orca quantum chemistry package. The script automates the conversion of SMILES representations of molecules into 3D coordinates, performs geometry optimization, and executes time-dependent density functional theory (TDDFT) calculations. 

## Features

- **3D Coordinate Generation**: Converts SMILES formulas into 3D coordinates using RDKit.
- **Query Generation**: Creates and writes Orca input files for different computational tasks.
- **Output Parsing**: Parses output files from Orca to extract relevant information such as optimized coordinates and spectra data.
- **Data Visualization**: Generates and saves absorbance spectra graphs based on TDDFT calculations.
- **Configurable Input**: Uses YAML configuration files to specify molecular properties and calculation parameters.

## Requirements

- **Julia**: The script is written in Julia and requires a Julia environment.
- **Packages**: The following Julia packages are needed:
  - `FilePathsBase`
  - `LaTeXStrings`
  - `DataFrames`
  - `PyCall`
  - `Plots`
  - `YAML`
  - `CSV`
- **RDKit**: The RDKit library must be installed and accessible through Python.

## Installation

1. **Install Julia**: Download and install Julia from [the official website](https://julialang.org/downloads/).
2. **Install Required Packages**: Use the Julia package manager to install the required packages:
   ```julia
   using Pkg
   Pkg.add("FilePathsBase")
   Pkg.add("LaTeXStrings")
   Pkg.add("DataFrames")
   Pkg.add("PyCall")
   Pkg.add("Plots")
   Pkg.add("YAML")
   Pkg.add("CSV")
   ```
3. **The RDKit library must be installed and accessible through Python** You can install RDKit using conda:
```bash
conda install -c conda-forge rdkit
```
Confirm that RDKit is installed correctly by opening a Python session and importing the library:
```python
import rdkit
print(rdkit.__version__)
```
## Usage

### Create a YAML Configuration File

Specify the molecules and calculation parameters in a YAML file. An example is provided below:

```yaml
# Skye Goetz (CalPoly) 10/22/2024

Molecules : 
# [[ list of molecules and attributes ]]
  # smiles formulas REQUIRE explicit hydrogens

 - SmilesFormula : "[CH3][CH3]"
   MolecularCharge : 0 
   MolecularMultiplicity : 1

TightOpt : 
# [[ TightOpt Query Parameters ]]

 Functional : PBE0
 BasisSet : Def2-SVPD
 SolventModel : SMD
 SolventToUse : ethanol

ScanAngle : 
# [[ Scan Angle Query Parameters ]]

 Functional : PBE0
 BasisSet : Def2-SVPD
 SolventModel : SMD
 SolventToUse : ethanol

 ScanType : D
 Angles : 
 # [[ List Of Angles For Geometry Scan ]]

   - 0
   - 1
   - 4
   - 6

 Advanced : 

   ScanRange :
   # [[ In Degrees ]]
     
     - 0
     - 90

   NumberOfScans : 10

TDDFT : 
# [[ TDDFT Query Parameters ]]

 Functional : CAM-B3LYP
 BasisSet : Def2-TZVPD
 SolventModel : CPCM
 SolventToUse : ethanol

 Advanced : 

   Nroots : 10
   Maxdim : 100

AdvancedSettings : 
# [[ Configures OrcaFunctionalHub ]]
 
 NumberOfGeometryScansPerTDDFTCalculation : 1
```

## Run the Script

Execute the Julia script with the path to your YAML file:

```bash
julia OrcaFunctionalHub.jl path/to/config.yaml
```

## Output

The script generates several output files, including:

- Orca input files for each calculation.
- XYZ files containing molecular coordinates.
- Output files containing results from Orca.
- Graphs of absorbance spectra saved as PNG files.

### Author

Skye Goetz (CalPoly)  
Date: 10/22/2024
