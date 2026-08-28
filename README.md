

# USPEX-GP

## Introduction

USPEX-GP is an open-source computational materials science package designed for crystal structure prediction using evolutionary algorithms (USPEX-9.4.4). It supports various computational codes (VASP, Quantum ESPRESSO, LAMMPS, etc.) and is designed for high-performance computing environments. This repository contains a hybrid codebase with both MATLAB main scripts and Python helper tools/extensions.

## Features

* **Crystal Structure Prediction**: Predicts stable crystal structures given only chemical composition.
* **Evolutionary Algorithm (EA)**: Uses a robust genetic algorithm with local optimization.
* **Particle Swarm Optimization (PSO)**: Alternative optimization algorithm included.
* **Variable Composition Search**: Supports searching across different chemical compositions.
* **Molecule & Surface Predictions**: Specialized modes for molecular crystals and surface structures.
* **Metadynamics**: Enhanced sampling capabilities for complex energy landscapes.
* **Multi-Code Support**: Interfaces with popular DFT and MD codes.
* **Hybrid MATLAB/Python Architecture**: Combines MATLAB's numerical power with Python's flexibility.

## Installation

### Prerequisites

* **OCTAVE** with toolboxes:
  * Optimization Toolbox
  * Statistics and Machine Learning Toolbox
* **Python** (3.7 or newer)
* **Cython** (for compiling core extensions)
* **Numerical Python (NumPy)**
* **C Compiler** (e.g., GCC, MinGW) compatible with your Python/MATLAB setup.

### Steps

1. **Clone the repository**:
   ```bash
   git clone <repository_url>
   cd uspex-gp
   ```

2. **Setup Python Environment**:
   It is recommended to use a virtual environment.
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install numpy cython
   ```

3. **Compile Cython Extensions**:
   Navigate to the `FunctionFolder/sys` directory and compile the extensions (refer to specific `.pyx` files or local build instructions if provided).
   ```bash
   python setup.py build_ext --inplace
   ```

4. **MATLAB Configuration**:
   Add the `FunctionFolder` and its subdirectories to your MATLAB path. You can do this via the MATLAB interface or by adding a path file in your startup script.

## Usage

1. **Configure Your Calculation**:
   Modify the `INPUT.txt` file (or equivalent configuration) in your working directory. This file controls all USPEX parameters (population size, number of generations, fitness function, etc.).

2. **Select Fitness Calculator**:
   Ensure your `calc.exe` (or similar) script is correctly set up to call your chosen DFT/MD code (e.g., VASP, Quantum ESPRESSO) based on the generated structure files.

3. **Launch USPEX**:
   Run the main USPEX script from OCTAVE:
   ```octave
   >> USPEX
   ```
   Or use the Python wrapper if available:
   ```bash
   python run_uspex.py
   ```

4. **Monitor & Analyze**:
   Check generation folders (`generation_X`) for output files. Use the visualization scripts in `FunctionFolder/USPEX/` to analyze results.

## File Structure

* **FunctionFolder/**: Main source code.
    * **AbinitCode/**: Interfaces and I/O for various simulation codes (VASP, QE, LAMMPS, etc.).
    * **USPEX/**: Core evolutionary algorithm modules and structure prediction logic (`.m` files).
    * **PSO/**: Particle Swarm Optimization implementation.
    * **METADYNAMICS/**: Metadynamics-related functions.
    * **sys/**: System utility functions, core Cython extensions (`.pyx`, `.c`), and Python helpers.
    * **Tool/**: External tools and scripts.
    * **spacegroup/**: Symmetry analysis tools.
* **README.md**: This file.

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## License

This project is licensed under the [LICENSE NAME] License - see the [LICENSE](LICENSE) file for details. As an open-source project on Gitee, please ensure you comply with the specific terms of its license.

## Contact

* **Author**: Fenggo
* **Gitee Profile**: [https://gitee.com/fenggo](https://gitee.com/fenggo)
* **Project Home**: [https://gitee.com/fenggo/uspex-gp](https://gitee.com/fenggo/uspex-gp)

---
*Note: This README was generated based on the repository structure. Please verify specific details with the code comments and documentation within the repository.*