# ActiveNematics2D
[![License](https://img.shields.io/badge/License-MIT-blue.svg)] (LICENSE)

## Overview
This repository contains a C-based simulation of active nematic systems, a computational model for studying the dynamics of active matter, such as liquid crystals or biological systems with self-driven components. The code implements a lattice-based numerical simulation to model the evolution of the Q-tensor (describing nematic order), velocity field, and pressure field, using Euler updates and finite-difference methods. Key features include GPU-acceleratable algorithms (aligned with CUDA expertise from my resume) and efficient handling of large-scale datasets for physical systems.

## Features
- **2D Lattice Simulation**: Models active nematic systems on a 64x64 grid with periodic boundary conditions.
- **Q-Tensor Dynamics**: Computes the nematic order parameter (Q-tensor) and its evolution, capturing orientational order and defects.
- **Velocity and Pressure Fields**: Solves for incompressible flow using a pressure-Poisson equation, ensuring physical accuracy.
- **Numerical Methods**: Implements finite-difference schemes (Laplacian, divergence, upwind advection) for robust computation.
- **Data Output**: Saves simulation results (Q-tensor, velocity, pressure, saddle-splay energy, etc.) to text files for analysis and visualization.
- **Modular Design**: Structured functions for easy extension, e.g., adding GPU acceleration or new boundary conditions.


## Prerequisites
- **Operating System**: Linux, macOS, or Windows (with a C compiler like `gcc`).
- **Compiler**: GCC or any C99-compatible compiler.
- **Dependencies**: Standard C libraries (`stdlib.h`, `math.h`, `stdio.h`, `string.h`, `unistd.h`, `sys/types.h`, `sys/stat.h`, `stdbool.h`, `time.h`).
- **Optional for Visualization**: Python with NumPy/Matplotlib for post-processing simulation outputs (used in `create_plot.py` file).

## Installation
1. **Clone the Repository**:
```bash
git clone https://github.com/tomidiy/active-nematic-simulation.git
cd active-nematic-simulation
```

2. **Compile the Code**:
```bash
gcc -o active_nematic main.c -lm
```
- The `-lm` flag links the math library for functions like `sqrt` and `cos`.
- Ensure `gcc` is installed (sudo apt install gcc on Ubuntu or equivalent).

3. Create Output Directories:
- The code automatically creates `Results/` and `Images/` directories for simulation outputs. Ensure write permissions in the working directory.

## Usage
- **Run the Simulation**:
```bash
./active_nematic
```
-- The program initializes a 64x64 grid with random Q-tensor orientations and evolves the system for up to 10^7 steps or until convergence (`udiff_thresh` or `max_t` is reached).
-- Output files are saved in `Results/test/` with subdirectories for Q-tensor (`Q/`), velocity (`u/`), pressure (`p/`), saddle-splay (`ss/`), strain rate (`E/`), vorticity (`omega/`), and nematic director (`n_dor/`).

- **Key Parameters** (defined in `main.c`):
-- `Lx`, `Ly`: Grid size (64x64).
-- `dt`: Time step (2e-4).
-- `K`: Elastic constant (256^2).
-- `zeta`: Activity parameter (100.0).
-- `gamma`: Rotational viscosity (10 * 256).
-- `lambda`: Flow-alignment parameter (0.1).
-- `save_every_n_steps`: Frequency of saving outputs (every 100 steps).
-- Modify these in `main.c` to adjust simulation behavior.

- **Output Files**:
-- `Results/test/Q/Q_XXXXXXXXXX.txt`: Q-tensor components at each step.
-- `Results/test/u/u_XXXXXXXXXX.txt`: Velocity field components.
-- `Results/test/p_mean/AverageP_100.0.txt`: Average pressure over time.
-- Additional files for saddle-splay, vorticity, strain rate, and nematic director.
-- Use Python/Matplotlib to visualize outputs. create_plot.py can be used for visualization.

## File Structure

active-nematic-simulation/
├── main.c               # Main simulation code
├── Results/             # Output directory for simulation data
│   ├── test/            # Subdirectory for run labeled "test"
│   │   ├── Q/           # Q-tensor outputs
│   │   ├── u/           # Velocity field outputs
│   │   ├── p/           # Pressure field outputs
│   │   ├── ss/          # Saddle-splay energy outputs
│   │   ├── E/           # Strain rate tensor outputs
│   │   ├── omega/       # Vorticity outputs
│   │   ├── n_dor/       # Nematic director and degree of order
│   │   ├── p_mean/      # Average pressure outputs
│   │   └── test_consts.txt # Simulation parameters
├── Images/              # Directory for plots (empty by default)
└── README.md             # This file

## Key Functions
- **Laplacian()**: Computes the Laplacian of a scalar field for each lattice site.
- **Laplacian_vector()**: Computes the Laplacian of a vector field (e.g., Q-tensor, velocity).
- **div_vector()**: Calculates the divergence of the velocity field.
- **H_S_from_Q()**: Computes molecular field (H) and co-rotation tensor (S) for Q-tensor dynamics.
- **calculate_Pi()**: Calculates elastic and active stress tensors.
- **relax_pressure()**: Solves the pressure-Poisson equation to ensure incompressibility.
- **update_step()**: Performs Euler updates for Q-tensor, velocity, and pressure fields.
- **run_active_nematic_sim()**: Main simulation loop, handling initialization and output.

## Future Improvements
- **Customizable Boundary Conditions**: Implement apply_Q_boundary_conditions(), apply_u_boundary_conditions(), and apply_p_boundary_conditions() for flexible configurations.
- **Parallelization**: Extend with MPI for distributed computing on larger grids.



## Contributing
Contributions are welcome! To contribute:
- Fork the repository.
- Create a new branch (`git checkout -b feature-branch`).
- Make changes and commit (`git commit -m "Add feature"`).
- Push to your fork (`git push origin feature-branch`).
- Open a pull request with a clear description of changes.
Please ensure code follows the existing style (e.g., consistent indentation, clear comments) and includes tests if applicable.

## License
This project is licensed under the MIT License (LICENSE).


