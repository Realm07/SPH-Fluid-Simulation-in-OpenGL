# Real-Time 2D SPH Fluid Simulation

[![Language](https://img.shields.io/badge/Language-C%2B%2B-blue.svg)](https://isocpp.org/)
[![API](https://img.shields.io/badge/Graphics-OpenGL-blue.svg)](https://www.opengl.org/)
[![GUI](https://img.shields.io/badge/GUI-Dear%20ImGui-orange.svg)](https://github.com/ocornut/imgui)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A high-performance, interactive 2D fluid simulation built with modern C++ and OpenGL. This project implements the **Smoothed Particle Hydrodynamics (SPH)** method to simulate complex fluid behaviors in real-time, complete with a dynamic GUI for parameter tuning and live experimentation.

This work is heavily inspired by the fantastic tutorials and projects by **[Sebastian Lague](https://www.youtube.com/c/SebastianLague)**.

![Project Screenshot](https://i.ibb.co/MD3xQyzL/Cost-Effective-Interactive-2-D-Fuild-Simulation.png)

## About The Project

This project explores the fascinating world of fluid dynamics by implementing a particle-based simulation from scratch. It serves as a practical application of physics principles, numerical methods, and performance optimization techniques in C++.

By approximating the continuous fluid as a large set of discrete particles, we can simulate complex phenomena like pressure waves, viscosity, and surface tension. The goal was to create a simulation that is not only visually engaging but also highly interactive, allowing users to build an intuition for how different physical parameters affect fluid behavior.

### Key Features

*   **SPH Physics Core:** Simulates pressure, viscosity, and near-pressure (for surface tension) forces based on the principles of Smoothed Particle Hydrodynamics.
*   **High Performance:** Achieves real-time frame rates with thousands of particles by using an optimized **Spatial Grid Hashing** system with **Radix Sort** for extremely fast neighbor lookups.
*   **Modern C++:** Leverages C++17 features, including parallel execution policies (`std::execution::par`) for updating particle positions.
*   **Fully Interactive GUI:** Powered by **Dear ImGui**, allowing live manipulation of every core simulation parameter.
*   **Dynamic Obstacle Editor:** Add, remove, and resize obstacles in real-time to see how the fluid reacts.
*   **Rich Visualization:**
    *   Particles are colored by velocity magnitude (from blue for slow to red for fast) for at-a-glance analysis.
    *   Multiple debug overlays to visualize the smoothing radius, density fields, the spatial grid, and force vectors.

## How It Works

The simulation approximates a continuous fluid by representing it as a large number of discrete particles. The core physics loop is driven by the Navier-Stokes equations for incompressible flow.

1.  **Neighbor Search:** For each particle, an efficient neighbor search is performed using a spatial hashing grid. This avoids the O(n²) complexity of checking every other particle and is the key to real-time performance.
2.  **Density Calculation:** The density at each particle's position is calculated by summing the mass of its neighbors, weighted by a `SmoothingKernel` function. A separate `NearDensityKernel` is used to compute forces that create surface tension.
3.  **Pressure Calculation:** Density is converted into pressure. The pressure gradient creates forces that push particles from high-density areas to low-density areas, mimicking fluid behavior.
4.  **Force Calculation:** In addition to pressure, forces for viscosity (internal friction) and external forces (like gravity and mouse interaction) are computed.
5.  **Integration:** The final forces are used to update each particle's velocity and position over time using a simple Euler integration step.

## Tech Stack

*   **Language:** C++17
*   **Graphics API:** OpenGL (via GLEW)
*   **GUI Library:** [Dear ImGui](https://github.com/ocornut/imgui)
*   **Windowing & Input:** [GLFW](https://www.glfw.org/)

## Getting Started

To get a local copy up and running, you will need a C++ development environment and the required libraries (GLEW, GLFW).

### Prerequisites

*   A C++ compiler (MSVC, GCC, Clang)
*   [GLEW](http://glew.sourceforge.net/)
*   [GLFW](https://www.glfw.org/)
*   The project source code includes the necessary headers for [Dear ImGui](https://github.com/ocornut/imgui).

### Building
1.  **Clone the repo:**
    ```sh
    git clone https://github.com/Realm07/SPH-Fluid-Simulation-in-OpenGL.git
    cd SPH-Fluid-Simulation-in-OpenGL
    ```
2.  **Set up your build environment:** Ensure your IDE (like Visual Studio) or build system (like Make with g++) is configured to link against OpenGL, GLEW, and GLFW libraries.
3.  **Compile the source files:** Compile `main.cpp` along with the ImGui implementation files (`imgui.cpp`, `imgui_draw.cpp`, `imgui_tables.cpp`, `imgui_widgets.cpp`, `imgui_impl_glfw.cpp`, `imgui_impl_opengl3.cpp`).

## Usage

Once running, you can interact with the simulation using the GUI windows.

### Config Window
*   **Gravity:** Controls the strength of the downward gravitational force.
*   **Smoothing Radius:** The influence radius of each particle, critical for density and force calculations.
*   **Pressure Multiplier:** A stiffness constant. Higher values create stronger, more reactive pressure forces.
*   **Target Density:** The density the fluid strives to maintain at rest.
*   **Viscosity Strength:** Controls the fluid's "thickness" or internal friction.
*   **Near Pressure Mult:** Adjusts the strength of the surface tension effect.
*   **Mouse Interaction:** Tweak the radius and strength of the interactive force applied with the mouse cursor.
*   **Visualization Toggles:** Enable or disable various debug overlays.
*   **Reset Buttons:** Instantly reset the particles to an ordered grid or a random distribution.

### Obstacle Editor
*   Enable the **"Obstacle Editor"** checkbox to open its window.
*   Click **"Add New Obstacle"** to create a new rectangular body.
*   Use the sliders to adjust an obstacle's **Center**, **Width**, and **Height** and see the fluid react in real-time.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

This means you are free to use, modify, and distribute this code for any purpose, including commercial projects.

### Attribution (A Friendly Request)

While the MIT license does not require it, if you find this project useful for your own work, I would greatly appreciate a credit or a link back to this repository. Thank you!

## Acknowledgements

*   **[Sebastian Lague](https://www.youtube.com/c/SebastianLague)** for the invaluable educational content and inspiration.
*   **[Dear ImGui](https://github.com/ocornut/imgui)** for the fantastic and easy-to-use GUI library.
*   The communities behind **OpenGL**, **GLEW**, and **GLFW**.
#### [Application Download](https://github.com/Realm07/SPH-Fluid-Simulation-in-OpenGL/releases)
# Please Note:
When running the application, windows will prompt it as untrusted and malicious. That is because my code isn't signed. To get it signed I'll need to pay ~350$ an year. **The code is clean.** You can view it in the repository under `src/`. Alternatively you can build and compile it yourself.
