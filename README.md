# Experimental-Modal-Analisys

# Structural Modal Analysis: Bridging Numerical Simulation and Real-World Data

This repository features a comprehensive structural dynamics study focused on characterizing a metallic plate. The primary goal was to thoroughly determine its dynamic properties: **natural frequencies**, **damping factors**, and the complex **mode shapes**. The project elevates standard analysis by robustly combining the power of **Finite Element Method (FEM)** simulations with experimental data gathered via impact hammer testing, creating a validated workflow from theory to practice.

## Key Methodologies and Project Components

* **Finite Element Modeling (FEM):** A detailed model of the plate was developed, strictly adhering to its precise geometries and boundary conditions. This numerical model served as a crucial predictive tool, providing an initial baseline for modal frequency prediction.

* **Experimental Data Acquisition:** Using an instrumented hammer, the structure was excited at strategic points, capturing the resulting response with accelerometers. The raw time-domain data was then carefully processed to derive the **Frequency Response Functions (FRFs)**, which form the foundation of the experimental analysis.

* **Advanced Modal Identification Techniques:** Several powerful methodologies were employed to extract dynamic parameters from the measured FRFs:
    * **SDOF Analysis (Single Degree of Freedom):** Used for the analysis of isolated modes, applying the circle fit approximation to FRF curves near resonance peaks for a quick estimation of frequency and damping.
    * **MDOF Analysis (Multiple Degree of Freedom):** To handle closely coupled or complex modes, more robust MDOF techniques were implemented, including:
        * **ITD Method (Ibrahim's Time Domain):** A time-domain technique that accurately derives modal parameters from the free decay response.
        * **Prony's Method:** Another time-domain approach that models the signal as a sum of decaying sinusoids, directly yielding frequencies, damping ratios, and amplitudes.

* **Bimodal Validation and Correlation:** The project culminates in a rigorous validation phase. We compared the mode shapes derived from the FEM model with those extracted experimentally, quantifying the agreement using the **Modal Assurance Criterion (MAC)**. The MAC matrix provides an essential numerical and visual indicator to confirm the consistency between the simulation and measurement results.

* **Data Visualization:** The repository includes scripts for generating graphs and animations that make the complex dynamic behavior of the plate easy to understand, including FRF curves, Nyquist plots (modal circles), and animations of the identified mode shapes.

## Tools Employed

* **MATLAB:** The primary platform for data processing, modal analysis, and visualization, leveraging its powerful capabilities in signal processing and system identification.
* **FEM Software (Abaqus):** Used for the creation and execution of the initial finite element simulation.

## Intended Audience

This repository is a valuable resource for **engineers**, **researchers**, and **students** in mechanical, civil, and aerospace engineering fields. It offers a complete, practical overview for anyone looking to learn or actively apply the principles of experimental modal analysis, structural dynamics, and system identification techniques.
