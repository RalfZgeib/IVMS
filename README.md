# IVMS: Low-Cost, Load Cell-Based In-Situ Viscosity Measurement System for Extrusion-Based Bioprinting

This repository contains the code, hardware designs, and documentation for the **In-situ Viscosity Measurement System (IVMS)**, a low-cost, modular platform for real-time measurement of bioink viscosity during extrusion-based bioprinting. IVMS was developed to enable continuous, accurate rheological characterization of shear-thinning hydrogels (such as GelMA) directly during the bioprinting process, supporting rapid optimization of printing parameters and improved reproducibility in tissue engineering research.

---

## Table of Contents

- [Overview](#overview)
- [System Features](#system-features)
- [Hardware Setup](#hardware-setup)
- [Software Setup](#software-setup)
- [Usage Instructions](#usage-instructions)
- [Data Analysis](#data-analysis)
- [Validation](#validation)
- [File Structure](#file-structure)
- [Dependencies](#dependencies)
- [Citing This Work](#citing-this-work)
- [License](#license)
- [Contact](#contact)

---

## Overview

**IVMS** integrates an S-beam load cell into a custom, 3D-printed extrusion head on an open-source bioprinter. The system records extrusion forces during printing, converts them to pressure, and—using MATLAB-based scripts—applies rheological corrections (including friction, Bagley, wall-slip, and Rabinowitsch corrections) to estimate true, shear-rate-dependent viscosity in real time. The approach is validated against standard parallel plate rheometry across multiple GelMA concentrations and temperatures, demonstrating close agreement and robust performance for in-situ monitoring.

For full experimental details, see the associated publication in *ACS Biomaterials Science & Engineering*:

> **Development of a Low-Cost, Load Cell-Based, In-Situ Viscosity Measurement System for Extrusion-Based Bioprinting**  
> Ralf Zgeib, Xiao Zhao, Ahmadreza Zaeri, Fucheng Zhang, Kai Cao, Robert Chang  
> *ACS Biomater. Sci. Eng.*, 2025  
> [DOI: to be added upon publication]

---

## System Features

- **Real-time viscosity measurement** during extrusion-based bioprinting
- Supports both in-air and support bath printing paradigms
- Modular, open-source hardware design with 3D-printable components
- Accurate force-to-viscosity conversion with comprehensive rheological corrections
- MATLAB-based data acquisition and analysis pipeline
- Validated using GelMA hydrogels at multiple concentrations and temperatures

---

## Hardware Setup

### Materials

- **S-beam load cell** (e.g., DYLY106 Micro Size 5 kg)
- **NEMA 17 stepper motor**
- **3 mL syringes** (BD Disposable Syringes with Luer-Lok Tips)
- **3D-printed parts** (STL files in `/hardware/`)
- **Fasteners:** M6x10, M3x20 screws, M3 nuts, etc.
- **Threaded rod:** M3x0.5x300mm
- **High vacuum grease, cyanoacrylate adhesive**
- **FDM 3D printer** (for fabricating custom parts)
- **Basic tools** (screwdrivers, Allen wrenches, pliers, sandpaper)

### Assembly

1. **Print all 3D parts** from the `/hardware/` directory.
2. **Assemble the extrusion head** according to the instructions in `/docs/AssemblyGuide.md`.
3. **Mount the load cell** above the syringe using the provided coupler.
4. **Connect the load cell** to the signal conditioning electronics (see `/electronics/`).
5. **Integrate the extrusion system** with your open-source bioprinter.

For detailed step-by-step instructions and diagrams, see `/docs/AssemblyGuide.md`.

---

## Software Setup

- **MATLAB** (tested with R2022b and later)
- Required custom scripts: see `/code/`
- Load cell calibration and signal acquisition scripts: `/code/calibration/`
- Data processing and viscosity calculation scripts: `/code/analysis/`

### Installation

1. Clone this repository:
git clone https://github.com//IVMS.git
cd IVMS
2. Install MATLAB and required toolboxes.
3. Open `/code/IVMS_main.m` in MATLAB.

---

## Usage Instructions

1. **Calibrate the load cell** using the provided calibration script and known weights.
2. **Prepare your bioink** (e.g., GelMA at desired concentration and temperature).
3. **Load the bioink** into the syringe and install it in the extrusion head.
4. **Set extrusion parameters** (flow rate, nozzle geometry) via the bioprinter’s g-code or control software.
5. **Start the MATLAB acquisition script** to record extrusion force data in real time.
6. **Perform a print**; force data will be logged automatically.
7. **Process the data** using `/code/analysis/IVMS_viscosity_calc.m` to obtain shear rate and viscosity profiles.

---

## Data Analysis

- **Force data** is converted to pressure using the known syringe and nozzle geometry.
- **Shear rate** is calculated based on flow rate and nozzle dimensions.
- **Rheological corrections** (Bagley, wall-slip, Rabinowitsch) are applied for accurate viscosity estimation.
- **Viscosity curves** are fitted to the Carreau-Yasuda model for non-Newtonian fluids.
- **Results** can be compared to offline rheometer data for validation.

Sample data and analysis scripts are provided in `/data/` and `/code/analysis/`.

---

## Validation

The IVMS system was validated using GelMA hydrogels at 5%, 7%, and 10% w/v, and at 25°C, 30°C, and 37°C. IVMS measurements matched rheometer results within acceptable error margins, capturing the expected shear-thinning behavior and temperature/concentration dependencies.

For detailed validation results, see the `Validation` section in the associated publication and `/docs/ValidationReport.pdf`.

---

## File Structure

/hardware/ # 3D printable part files (STL, CAD)
/electronics/ # Circuit diagrams and BOM for load cell signal conditioning
/code/ # MATLAB scripts for calibration, data acquisition, and analysis
/data/ # Sample raw and processed data files
/docs/ # Assembly guide, user manual, validation report
LICENSE
README.md


---

## Dependencies

- MATLAB (R2022b or later)
- Instrument Control Toolbox (for data acquisition)
- [Optional] Python 3.x for additional data processing (see `/code/python/`)
- Standard open-source 3D printing and CAD software (e.g., Fusion 360)

---

## Citing This Work

If you use IVMS in your research, please cite:

> Zgeib, R., Zhao, X., Zaeri, A., Zhang, F., Cao, K., Chang, R.  
> Development of a Low-Cost, Load Cell-Based, In-Situ Viscosity Measurement System for Extrusion-Based Bioprinting.  
> *ACS Biomaterials Science & Engineering*, 2025.  
> [DOI: to be added]

---

## License

This project is released under the MIT License. See `LICENSE` for details.

---

## Contact

For questions, feedback, or contributions, please contact:

- **Robert Chang** (corresponding author): rchang6@stevens.edu
- Or open an issue in this repository.

---

**Acknowledgments:**  
This work was supported by the BMBM Lab, Department of Mechanical Engineering, Stevens Institute of Technology. For further details, see the full methods in the published manuscript.

---

*This README and all repository contents are intended to support reproducibility and open science in extrusion-based bioprinting research. Please contribute improvements and report issues via GitHub.*
