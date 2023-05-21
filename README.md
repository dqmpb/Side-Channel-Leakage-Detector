# Side-Channel-Leakage-Detector

# Program Summary

This program is a Python-based tool for side-channel leakage analysis, specifically designed for analyzing power consumption traces in cryptographic implementations. It uses Python 3.9 and the following libraries: h5py, scipy, matplotlib, numpy, math, and tqdm.

## Datasets

The program was tested and developed using the following public datasets:

1. **ASCAD**: This dataset contains power traces from ATMEGA and STM32 microcontrollers with different masking schemes and jitters. The program uses only the fixed key measurements from this dataset. The dataset has 50,000 traces with 700 trace points each. https://github.com/ANSSI-FR/ASCAD

2. **AES HD**: This dataset contains traces measured on a Xilinx Virtex-5 FPGA. All samples in this dataset are encrypted with the same key. The dataset has 45,000 traces with 1250 trace points each. https://github.com/AISyLab/AES_HD_Ext

3. **AES RD**: This dataset contains traces with random delays, simulating different clock cycles. The samples were obtained from an 8-bit Atmel AVR microcontroller. The dataset requires pre-processing for an effective analysis. The dataset has 50,000 traces with 3500 trace points each. https://github.com/ikizhvatov/randomdelays-traces

The program would also work on any file that follows the ASCAD's HDF5 format or the MATLAB format of the AES_RD.

## Features

1. **HDF5 File Support**: The program accepts HDF5 files for analysis and also provides an option to convert MATLAB files to HDF5 format.

2. **Leakage Detection Methods**: The program offers two statistical methods for detecting leakage in the power traces:

   - *Welch's T-test*: The test compares the values of the traces that correspond to a specific intermediate value (S-box output) to the values of the traces that correspond to all other possible intermediate values.
   
   - *Pearson's Chi2-test*: The test compares the observed distribution of the traces that correspond to a specific intermediate value (S-box output) to the expected distribution assuming no difference between the values.

3. **Visualization and Data Export**: The program generates and saves various plots and CSV files for the analysis results, including t-values, p-values, chi2 values, and degrees of freedom in the chi2 test.

## Installation Guide for Side-Channel-Leakage-Detector

**Prerequisites**

    Python 3.x installed on your system.
    Git installed on your system (optional).

**Step 1: Obtain the Files**

You can either clone the repository or directly download the necessary files from GitHub:

 Clone the repository:

      git clone https://github.com/dqmpb/Side-Channel-Leakage-Detector.git

      cd Side-Channel-Leakage-Detector

 *Or download the leakage_analyser.py and requirements.txt files directly from the GitHub repository.*

**Step 2: Install Dependencies (Optional: Create a Virtual Environment)**

   Install the required dependencies using pip:

      pip install -r requirements.txt

   You can also install the dependencies listed in requirements.txt individually:

      pip install matplotlib numpy scipy

**Step 3: Run the Leakage Analyser

      python leakage_analyser.py

That's it! You have successfully installed and set up the Side-Channel-Leakage-Detector project on your local machine. You may also consider creating a virtual environment to isolate dependencies specific to this project from your global Python installation.


## Usage

1. Run the program and choose an option:
   - `[0]` Submit HDF5 file for analysis
   - `[1]` Convert MATLAB file to HDF5 file
   - `[2]` Exit

2. Select a trace group for analysis.

3. Choose a leakage detection test:
   - `[0]` Leakage detection with Welch's T-test
   - `[1]` Leakage detection with Pearson's chi2-test
   - `[2]` Run both tests - option `[0]` + option `[1]`
   - `[3]` Exit

4. Select a target byte (0-15).

5. The program will run the selected test(s) and save the generated diagrams and CSV files.
