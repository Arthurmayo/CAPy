# ActivityAnalysisGUI

**ActivityAnalysisGUI** is a Python-based graphical user interface (GUI) application designed for analyzing activity data using actograms and Fourier analysis. This tool is particularly useful for researchers and analysts working with time-series event count data, enabling them to visualize patterns, detect key activity periods, and compare experimental results with clock data.

## Features

- **Data File Management**
  - **Main Data File:** Load primary CSV files containing event counts.
  - **Onset/Offset Clock Data Files:** Optionally load additional CSV files for comparative analyses.
  - **CSV Save Path:** Specify where to save the analysis results.

- **Tau Selection**
  - **Calculated Tau:** Automatically determine the dominant period using Fourier analysis.
  - **Manual Tau Entry:** Manually input the tau value if preferred.

- **Activity Parameters**
  - **Threshold Percentile:** Define the percentile for activity thresholding.
  - **Inactivity/Activity Hours:** Configure hours of inactivity (`N_hours_inactivity`) and activity (`M_hours_activity`).

- **Task Selection**
  - **Generate Actograms:** Create single or double actograms to visualize activity patterns.
  - **Plot Fourier Analysis:** Display the Fourier power spectrum to identify dominant periods.
  - **Save to CSV:** Export analysis results to a CSV file.
  - **Perform Comparisons:** Compare detected onset and offset times with experimental clock data.

- **Interactive Plots**
  - **Draggable Markers:** Adjust onset, offset, acrophase, and bathyphase times directly on the actograms.
  - **Dynamic Legends:** View amplitudes of dominant periods in plot legends.

- **Session Management**
  - **Save Session:** Preserve current analysis settings and marker positions.
  - **Load Session:** Restore previously saved sessions for continued analysis.

## Installation

### Prerequisites

- **Python 3.x**: Ensure you have Python installed. You can download it from [python.org](https://www.python.org/downloads/).

### Dependencies

Install the required Python packages using `pip`:

```bash
pip install PyQt5 pandas matplotlib scipy
