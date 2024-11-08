# CAPy: Circadian Analysis Python

CAPy is a graphical user interface (GUI) application designed for analyzing circadian activity data. Built with PyQt5, this tool allows researchers and enthusiasts to process, visualize, and interpret activity data with ease.

## Table of Contents

- [Features](#features)
- [Screenshots](#screenshots)
- [Installation](#installation)
- [Usage](#usage)
- [Prerequisites](#prerequisites)
- [License](#license)
  
## Features

- **Data Import**: Load main activity data along with optional onset and offset clock data files.
- **Session Management**: Save and load sessions to preserve your work.
- **Analysis Options**:
  - Generate single and double actograms.
  - Perform Fourier analysis to determine dominant period.
  - Calculate onset, offset, acrophase, and bathyphase times.
  - Save results to CSV for further analysis.
- **Interactive Plots**:
  - Drag and adjust markers on actograms to refine analysis.
  - Real-time updates to results upon adjustments.
- **Customizable Parameters**:
  - Set threshold percentiles for activity detection.
  - Define inactivity and activity durations (N and M hours).
  - Choose to use calculated free running period (tau) or manually set values.


## Screenshots

*(Screenshots can be added here to showcase the GUI and features.)*

## Installation
  We recomend using anaconda for ease of instillation
1. **Clone the Repository**

    ```
    git clone https://github.com/yourusername/CAPy.git
    cd CAPy
    ```

2. **Create a Virtual Environment with Dependencies**

    ```
    conda env create -f CAPy_env.yaml
    conda activate CAPy
    ```

3. **Run the Application**

    ```
    python gui.py
    ```

## Usage

1. **Launch the Application**

    Run `python gui.py` in your terminal.

2. **Load Data Files**

    - Click on **"Select Main Data File"** to load your primary activity CSV file.
    - Optionally, load onset and offset clock data files for comparisons.

3. **Configure Analysis Parameters**

    - Choose to use the calculated tau or manually enter a tau value.
    - Set threshold percentile, N hours of inactivity, and M hours of activity.
    - Select tasks such as generating actograms, plotting Fourier analysis, saving results to CSV, or performing comparisons.

4. **Run Analysis**

    Click on **"Run Analysis"** to perform the selected tasks.

5. **Interact with Results**

    - View generated plots in the **"Plots"** tab.
    - Drag markers on actograms to adjust onset, offset, acrophase, and bathyphase times.
    - Click **"Update Results with Adjusted Times"** to refresh the analysis based on adjustments.

6. **Save and Load Sessions**

    - Use **"Save Session"** to save your current state.
    - Use **"Load Session"** to resume work from a previous session.

## Prerequisites

- **Python 3.6 or higher**
- **Packages**:
  - PyQt5
  - pandas
  - matplotlib
  - numpy
  - scipy

*Note: Ensure all dependencies are installed as per the [Installation](#installation) section.*

## License

This project is licensed under the [MIT License](LICENSE).

