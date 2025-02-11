import sys
import os
from PyQt5 import QtCore, QtGui, QtWidgets
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import json
from activity_analysis import ActivityAnalysis
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from converter import convert_awd_to_csv

class ActivityAnalysisGUI(QtWidgets.QMainWindow):
    """
    A PyQt-based GUI for loading, analyzing, and plotting circadian activity data.
    """
    def __init__(self):
        """
        Initialize the main GUI application window and all its widgets.
        """
        super().__init__()

        # Initialize variables to hold file paths and other inputs
        self.file_path = None
        self.csv_save_path = None  

        self.use_calculated_tau = True
        self.manual_tau = None

        # Initialize variables for task selection checkboxes
        self.perform_single_actogram = False
        self.perform_double_actogram = False
        self.plot_fourier = False
        self.create_dataframe = False

        # Initialize variables for new input fields
        self.threshold_percentile = 30.0  
        self.N_hours_inactivity = 12.0
        self.M_hours_activity = 12.0

        # Initialize variables for label options in actograms
        self.label_onset = False
        self.label_offset = False
        self.label_acrophase = False
        self.label_bathyphase = False

        # Initialize variables for activity profile options
        self.display_com = False
        self.display_sem = False
        self.base_on_tau = True  

        self.marker_dict = {}  

        # NEW: A boolean to remember if user wants a 12-hour shift
        self.shift_12h = False

        self.initUI()

    def initUI(self):
        """
        Build and lay out all GUI elements (buttons, tabs, panels, etc.).
        """
        # Set window title and size
        self.setWindowTitle('CAPy')
        self.resize(1200, 700)

        icon_path = r'capy.jpg'
        if os.path.exists(icon_path):
            self.setWindowIcon(QtGui.QIcon(icon_path))

        # Main layout
        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QtWidgets.QVBoxLayout(central_widget)

        # 1. Convert/Clean AWD to CSV Button at the Very Top
        awd_layout = QtWidgets.QHBoxLayout()
        convert_awd_button = QtWidgets.QPushButton("Convert/Clean AWD to CSV")
        convert_awd_button.clicked.connect(self.convert_awd_file)
        awd_layout.addWidget(convert_awd_button)
        main_layout.addLayout(awd_layout)

        # Create tabs
        tabs = QtWidgets.QTabWidget()
        main_layout.addWidget(tabs)

        # Main Tab
        main_tab = QtWidgets.QWidget()
        tabs.addTab(main_tab, "Main")
        main_layout_tab = QtWidgets.QGridLayout(main_tab)

        # Session Layout (Save/Load)
        session_layout = QtWidgets.QHBoxLayout()
        save_session_button = QtWidgets.QPushButton("Save Session")
        save_session_button.clicked.connect(self.save_session)
        load_session_button = QtWidgets.QPushButton("Load Session")
        load_session_button.clicked.connect(self.load_session)
        session_layout.addWidget(save_session_button)
        session_layout.addWidget(load_session_button)
        main_layout_tab.addLayout(session_layout, 0, 0, 1, 3)

        # File selection area
        file_group = QtWidgets.QGroupBox("Data Files")
        file_layout = QtWidgets.QGridLayout()
        file_group.setLayout(file_layout)

        self.main_file_label = QtWidgets.QLabel("No file selected")
        main_file_button = QtWidgets.QPushButton("Select Main Data File")
        main_file_button.clicked.connect(self.select_main_data_file)
        file_layout.addWidget(main_file_button, 0, 0)
        file_layout.addWidget(self.main_file_label, 0, 1)

        self.csv_save_label = QtWidgets.QLabel("No CSV save path selected")
        csv_save_button = QtWidgets.QPushButton("Select CSV Save Path")
        csv_save_button.clicked.connect(self.select_csv_save_path)
        file_layout.addWidget(csv_save_button, 1, 0)
        file_layout.addWidget(self.csv_save_label, 1, 1)

        main_layout_tab.addWidget(file_group, 1, 0, 1, 3)

        # Data Range Selection area
        data_range_group = QtWidgets.QGroupBox("Data Range Selection")
        data_range_layout = QtWidgets.QVBoxLayout()
        data_range_group.setLayout(data_range_layout)

        self.data_range_all_radio = QtWidgets.QRadioButton("Use All Data")
        self.data_range_all_radio.setChecked(True)
        self.data_range_partial_radio = QtWidgets.QRadioButton("Enter Data Range")
        data_range_layout.addWidget(self.data_range_all_radio)
        data_range_layout.addWidget(self.data_range_partial_radio)

        day_hour_layout = QtWidgets.QHBoxLayout()
        self.start_day_hour_entry = QtWidgets.QLineEdit()
        self.end_day_hour_entry = QtWidgets.QLineEdit()
        self.start_day_hour_entry.setPlaceholderText("Start (Day:Hour)")
        self.end_day_hour_entry.setPlaceholderText("End (Day:Hour)")
        self.start_day_hour_entry.setEnabled(False)
        self.end_day_hour_entry.setEnabled(False)

        day_hour_layout.addWidget(QtWidgets.QLabel("Start:"))
        day_hour_layout.addWidget(self.start_day_hour_entry)
        day_hour_layout.addWidget(QtWidgets.QLabel("End:"))
        day_hour_layout.addWidget(self.end_day_hour_entry)
        data_range_layout.addLayout(day_hour_layout)

        
        shift_button = QtWidgets.QPushButton("Shift Data by 12h")
        shift_button.clicked.connect(self.toggle_shift_12h)
        data_range_layout.addWidget(shift_button)

        self.data_range_partial_radio.toggled.connect(self.update_data_range_entry_state)
        main_layout_tab.addWidget(data_range_group, 2, 0, 1, 3)

        # Tau selection area
        tau_group = QtWidgets.QGroupBox("Free Running Period (Tau) Selection")
        tau_layout = QtWidgets.QVBoxLayout()
        tau_group.setLayout(tau_layout)

        self.tau_radio_calculated = QtWidgets.QRadioButton("Use Calculated Tau")
        self.tau_radio_calculated.setChecked(True)
        self.tau_radio_calculated.toggled.connect(self.update_tau_entry_state)
        self.tau_radio_manual = QtWidgets.QRadioButton("Enter Tau Manually")

        tau_layout.addWidget(self.tau_radio_calculated)
        tau_layout.addWidget(self.tau_radio_manual)

        self.tau_method_combo = QtWidgets.QComboBox()
        self.tau_method_combo.addItems(["Fourier Analysis", "Chi-Squared", "Lomb-Scargle"])
        tau_layout.addWidget(QtWidgets.QLabel("Select Tau Calculation Method:"))
        tau_layout.addWidget(self.tau_method_combo)

        self.manual_tau_entry = QtWidgets.QLineEdit()
        self.manual_tau_entry.setEnabled(False)
        tau_layout.addWidget(QtWidgets.QLabel("Enter Tau Value (hours):"))
        tau_layout.addWidget(self.manual_tau_entry)

        main_layout_tab.addWidget(tau_group, 3, 0, 1, 3)

        # Input fields for threshold and activity hours
        input_group = QtWidgets.QGroupBox("Activity Parameters")
        input_layout = QtWidgets.QGridLayout()
        input_group.setLayout(input_layout)

        self.threshold_entry = QtWidgets.QLineEdit(str(self.threshold_percentile))
        self.N_inactivity_entry = QtWidgets.QLineEdit(str(self.N_hours_inactivity))
        self.M_activity_entry = QtWidgets.QLineEdit(str(self.M_hours_activity))

        input_layout.addWidget(QtWidgets.QLabel("Threshold Percentile:"), 0, 0)
        input_layout.addWidget(self.threshold_entry, 0, 1)
        input_layout.addWidget(QtWidgets.QLabel("N Hours of Inactivity:"), 1, 0)
        input_layout.addWidget(self.N_inactivity_entry, 1, 1)
        input_layout.addWidget(QtWidgets.QLabel("M Hours of Activity:"), 2, 0)
        input_layout.addWidget(self.M_activity_entry, 2, 1)

        main_layout_tab.addWidget(input_group, 4, 0, 1, 3)

        # Task selection
        task_group = QtWidgets.QGroupBox("Select Tasks")
        task_layout = QtWidgets.QVBoxLayout()
        task_group.setLayout(task_layout)

        self.checkbox_activity_profile = QtWidgets.QCheckBox("Generate Activity Profile")
        task_layout.addWidget(self.checkbox_activity_profile)

        self.checkbox_single_actogram = QtWidgets.QCheckBox("Generate Single Actogram")
        self.checkbox_single_actogram.stateChanged.connect(self.update_label_options_state)
        self.checkbox_double_actogram = QtWidgets.QCheckBox("Generate Double Actogram")
        self.checkbox_double_actogram.stateChanged.connect(self.update_label_options_state)
        self.checkbox_plot_fourier = QtWidgets.QCheckBox("Generate Periodogram")
        self.checkbox_save_csv = QtWidgets.QCheckBox("Save to CSV")

        task_layout.addWidget(self.checkbox_single_actogram)
        task_layout.addWidget(self.checkbox_double_actogram)
        task_layout.addWidget(self.checkbox_plot_fourier)
        task_layout.addWidget(self.checkbox_save_csv)

        main_layout_tab.addWidget(task_group, 5, 0, 3, 1)

        # Label options and Activity Profile Options
        self.label_group = QtWidgets.QGroupBox("Label Options")
        self.label_group.setVisible(False)
        label_layout = QtWidgets.QVBoxLayout()
        self.label_group.setLayout(label_layout)

        self.checkbox_label_onset = QtWidgets.QCheckBox("Label Onset")
        self.checkbox_label_offset = QtWidgets.QCheckBox("Label Offset")
        self.checkbox_label_acrophase = QtWidgets.QCheckBox("Label Acrophase")
        self.checkbox_label_bathyphase = QtWidgets.QCheckBox("Label Bathyphase")

        label_layout.addWidget(self.checkbox_label_onset)
        label_layout.addWidget(self.checkbox_label_offset)
        label_layout.addWidget(self.checkbox_label_acrophase)
        label_layout.addWidget(self.checkbox_label_bathyphase)

        self.activity_profile_options_group = QtWidgets.QGroupBox("Activity Profile Options")
        self.activity_profile_options_layout = QtWidgets.QVBoxLayout()
        self.activity_profile_options_group.setLayout(self.activity_profile_options_layout)

        self.checkbox_display_com = QtWidgets.QCheckBox("Display Center of Mass Line")
        self.checkbox_display_sem = QtWidgets.QCheckBox("Display SEM")
        self.checkbox_base_on_tau = QtWidgets.QCheckBox("Base Activity Profile on Tau (Unselect for 24 Hour Profile)")
        self.checkbox_base_on_tau.setChecked(True)

        self.activity_profile_options_layout.addWidget(self.checkbox_display_com)
        self.activity_profile_options_layout.addWidget(self.checkbox_display_sem)
        self.activity_profile_options_layout.addWidget(self.checkbox_base_on_tau)

        self.activity_profile_options_group.setVisible(False)
        main_layout_tab.addWidget(self.activity_profile_options_group, 5, 1, 3, 1)
        main_layout_tab.addWidget(self.label_group, 5, 2, 3, 1)

        self.checkbox_activity_profile.stateChanged.connect(self.toggle_shift_12h)

        # Run Analysis button
        run_button = QtWidgets.QPushButton("Run Analysis")
        run_button.clicked.connect(self.run_analysis)
        main_layout_tab.addWidget(run_button, 8, 0, 1, 3)

        # Plots tab
        plots_tab = QtWidgets.QWidget()
        tabs.addTab(plots_tab, "Plots")
        plots_layout = QtWidgets.QVBoxLayout(plots_tab)

        self.plots_tab_widget = QtWidgets.QTabWidget()
        plots_layout.addWidget(self.plots_tab_widget)

        export_markers_button = QtWidgets.QPushButton("Export Markers to CSV")
        export_markers_button.clicked.connect(self.export_markers_to_csv)
        plots_layout.addWidget(export_markers_button)

        # Results tab
        results_tab = QtWidgets.QWidget()
        tabs.addTab(results_tab, "Results")
        results_layout = QtWidgets.QVBoxLayout(results_tab)

        self.results_text = QtWidgets.QTextEdit()
        self.results_text.setReadOnly(True)
        results_layout.addWidget(self.results_text)

        # Status Bar
        self.status_bar = QtWidgets.QStatusBar()
        self.setStatusBar(self.status_bar)

    def toggle_shift_12h(self):
        """
        Toggle whether we shift the data by 12 hours before analysis.
        This simply sets a boolean flag (self.shift_12h); the actual shift is done in run_analysis().
        """
        self.shift_12h = not self.shift_12h
        if self.shift_12h:
            QtWidgets.QMessageBox.information(
                self,
                "12h Shift Enabled",
                "Data will be shifted by 12 hours when you run the analysis."
            )
        else:
            QtWidgets.QMessageBox.information(
                self,
                "12h Shift Disabled",
                "12-hour shift is now disabled."
            )

    def select_main_data_file(self):
        """
        Show a file dialog to select the main CSV data file and update self.file_path accordingly.
        """
        options = QtWidgets.QFileDialog.Options()
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select Main Data CSV File",
            "",
            "CSV Files (*.csv)",
            options=options
        )
        if file_path:
            self.file_path = file_path
            self.main_file_label.setText(self.file_path)

    def select_csv_save_path(self):
        """
        Show a file dialog to select or create a CSV file path where output data can be saved.
        """
        options = QtWidgets.QFileDialog.Options()
        default_filename = "mean_values.csv"
        file_path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Select CSV Save Path",
            default_filename,
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )
        if file_path:
            self.csv_save_path = file_path
            self.csv_save_label.setText(self.csv_save_path)

    def update_tau_entry_state(self):
        """
        Enable or disable the Tau entry fields depending on whether
        'Use Calculated Tau' or 'Enter Tau Manually' is selected.
        """
        if self.tau_radio_calculated.isChecked():
            self.tau_method_combo.setEnabled(True)
            self.manual_tau_entry.setEnabled(False)
        else:
            self.tau_method_combo.setEnabled(False)
            self.manual_tau_entry.setEnabled(True)

    def update_label_options_state(self):
        """
        Show or hide the label options group depending on whether
        single or double actogram generation is requested.
        """
        if self.checkbox_single_actogram.isChecked() or self.checkbox_double_actogram.isChecked():
            self.label_group.setVisible(True)
        else:
            self.label_group.setVisible(False)

    def add_plot(self, figure, title, markers=[]):
        """
        Add a new plot (figure) to the 'Plots' tab as a new tab.

        Parameters:
            figure (matplotlib.figure.Figure): The figure to display.
            title (str): The name of the tab.
            markers (list): A list of (marker, day, label_type) for any draggable markers present.
        """
        plot_widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(plot_widget)
        canvas = FigureCanvas(figure)
        toolbar = NavigationToolbar(canvas, self)
        layout.addWidget(toolbar)
        layout.addWidget(canvas)
        self.plots_tab_widget.addTab(plot_widget, title)

        for marker, day, label_type in markers:
            key = (day, label_type)
            if key not in self.marker_dict:
                self.marker_dict[key] = []
            self.marker_dict[key].append(marker)

    def on_marker_drag(self, day, label_type, new_time):
        """
        Callback for when a draggable marker is moved. Updates the underlying activity data.

        Parameters:
            day (int): The index of the day/period in which the marker is present.
            label_type (str): 'onset', 'offset', 'acrophase', or 'bathyphase'.
            new_time (float): The new (hour) position for the marker.
        """
        activity = self.activity
        if label_type == 'onset':
            new_bin = int((new_time % 24) * 60 / activity.bin_size_minutes)
            activity.onset_times[day] = new_bin
        elif label_type == 'offset':
            new_bin = int((new_time % 24) * 60 / activity.bin_size_minutes)
            activity.offset_times[day] = new_bin
        elif label_type == 'acrophase':
            hours = int(new_time % 24)
            minutes = int((new_time % 24 - hours) * 60)
            activity.acrophase_times[day] = (hours, minutes)
        elif label_type == 'bathyphase':
            hours = int(new_time % 24)
            minutes = int((new_time % 24 - hours) * 60)
            activity.bathophase_times[day] = (hours, minutes)

        key = (day, label_type)
        if key in self.marker_dict:
            for marker in self.marker_dict[key]:
                marker.set_data([new_time % 24], [day + 0.5])
                marker.figure.canvas.draw_idle()

    def update_data_range_entry_state(self):
        """
        Enable or disable the Start/End data range fields depending on the radio selection.
        """
        if self.data_range_partial_radio.isChecked():
            self.start_day_hour_entry.setEnabled(True)
            self.end_day_hour_entry.setEnabled(True)
        else:
            self.start_day_hour_entry.setEnabled(False)
            self.end_day_hour_entry.setEnabled(False)

    def run_analysis(self):
        """
        Main entry point for performing the circadian activity analysis.
        Loads data, optionally shifts by 12h, calculates tau, generates plots, etc.
        """
        self.results_text.clear()

        # Clear existing plots
        while self.plots_tab_widget.count():
            self.plots_tab_widget.removeTab(0)
        self.marker_dict = {}

        if not self.file_path:
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                "Please select the main data file before running the analysis."
            )
            return

        # Load data
        try:
            data = pd.read_csv(self.file_path, header=None)
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                f"Error loading main data file: {e}"
            )
            return

        # Extract basic info
        try:
            start_hour_str = data.iloc[2, 0]
            bin_size_minutes = int(float(data.iloc[3, 0])) // 4
            event_type = data.iloc[0, 0]
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                f"Error extracting data from main file: {e}"
            )
            return

        # Shift data by 12h if requested
        if self.shift_12h:
            try:
                hour_part, minute_part = start_hour_str.split(':')
                hour_val = int(hour_part)
                min_val = int(minute_part)
                hour_val = (hour_val + 12) % 24
                start_hour_str = f"{hour_val:02d}:{min_val:02d}"
            except:
                QtWidgets.QMessageBox.warning(
                    self,
                    "Shift Warning",
                    "Could not parse start hour string. 12h shift skipped."
                )

        # Convert the event counts
        raw_event_counts = data.iloc[7:, 0]
        event_counts = pd.to_numeric(raw_event_counts, errors='coerce').fillna(0).astype(int)

        # Handle partial range if selected
        if self.data_range_partial_radio.isChecked():
            start_str = self.start_day_hour_entry.text().strip()
            end_str = self.end_day_hour_entry.text().strip()

            def parse_day_hour(s):
                d, h = s.split(':')
                d = int(d)
                h = int(h)
                return d * 24 + h

            start_hour_total = parse_day_hour(start_str)
            end_hour_total = parse_day_hour(end_str)
            start_bin = int(start_hour_total * 60 / bin_size_minutes)
            end_bin = int(end_hour_total * 60 / bin_size_minutes)
            event_counts = event_counts.iloc[start_bin:end_bin].reset_index(drop=True)

        L = len(event_counts)
        data_modified = data.iloc[:7 + L, :].copy()

        # Update the start hour if we shifted
        data_modified.iloc[2, 0] = start_hour_str
        data_modified.iloc[7:, 0] = event_counts.values

        # Initialize ActivityAnalysis
        self.activity = ActivityAnalysis(
            data_modified,
            start_hour_str,
            bin_size_minutes,
            event_type
        )
        activity = self.activity

        # Determine tau
        self.plot_fourier = self.checkbox_plot_fourier.isChecked()
        if self.tau_radio_calculated.isChecked():
            selected_method = self.tau_method_combo.currentText()
            try:
                if selected_method == "Fourier Analysis":
                    if self.plot_fourier:
                        tau, fig = activity.fourier(plot=True)
                        self.add_plot(fig, "Fourier Analysis")
                    else:
                        tau = activity.fourier(plot=False)
                elif selected_method == "Chi-Squared":
                    if self.plot_fourier:
                        tau, fig = activity.chi_squared(plot=True)
                        self.add_plot(fig, "Chi-Squared Analysis")
                    else:
                        tau = activity.chi_squared(plot=False)
                elif selected_method == "Lomb-Scargle":
                    if self.plot_fourier:
                        tau, fig = activity.lomb_scargle(plot=True)
                        self.add_plot(fig, "Lomb-Scargle Analysis")
                    else:
                        tau = activity.lomb_scargle(plot=False)
                else:
                    QtWidgets.QMessageBox.critical(
                        self,
                        "Error",
                        f"Unknown tau calculation method: {selected_method}"
                    )
                    return
                self.results_text.append(f"Using calculated free running period: {tau:.2f} hours\n")
            except Exception as e:
                QtWidgets.QMessageBox.critical(
                    self,
                    "Error",
                    f"Error in {selected_method} tau calculation: {e}"
                )
                return
        else:
            try:
                tau = float(self.manual_tau_entry.text())
                activity.tau = tau
                self.results_text.append(f"Using manually entered free running period: {tau:.2f} hours\n")
                if self.plot_fourier:
                    _, fig = activity.fourier(plot=True)
                    self.add_plot(fig, "Fourier Analysis")
            except ValueError:
                QtWidgets.QMessageBox.critical(
                    self,
                    "Error",
                    "Invalid tau value. Please enter a numerical value."
                )
                return
            except Exception as e:
                QtWidgets.QMessageBox.critical(
                    self,
                    "Error",
                    f"Error in Fourier analysis: {e}"
                )
                return

        activity.tau = tau
        activity.bins_per_period = int((tau * 60) / activity.bin_size_minutes)

        # Get user inputs for threshold, inactivity, activity
        try:
            self.threshold_percentile = float(self.threshold_entry.text())
            self.N_hours_inactivity = float(self.N_inactivity_entry.text())
            self.M_hours_activity = float(self.M_activity_entry.text())
        except ValueError:
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                "Invalid input for threshold or N/M hours. Please enter numerical values."
            )
            return

        # Retrieve activity profile options
        self.display_com = self.checkbox_display_com.isChecked()
        self.display_sem = self.checkbox_display_sem.isChecked()
        self.base_on_tau = self.checkbox_base_on_tau.isChecked()

        # Perform analysis steps
        binary_activity = activity.calculate_binary_activity(self.threshold_percentile)
        convolution_result = activity.perform_convolution(
            binary_activity,
            N_hours_inactivity=self.N_hours_inactivity,
            M_hours_activity=self.M_hours_activity
        )
        activity.detect_onset_offset(convolution_result)
        activity.detect_onset_offset_res(convolution_result)
        activity.acro_bathy()
        onset_av, offset_av = activity.calculate_mean_on_off(
            activity.onset_times,
            activity.offset_times
        )
        df = activity.results_writer(
            activity.onset_res,
            activity.offset_res,
            activity.acro_res,
            activity.batho_res,
            tau,
            onset_av,
            offset_av
        )

        # Handle actogram and activity profile checkboxes
        self.label_onset = self.checkbox_label_onset.isChecked()
        self.label_offset = self.checkbox_label_offset.isChecked()
        self.label_acrophase = self.checkbox_label_acrophase.isChecked()
        self.label_bathyphase = self.checkbox_label_bathyphase.isChecked()

        if self.checkbox_single_actogram.isChecked():
            fig, markers = activity.single_actogram(
                label_onset=self.label_onset,
                label_offset=self.label_offset,
                label_acrophase=self.label_acrophase,
                label_bathyphase=self.label_bathyphase,
                on_drag=self.on_marker_drag,
                on_release_callback=self.regenerate_single_actogram
            )
            self.add_plot(fig, "Single Actogram", markers)

        if self.checkbox_double_actogram.isChecked():
            fig, markers = activity.double_actogram(
                label_onset=self.label_onset,
                label_offset=self.label_offset,
                label_acrophase=self.label_acrophase,
                label_bathyphase=self.label_bathyphase,
                on_drag=self.on_marker_drag,
                on_release_callback=self.regenerate_double_actogram
            )
            self.add_plot(fig, "Double Actogram", markers)

        if self.checkbox_activity_profile.isChecked():
            try:
                fig = activity.plot_activity_profile(
                    display_com=self.display_com,
                    display_sem=self.display_sem,
                    base_on_tau=self.base_on_tau
                )
                self.add_plot(fig, "Activity Profile")
            except Exception as e:
                QtWidgets.QMessageBox.critical(
                    self,
                    "Error",
                    f"Error generating Activity Profile: {e}"
                )
                return

        if self.checkbox_save_csv.isChecked():
            if not self.csv_save_path:
                QtWidgets.QMessageBox.warning(
                    self,
                    "CSV Save Path Not Set",
                    "Please select a CSV save path using the 'Select CSV Save Path' button."
                )
                return
            activity.csv_save(df, self.csv_save_path)

        # Display results
        self.results_text.append("\nOnset Times:")
        self.results_text.append('\n'.join(activity.onset_res))
        self.results_text.append("\n\nOffset Times:")
        self.results_text.append('\n'.join(activity.offset_res))
        self.results_text.append("\n\nAcrophase Times:")
        self.results_text.append('\n'.join(activity.acro_res))
        self.results_text.append("\n\nBathyphase Times:")
        self.results_text.append('\n'.join(activity.batho_res))
        self.results_text.append(f"\n\nDominant Period (Tau): {tau:.2f} hours")
        self.results_text.append(f"Onset Average: {onset_av}")
        self.results_text.append(f"Offset Average: {offset_av}")

        self.status_bar.showMessage("Analysis complete.", 5000)

    def regenerate_single_actogram(self):
        """
        Re-generate the Single Actogram plot after markers have been dragged.
        """
        try:
            for i in range(self.plots_tab_widget.count()):
                if self.plots_tab_widget.tabText(i) == "Single Actogram":
                    self.plots_tab_widget.removeTab(i)
                    break

            fig, markers = self.activity.single_actogram(
                label_onset=self.label_onset,
                label_offset=self.label_offset,
                label_acrophase=self.label_acrophase,
                label_bathyphase=self.label_bathyphase,
                on_drag=self.on_marker_drag,
                on_release_callback=self.regenerate_single_actogram
            )
            self.add_plot(fig, "Single Actogram", markers)
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed to regenerate Single Actogram: {e}")

    def regenerate_double_actogram(self):
        """
        Re-generate the Double Actogram plot after markers have been dragged.
        """
        try:
            fig, markers = self.activity.double_actogram(
                label_onset=self.label_onset,
                label_offset=self.label_offset,
                label_acrophase=self.label_acrophase,
                label_bathyphase=self.label_bathyphase,
                on_drag=self.on_marker_drag,
                on_release_callback=self.regenerate_double_actogram
            )
            self.add_plot(fig, "Double Actogram", markers)
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed to regenerate Double Actogram: {e}")

        for i in range(self.plots_tab_widget.count()):
            if self.plots_tab_widget.tabText(i) == "Double Actogram":
                self.plots_tab_widget.removeTab(i)
                break

    def save_session(self):
        """
        Save the current session (file paths, parameters, marker positions, etc.) to a JSON file.
        """
        if not hasattr(self, 'activity') or self.activity is None:
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                "No analysis has been run yet. Please run an analysis before saving the session."
            )
            return

        session_data = {
            'file_path': self.file_path,
            'csv_save_path': self.csv_save_path,
            'threshold_percentile': self.threshold_percentile,
            'N_hours_inactivity': self.N_hours_inactivity,
            'M_hours_activity': self.M_hours_activity,
            'tau_radio_calculated': self.tau_radio_calculated.isChecked(),
            'tau_calculation_method': (
                self.tau_method_combo.currentText() if self.tau_radio_calculated.isChecked() else None
            ),
            'manual_tau': self.manual_tau_entry.text(),
            'checkboxes': {
                'single_actogram': self.checkbox_single_actogram.isChecked(),
                'double_actogram': self.checkbox_double_actogram.isChecked(),
                'plot_fourier': self.checkbox_plot_fourier.isChecked(),
                'save_csv': self.checkbox_save_csv.isChecked(),
                'activity_profile': self.checkbox_activity_profile.isChecked(),
                'label_onset': self.checkbox_label_onset.isChecked(),
                'label_offset': self.checkbox_label_offset.isChecked(),
                'label_acrophase': self.checkbox_label_acrophase.isChecked(),
                'label_bathyphase': self.checkbox_label_bathyphase.isChecked(),
                'display_com': self.checkbox_display_com.isChecked(),
                'display_sem': self.checkbox_display_sem.isChecked(),
                'base_on_tau': self.checkbox_base_on_tau.isChecked(),
            },
            'marker_positions': {
                'onset_times': self.activity.onset_times,
                'offset_times': self.activity.offset_times,
                'acrophase_times': self.activity.acrophase_times,
                'bathophase_times': self.activity.bathophase_times,
            }
        }

        options = QtWidgets.QFileDialog.Options()
        file_path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save Session",
            "",
            "JSON Files (*.json)",
            options=options
        )
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    json.dump(session_data, f, indent=4)
                QtWidgets.QMessageBox.information(
                    self,
                    "Session Saved",
                    f"Session saved to {file_path}"
                )
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, "Error", f"Failed to save session: {e}")

    def load_session(self):
        """
        Load a previously saved session from a JSON file, restoring file paths, parameters, and markers.
        """
        options = QtWidgets.QFileDialog.Options()
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Load Session",
            "",
            "JSON Files (*.json)",
            options=options
        )
        if file_path:
            try:
                with open(file_path, 'r') as f:
                    session_data = json.load(f)

                self.file_path = session_data.get('file_path')
                self.main_file_label.setText(self.file_path if self.file_path else "No file selected")
                self.csv_save_path = session_data.get('csv_save_path')
                self.csv_save_label.setText(self.csv_save_path if self.csv_save_path else "No CSV save path selected")

                self.threshold_percentile = session_data.get('threshold_percentile', 30.0)
                self.N_hours_inactivity = session_data.get('N_hours_inactivity', 12.0)
                self.M_hours_activity = session_data.get('M_hours_activity', 12.0)
                self.threshold_entry.setText(str(self.threshold_percentile))
                self.N_inactivity_entry.setText(str(self.N_hours_inactivity))
                self.M_activity_entry.setText(str(self.M_hours_activity))

                self.tau_radio_calculated.setChecked(session_data.get('tau_radio_calculated', True))
                if self.tau_radio_calculated.isChecked():
                    self.tau_method_combo.setCurrentText(
                        session_data.get('tau_calculation_method', "Fourier Analysis")
                    )
                else:
                    self.tau_method_combo.setCurrentIndex(-1)

                self.manual_tau_entry.setText(session_data.get('manual_tau', ''))
                self.update_tau_entry_state()

                checkboxes = session_data.get('checkboxes', {})
                self.checkbox_single_actogram.setChecked(checkboxes.get('single_actogram', False))
                self.checkbox_double_actogram.setChecked(checkboxes.get('double_actogram', False))
                self.checkbox_plot_fourier.setChecked(checkboxes.get('plot_fourier', False))
                self.checkbox_save_csv.setChecked(checkboxes.get('save_csv', False))
                self.checkbox_activity_profile.setChecked(checkboxes.get('activity_profile', False))
                self.checkbox_label_onset.setChecked(checkboxes.get('label_onset', False))
                self.checkbox_label_offset.setChecked(checkboxes.get('label_offset', False))
                self.checkbox_label_acrophase.setChecked(checkboxes.get('label_acrophase', False))
                self.checkbox_label_bathyphase.setChecked(checkboxes.get('label_bathyphase', False))
                self.checkbox_display_com.setChecked(checkboxes.get('display_com', False))
                self.checkbox_display_sem.setChecked(checkboxes.get('display_sem', False))
                self.checkbox_base_on_tau.setChecked(checkboxes.get('base_on_tau', True))
                self.update_label_options_state()

                marker_positions = session_data.get('marker_positions', {})
                if not hasattr(self, 'activity') or self.activity is None:
                    if self.file_path:
                        try:
                            data = pd.read_csv(self.file_path, header=None)
                        except Exception as e:
                            QtWidgets.QMessageBox.critical(
                                self,
                                "Error",
                                f"Error loading main data file: {e}"
                            )
                            return

                        try:
                            start_hour = data.iloc[2, 0]
                            bin_size_minutes = int(float(data.iloc[3, 0])) // 4
                            event_type = data.iloc[0, 0]
                        except Exception as e:
                            QtWidgets.QMessageBox.critical(
                                self,
                                "Error",
                                f"Error extracting data from main file: {e}"
                            )
                            return

                        self.activity = ActivityAnalysis(data, start_hour, bin_size_minutes, event_type)

                self.activity.onset_times = marker_positions.get('onset_times', [])
                self.activity.offset_times = marker_positions.get('offset_times', [])
                self.activity.acrophase_times = marker_positions.get('acrophase_times', [])
                self.activity.bathophase_times = marker_positions.get('bathophase_times', [])

                self.run_analysis()
                QtWidgets.QMessageBox.information(
                    self,
                    "Session Loaded",
                    f"Session loaded from {file_path}"
                )
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, "Error", f"Failed to load session: {e}")

    def export_markers_to_csv(self):
        """
        Export current (and updated) onset, offset, acrophase, and bathyphase markers to CSV.
        Also appends a summary section at the end with average onset, offset, and the dominant period (tau).
        """
        if not hasattr(self, 'activity') or self.activity is None:
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                "No analysis has been run yet. Please run an analysis before exporting markers."
            )
            return

        options = QtWidgets.QFileDialog.Options()
        default_filename = "updated_markers.csv"
        file_path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export Markers to CSV",
            default_filename,
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )
        if not file_path:
            return

        num_periods = max(
            len(self.activity.onset_times),
            len(self.activity.offset_times),
            len(self.activity.acrophase_times),
            len(self.activity.bathophase_times)
        )
        
        onset_list = []
        offset_list = []
        acrophase_list = []
        bathyphase_list = []
        
        for i in range(num_periods):
            # Onset
            if i < len(self.activity.onset_times):
                onset_bin = self.activity.onset_times[i]
                onset_hour = (onset_bin * self.activity.bin_size_minutes) // 60
                onset_minute = (onset_bin * self.activity.bin_size_minutes) % 60
                onset_str = f"{onset_hour:02d}:{onset_minute:02d}"
            else:
                onset_str = ""
            onset_list.append(onset_str)

            # Offset
            if i < len(self.activity.offset_times):
                offset_bin = self.activity.offset_times[i]
                offset_hour = (offset_bin * self.activity.bin_size_minutes) // 60
                offset_minute = (offset_bin * self.activity.bin_size_minutes) % 60
                offset_str = f"{offset_hour:02d}:{offset_minute:02d}"
            else:
                offset_str = ""
            offset_list.append(offset_str)

            # Acrophase
            if i < len(self.activity.acrophase_times):
                acro = self.activity.acrophase_times[i]
                if isinstance(acro, (list, tuple)) and len(acro) == 2:
                    acro_str = f"{acro[0]:02d}:{acro[1]:02d}"
                else:
                    acro_str = str(acro)
            else:
                acro_str = ""
            acrophase_list.append(acro_str)

            # Bathyphase
            if i < len(self.activity.bathophase_times):
                batho = self.activity.bathophase_times[i]
                if isinstance(batho, (list, tuple)) and len(batho) == 2:
                    batho_str = f"{batho[0]:02d}:{batho[1]:02d}"
                else:
                    batho_str = str(batho)
            else:
                batho_str = ""
            bathyphase_list.append(batho_str)

        df = pd.DataFrame({
            "Period": [f"Period {i+1}" for i in range(num_periods)],
            "Onset": onset_list,
            "Offset": offset_list,
            "Acrophase": acrophase_list,
            "Bathyphase": bathyphase_list
        })

        onset_av, offset_av = self.activity.calculate_mean_on_off(
            self.activity.onset_times,
            self.activity.offset_times
        )
        tau_str = f"{self.activity.tau:.2f}" if self.activity.tau is not None else ""
        summary_data = {
            "Metric": ["Onset Average", "Offset Average", "Dominant Period (Tau)"],
            "Value": [onset_av, offset_av, tau_str]
        }
        df_summary = pd.DataFrame(summary_data)

        try:
            with open(file_path, 'w', newline='') as f:
                df.to_csv(f, index=False)
                f.write("\n")
                df_summary.to_csv(f, index=False)
            QtWidgets.QMessageBox.information(
                self,
                "Export Successful",
                f"Markers exported to {file_path}"
            )
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self,
                "Export Failed",
                f"Failed to export markers: {e}"
            )

    def convert_awd_file(self):
        """
        Open a file dialog to select an AWD file, then convert it to CSV
        using the convert_awd_to_csv function.
        """
        options = QtWidgets.QFileDialog.Options()
        awd_file, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select AWD File",
            "",
            "AWD Files (*.awd);;All Files (*)",
            options=options
        )
        if not awd_file:
            return

        default_csv_name = os.path.splitext(os.path.basename(awd_file))[0] + ".csv"
        csv_file, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save Converted CSV File",
            default_csv_name,
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )
        if not csv_file:
            return

        try:
            convert_awd_to_csv(awd_file, csv_file)
            QtWidgets.QMessageBox.information(
                self,
                "Conversion Successful",
                f"Converted file saved as {csv_file}"
            )
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self,
                "Conversion Failed",
                f"Failed to convert file: {e}"
            )

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    gui = ActivityAnalysisGUI()
    gui.show()
    sys.exit(app.exec_())
