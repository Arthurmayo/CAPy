import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from scipy.signal import convolve, lombscargle
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import os

class ActivityAnalysis:
    """
    A class that processes circadian activity data and generates actograms
    with an overlay for light data.
    """

    def __init__(self, data, start_hour, bin_size_minutes, event_type):
        """
        Initialize the ActivityAnalysis class with necessary data and parameters.

        Parameters:
            data (pd.DataFrame): The loaded CSV data.
            start_hour (str): The start hour in 'HH:MM' format.
            bin_size_minutes (int): The size of each bin in minutes.
            event_type (str): The type of event being analyzed.
        
        Raises:
            ValueError: If parsing the start datetime fails.
        """
        self.data = data
        self.event_type = event_type
        self.start_hour = start_hour
        self.bin_size_minutes = bin_size_minutes

        try:
            self.start_datetime = datetime.strptime(
                f"{self.data.iloc[1, 0]} {self.data.iloc[2, 0]}",
                "%d-%b-%Y %H:%M"
            )
        except Exception as e:
            raise ValueError(f"Error parsing start datetime: {e}")

        # Process event counts from row 8 onward (column 0)
        raw_event_counts = self.data.iloc[7:, 0]
        self.event_counts = pd.to_numeric(raw_event_counts, errors='coerce')
        if self.event_counts.isnull().any():
            print("Warning: Some event counts could not be converted; setting them to 0.")
            self.event_counts = self.event_counts.fillna(0).astype(int)
        else:
            self.event_counts = self.event_counts.astype(int)

        # Process light values from row 8 onward (column 1); if missing, use zeros.
        if self.data.shape[1] > 1:
            raw_light = self.data.iloc[7:, 1]
            self.light_values = pd.to_numeric(raw_light, errors='coerce').fillna(0).astype(int)
        else:
            self.light_values = np.zeros_like(self.event_counts)

        # Determine bins per day.
        self.bins_per_day = (60 // bin_size_minutes) * 24
        start_hour_int, start_minute_int = map(int, start_hour.split(':'))
        start_minutes = start_hour_int * 60 + start_minute_int
        bins_to_pad_start = start_minutes // bin_size_minutes

        # Pad event counts.
        self.padded_event_counts = np.pad(self.event_counts, (bins_to_pad_start, 0), 'constant')
        self.num_days = int(np.ceil(len(self.padded_event_counts) / self.bins_per_day))
        self.padded_event_counts = np.pad(
            self.padded_event_counts,
            (0, self.num_days * self.bins_per_day - len(self.padded_event_counts)),
            'constant'
        )

        # Pad light values in the same way.
        self.padded_light = np.pad(self.light_values, (bins_to_pad_start, 0), 'constant')
        self.padded_light = np.pad(
            self.padded_light,
            (0, self.num_days * self.bins_per_day - len(self.padded_light)),
            'constant'
        )

        # Initialize placeholders for results.
        self.onset_times = []
        self.offset_times = []
        self.acrophase_times = []
        self.bathophase_times = []
        self.onset_res = []
        self.offset_res = []
        self.acro_res = []
        self.batho_res = []
        self.draggable_markers = []
        self.tau = None
        self.bins_per_period = None

    def lomb_scargle(self, plot=False):
        """
        Perform Lomb-Scargle periodogram analysis to estimate tau.
        
        Parameters:
            plot (bool): If True, return a figure of the periodogram.
        
        Returns:
            float or (float, matplotlib.figure.Figure): The estimated tau (in hours).
        """
        min_period = 10
        max_period = 36
        frequencies = np.linspace(1 / max_period, 1 / min_period, 1000)
        angular_frequencies = 2 * np.pi * frequencies
        time = np.arange(len(self.event_counts)) * self.bin_size_minutes / 60.0
        y = self.event_counts - np.mean(self.event_counts)
        power = lombscargle(time, y, angular_frequencies)
        amplitude = np.sqrt(power)
        best_frequency = frequencies[np.argmax(amplitude)]
        tau = 1 / best_frequency
        self.tau = tau
        best_amplitude = amplitude.max()
        if plot:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(1 / frequencies, amplitude, 'k-', label='Lomb-Scargle Amplitude')
            ax.axvline(x=tau, color='r', linestyle='--', 
                       label=f'Dominant Tau: {tau:.2f} h\nAmplitude: {best_amplitude:.2f}')
            ax.set_xlabel('Tau (hours)')
            ax.set_ylabel('Amplitude')
            ax.set_title('Lomb-Scargle Periodogram')
            ax.legend()
            ax.grid(True)
            plt.tight_layout()
            return tau, fig
        else:
            return tau

    def fourier(self, plot=False):
        """
        Perform Fourier analysis (FFT) to estimate tau.
        
        Parameters:
            plot (bool): If True, return a figure of the amplitude spectrum.
        
        Returns:
            float or (float, matplotlib.figure.Figure): The estimated tau (in hours).
        """
        fft_values = np.fft.fft(self.event_counts)
        frequencies = np.fft.fftfreq(len(self.event_counts), d=self.bin_size_minutes * 60)
        amplitude_spectrum = np.abs(fft_values)
        nonzero = frequencies != 0
        frequencies = frequencies[nonzero]
        amplitude_spectrum = amplitude_spectrum[nonzero]
        period_in_seconds = 1 / frequencies
        period_in_hours = period_in_seconds / 3600
        valid = (period_in_hours >= 10) & (period_in_hours <= 36)
        period_in_hours = period_in_hours[valid]
        amplitude_spectrum = amplitude_spectrum[valid]
        if len(period_in_hours) == 0:
            raise ValueError("No valid periods found in Fourier analysis.")
        dominant_index = np.argmax(amplitude_spectrum)
        tau = period_in_hours[dominant_index]
        self.tau = tau
        best_amp = amplitude_spectrum[dominant_index]
        if plot:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(period_in_hours, amplitude_spectrum, 'k-', label='Amplitude Spectrum')
            ax.axvline(x=tau, color='r', linestyle='--', 
                       label=f'Dominant Tau: {tau:.2f} h\nAmplitude: {best_amp:.2f}')
            ax.set_xlabel('Period (hours)')
            ax.set_ylabel('Amplitude')
            ax.set_title('Fourier Amplitude Spectrum')
            ax.legend()
            ax.grid(True)
            plt.tight_layout()
            return tau, fig
        else:
            return tau

    def chi_squared(self, plot=False):
        """
        Perform Chi-Squared periodogram analysis to estimate tau.
        
        Parameters:
            plot (bool): If True, return a figure of the periodogram.
        
        Returns:
            float or (float, matplotlib.figure.Figure): The estimated tau (in hours).
        """
        possible_taus = np.linspace(10, 36, 260)
        chi_sq_values = []
        time = np.arange(len(self.event_counts)) * (self.bin_size_minutes / 60.0)
        for tau_candidate in possible_taus:
            bins_per_cycle = int(np.round(tau_candidate * 60 / self.bin_size_minutes))
            if bins_per_cycle <= 1:
                chi_sq_values.append(np.nan)
                continue
            num_cycles = len(self.event_counts) // bins_per_cycle
            if num_cycles < 1:
                chi_sq_values.append(np.nan)
                continue
            truncated = self.event_counts[:num_cycles * bins_per_cycle].to_numpy()
            folded = truncated.reshape((num_cycles, bins_per_cycle))
            observed = np.mean(folded, axis=0)
            expected = np.full_like(observed, np.mean(observed))
            chi, _ = chisquare(f_obs=observed, f_exp=expected)
            chi_sq_values.append(chi)
        chi_sq_values = np.array(chi_sq_values)
        max_idx = np.nanargmax(chi_sq_values)
        tau = possible_taus[max_idx]
        self.tau = tau
        if plot:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(possible_taus, chi_sq_values, 'k-', label='Chi-Squared')
            ax.set_xlabel('Tau (hours)')
            ax.set_ylabel('Chi-Squared')
            ax.set_title('Chi-Squared Periodogram')
            ax.axvline(tau, color='r', linestyle='--', label=f"Max Chi-Sq Tau: {tau:.2f} h")
            ax.legend()
            return tau, fig
        else:
            return tau

    def calculate_binary_activity(self, threshold_percentile):
        """
        Convert padded event counts into a binary array (-1, +1) based on a threshold.
        
        Parameters:
            threshold_percentile (float): The percentile threshold.
        
        Returns:
            np.ndarray: Binary activity array.
        """
        non_zero = self.padded_event_counts[self.padded_event_counts > 0]
        threshold = np.percentile(non_zero, threshold_percentile) if len(non_zero) > 0 else 0
        return np.where(self.padded_event_counts > threshold, 1, -1)

    def perform_convolution(self, binary_activity, N_hours_inactivity=12, M_hours_activity=12):
        """
        Convolve the binary activity with a template to aid in onset/offset detection.
        
        Parameters:
            binary_activity (np.ndarray): The binary activity array.
            N_hours_inactivity (float): Hours of inactivity (negative).
            M_hours_activity (float): Hours of activity (positive).
        
        Returns:
            np.ndarray: The convolved signal.
        """
        bins_per_hour = 60 // self.bin_size_minutes
        N_bins = int(N_hours_inactivity * bins_per_hour)
        M_bins = int(M_hours_activity * bins_per_hour)
        if N_bins <= 0 or M_bins <= 0:
            raise ValueError("N_bins and M_bins must be positive integers.")
        template = np.concatenate((np.full(N_bins, -1), np.full(M_bins, 1)))
        return convolve(binary_activity, template, mode='valid')

    def detect_onset_offset(self, convolution_result):
        """
        Detect onset and offset bin indices from the convolved signal.
        
        Parameters:
            convolution_result (np.ndarray): The convolved signal.
        
        Returns:
            (list, list): Onset and offset bin indices.
        """
        self.onset_times = []
        self.offset_times = []
        total_periods = len(convolution_result) // self.bins_per_period
        for period in range(total_periods):
            period_data = convolution_result[period * self.bins_per_period:(period + 1) * self.bins_per_period]
            self.onset_times.append(np.argmax(period_data))
            self.offset_times.append(np.argmin(period_data))
        return self.onset_times, self.offset_times

    def detect_onset_offset_res(self, convolution_result):
        """
        Convert onset/offset bin indices into "HH:MM" formatted strings.
        
        Parameters:
            convolution_result (np.ndarray): The convolved signal.
        
        Returns:
            (list, list): Onset and offset times as strings.
        """
        self.onset_res = []
        self.offset_res = []
        total_periods = len(self.onset_times)
        for period in range(total_periods):
            try:
                onset_bin = self.onset_times[period]
                offset_bin = self.offset_times[period]
                onset_hour = (onset_bin * self.bin_size_minutes) // 60
                onset_minute = (onset_bin * self.bin_size_minutes) % 60
                offset_hour = (offset_bin * self.bin_size_minutes) // 60
                offset_minute = (offset_bin * self.bin_size_minutes) % 60
                self.onset_res.append(f"{int(onset_hour):02d}:{int(onset_minute):02d}")
                self.offset_res.append(f"{int(offset_hour):02d}:{int(offset_minute):02d}")
            except IndexError:
                self.onset_res.append("00:00")
                self.offset_res.append("00:00")
        return self.onset_res, self.offset_res

    def acro_bathy(self):
        """
        Fit a sine curve to each period and determine acrophase and bathyphase.
        
        Returns:
            (list, list, list, list): Acrophase times, bathyphase times (as tuples) and their string representations.
        
        Raises:
            ValueError: If tau is not determined.
        """
        if self.tau is None:
            raise ValueError("Tau has not been determined yet.")
        time_axis = np.arange(len(self.event_counts)) * self.bin_size_minutes / 60.0

        def sine_func(t, A, phi, B):
            return A * np.sin(2 * np.pi * t / self.tau + phi) + B

        bins_per_period = self.bins_per_period
        self.acrophase_times = []
        self.bathophase_times = []
        self.acro_res = []
        self.batho_res = []
        total_periods = len(self.onset_times)
        for period in range(total_periods):
            try:
                start = period * bins_per_period
                end = (period + 1) * bins_per_period
                period_data = self.padded_event_counts[start:end]
                time_period = time_axis[start:end]
                if len(period_data) == 0:
                    raise ValueError("No data for this period.")
                initial_guess = [np.std(period_data), 0, np.mean(period_data)]
                params, _ = curve_fit(sine_func, time_period, period_data, p0=initial_guess)
                A_fit, phi_fit, B_fit = params
                fitted = sine_func(time_period, A_fit, phi_fit, B_fit)
                acrophase_bin = np.argmax(fitted)
                shift_bins = int((self.tau / 2) / (self.bin_size_minutes / 60))
                bathyphase_bin = (acrophase_bin + shift_bins) % int(bins_per_period)
                acrophase_hour = (acrophase_bin * self.bin_size_minutes) // 60
                acrophase_minute = (acrophase_bin * self.bin_size_minutes) % 60
                bathyphase_hour = (bathyphase_bin * self.bin_size_minutes) // 60
                bathyphase_minute = (bathyphase_bin * self.bin_size_minutes) % 60
                self.acrophase_times.append((int(acrophase_hour), int(acrophase_minute)))
                self.bathophase_times.append((int(bathyphase_hour), int(bathyphase_minute)))
            except Exception:
                self.acrophase_times.append((0, 0))
                self.bathophase_times.append((12, 0))
            self.acro_res.append(f"{self.acrophase_times[period][0]:02d}:{self.acrophase_times[period][1]:02d}")
            self.batho_res.append(f"{self.bathophase_times[period][0]:02d}:{self.bathophase_times[period][1]:02d}")
        return self.acrophase_times, self.bathophase_times, self.acro_res, self.batho_res

    def single_actogram(self, label_onset=False, label_offset=False,
                        label_acrophase=False, label_bathyphase=False,
                        on_drag=None, on_release_callback=None):
        """
        Generate and return a single-plotted actogram with a light-data overlay.
        
        Parameters:
            label_onset (bool): Whether to mark onsets.
            label_offset (bool): Whether to mark offsets.
            label_acrophase (bool): Whether to mark acrophase.
            label_bathyphase (bool): Whether to mark bathyphase.
            on_drag (callable): Callback for marker dragging.
            on_release_callback (callable): Callback after marker drag.
        
        Returns:
            (matplotlib.figure.Figure, list): The figure and a list of marker info.
        
        Raises:
            ValueError: If tau is not set.
        """
        if self.tau is None:
            raise ValueError("Tau has not been determined yet.")
        bins_per_period = self.bins_per_period
        events_per_period = [
            self.padded_event_counts[i * bins_per_period:(i + 1) * bins_per_period]
            for i in range(len(self.padded_event_counts) // bins_per_period)
        ]
        num_periods = len(events_per_period)
        fig_height = max(6, num_periods * 0.3)
        fig, ax = plt.subplots(figsize=(12, fig_height))
        cax = ax.imshow(events_per_period, aspect='auto', cmap='Greys', 
                        interpolation='none', extent=[0, self.tau, num_periods, 0])
        
        light_per_period = [
            self.padded_light[i * bins_per_period:(i + 1) * bins_per_period]
            for i in range(len(self.padded_light) // bins_per_period)
        ]
        light_cmap = LinearSegmentedColormap.from_list("light_cmap", ["white", "yellow"])
        ax.imshow(light_per_period, aspect='auto', cmap=light_cmap,
                  interpolation='none', extent=[0, self.tau, num_periods, 0],
                  alpha=0.5)
        ax.set_xticks(np.arange(0, self.tau + 1, 2))
        ax.set_xticklabels([f"{int(tick % 24):02d}" for tick in np.arange(0, self.tau + 1, 2)])
        ax.set_yticks(np.arange(0.5, num_periods, 1))
        ax.set_yticklabels([f"Period {i+1}" for i in range(num_periods)])
        ax.set_xlabel('Time of Day')
        ax.set_ylabel('Period')
        ax.set_title(f'Single-Plotted Actogram for {self.event_type}')
        plt.colorbar(cax, ax=ax, orientation='vertical', label='Event Count')

        legend_elements = []
        markers = []
        draggables = []

        if label_onset:
            for period, onset_bin in enumerate(self.onset_times):
                onset_hour = (onset_bin * self.bin_size_minutes) / 60
                marker, = ax.plot(onset_hour, period + 0.5, 'ro', picker=5)
                draggable = DraggableMarker(marker, self.onset_times, period, self.bin_size_minutes,
                                            on_drag=on_drag, day=period, label_type='onset',
                                            on_release_callback=on_release_callback)
                draggables.append(draggable)
                markers.append((marker, period, 'onset'))
            legend_elements.append(Line2D([0], [0], marker='o', color='w', label='Onset',
                                          markerfacecolor='r', markersize=8))

        if label_offset:
            for period, offset_bin in enumerate(self.offset_times):
                offset_hour = (offset_bin * self.bin_size_minutes) / 60
                marker, = ax.plot(offset_hour, period + 0.5, 'bo', picker=5)
                draggable = DraggableMarker(marker, self.offset_times, period, self.bin_size_minutes,
                                            on_drag=on_drag, day=period, label_type='offset',
                                            on_release_callback=on_release_callback)
                draggables.append(draggable)
                markers.append((marker, period, 'offset'))
            legend_elements.append(Line2D([0], [0], marker='o', color='w', label='Offset',
                                          markerfacecolor='b', markersize=8))

        if label_acrophase:
            for period, acrophase in enumerate(self.acrophase_times):
                acrophase_time_in_hours = acrophase[0] + acrophase[1] / 60.0
                marker, = ax.plot(acrophase_time_in_hours, period + 0.5, 'g^', picker=5)
                draggable = DraggableMarker(marker, self.acrophase_times, period, self.bin_size_minutes,
                                            is_tuple=True, on_drag=on_drag, day=period,
                                            label_type='acrophase', on_release_callback=on_release_callback)
                draggables.append(draggable)
                markers.append((marker, period, 'acrophase'))
            legend_elements.append(Line2D([0], [0], marker='^', color='w', label='Acrophase',
                                          markerfacecolor='g', markersize=8))

        if label_bathyphase:
            for period, bathyphase in enumerate(self.bathophase_times):
                bathyphase_time_in_hours = bathyphase[0] + bathyphase[1] / 60.0
                marker, = ax.plot(bathyphase_time_in_hours, period + 0.5, 'mv', picker=5)
                draggable = DraggableMarker(marker, self.bathophase_times, period, self.bin_size_minutes,
                                            is_tuple=True, on_drag=on_drag, day=period,
                                            label_type='bathyphase', on_release_callback=on_release_callback)
                draggables.append(draggable)
                markers.append((marker, period, 'bathyphase'))
            legend_elements.append(Line2D([0], [0], marker='v', color='w', label='Bathyphase',
                                          markerfacecolor='m', markersize=8))

        plt.tight_layout()
        if legend_elements:
            ax.legend(handles=legend_elements)
            self.draggable_markers.extend(draggables)        
        return fig, markers

    def double_actogram(self, label_onset=False, label_offset=False,
                        label_acrophase=False, label_bathyphase=False,
                        on_drag=None, on_release_callback=None):
        """
        Generate and return a double-plotted actogram with a light-data overlay.
        
        Parameters:
            label_onset (bool): Whether to mark onsets.
            label_offset (bool): Whether to mark offsets.
            label_acrophase (bool): Whether to mark acrophase.
            label_bathyphase (bool): Whether to mark bathyphase.
            on_drag (callable): Callback for marker dragging.
            on_release_callback (callable): Callback after marker drag.
        
        Returns:
            (matplotlib.figure.Figure, list): The figure and a list of marker info.
        
        Raises:
            ValueError: If tau is not set.
        """
        if self.tau is None:
            raise ValueError("Tau has not been determined yet.")
        self.bins_per_period = int((self.tau * 60) / self.bin_size_minutes)
        num_periods = len(self.padded_event_counts) // self.bins_per_period
        double_plotted_events = np.array([
            np.concatenate([
                self.padded_event_counts[i * self.bins_per_period:(i + 1) * self.bins_per_period],
                self.padded_event_counts[(i + 1) * self.bins_per_period:(i + 2) * self.bins_per_period]
                if (i + 1) < num_periods else np.zeros(self.bins_per_period)
            ])
            for i in range(num_periods)
        ])
        
        double_plotted_light = np.array([
            np.concatenate([
                self.padded_light[i * self.bins_per_period:(i + 1) * self.bins_per_period],
                self.padded_light[(i + 1) * self.bins_per_period:(i + 2) * self.bins_per_period]
                if (i + 1) < num_periods else np.zeros(self.bins_per_period)
            ])
            for i in range(num_periods)
        ])

        light_cmap = LinearSegmentedColormap.from_list("light_cmap", ["white", "yellow"])
        fig_height = max(6, double_plotted_events.shape[0] * 0.3)
        fig, ax = plt.subplots(figsize=(16, fig_height))
        ax.set_xticks(np.arange(0, 2 * self.tau + 1, 2))
        ax.set_xticklabels([f"{int(tick % 24):02d}" for tick in np.arange(0, 2 * self.tau + 1, 2)])
        cax = ax.imshow(double_plotted_events, aspect='auto', cmap='Greys',
                        interpolation='none', extent=[0, 2 * self.tau, double_plotted_events.shape[0], 0])
        ax.imshow(double_plotted_light, aspect='auto', cmap=light_cmap,
                  interpolation='none', extent=[0, 2 * self.tau, double_plotted_events.shape[0], 0],
                  alpha=0.5)
        ax.set_yticks(np.arange(0.5, double_plotted_events.shape[0], 1))
        ax.set_yticklabels([f"Period {i+1}" for i in range(double_plotted_events.shape[0])])
        ax.set_xlabel('Time of Day (Double Plotted)')
        ax.set_ylabel('Period')
        ax.set_title(f'Double-Plotted Actogram for {self.event_type}')
        ax.axvline(x=self.tau, color='red', linestyle='-', label=f'{self.tau:.2f} Hours')
        plt.colorbar(cax, ax=ax, orientation='vertical', label='Event Count')

        legend_elements = []
        markers = []
        draggables = []

        def create_double_markers(times_list, marker_style, color, label, is_tuple=False, label_type=None):
            for period in range(len(times_list)):
                time_value_period_i = times_list[period]
                if is_tuple:
                    time_in_hours_period_i = time_value_period_i[0] + time_value_period_i[1] / 60.0
                else:
                    time_in_hours_period_i = (time_value_period_i * self.bin_size_minutes) / 60.0

                marker1, = ax.plot(time_in_hours_period_i % self.tau, period + 0.5,
                                   marker_style, color=color, picker=5)
                draggable1 = DraggableMarker(marker1, times_list, period,
                                             self.bin_size_minutes, is_tuple=is_tuple,
                                             on_drag=on_drag, day=period, label_type=label_type,
                                             on_release_callback=on_release_callback)
                draggables.append(draggable1)
                markers.append((marker1, period, label_type))

                if period + 1 < len(times_list):
                    time_value_period_i_plus_1 = times_list[period + 1]
                    if is_tuple:
                        time_in_hours_period_i_plus_1 = time_value_period_i_plus_1[0] + time_value_period_i_plus_1[1] / 60.0
                    else:
                        time_in_hours_period_i_plus_1 = (time_value_period_i_plus_1 * self.bin_size_minutes) / 60.0

                    marker2, = ax.plot((time_in_hours_period_i_plus_1 % self.tau) + self.tau,
                                       period + 0.5, marker_style, color=color, picker=5)
                    draggable2 = DraggableMarker(marker2, times_list, period + 1,
                                                 self.bin_size_minutes, is_tuple=is_tuple,
                                                 on_drag=on_drag, day=period + 1, label_type=label_type,
                                                 on_release_callback=on_release_callback)
                    draggables.append(draggable2)
                    markers.append((marker2, period + 1, label_type))

                if period == 0:
                    legend_elements.append(Line2D([0], [0],
                                                  marker=marker_style, color='w', label=label,
                                                  markerfacecolor=color, markersize=8))

        # Markers for each label
        if label_onset:
            create_double_markers(self.onset_times, 'o', 'red', 'Onset', label_type='onset')
        if label_offset:
            create_double_markers(self.offset_times, 'o', 'blue', 'Offset', label_type='offset')
        if label_acrophase:
            create_double_markers(self.acrophase_times, '^', 'green', 'Acrophase', is_tuple=True, label_type='acrophase')
        if label_bathyphase:
            create_double_markers(self.bathophase_times, 'v', 'magenta', 'Bathyphase', is_tuple=True, label_type='bathyphase')

        plt.tight_layout()
        if legend_elements:
            ax.legend(handles=legend_elements)
            self.draggable_markers.extend(draggables)
        return fig, markers

    def plot_activity_profile(self, display_com=False, display_sem=False, base_on_tau=True):
        """
        Plot the average activity profile over one period.
        
        Parameters:
            display_com (bool): Whether to display the center-of-mass line.
            display_sem (bool): Whether to display SEM shading.
            base_on_tau (bool): If True, use tau; otherwise use 24 hours.
        
        Returns:
            matplotlib.figure.Figure: The activity profile figure.
        """
        if base_on_tau and self.tau is not None:
            bins_per_day = self.bins_per_period
            period_hours = self.tau
        else:
            bins_per_day = int(24 * 60 / self.bin_size_minutes)
            period_hours = 24
        total_periods = len(self.padded_event_counts) // bins_per_day
        truncated = self.padded_event_counts[:total_periods * bins_per_day]
        events_by_period = truncated.reshape((total_periods, bins_per_day))
        average_activity = np.mean(events_by_period, axis=0)
        sem_activity = np.std(events_by_period, axis=0) / np.sqrt(total_periods)
        time_axis = np.linspace(0, period_hours, bins_per_day, endpoint=False)
        com = np.sum(time_axis * average_activity) / np.sum(average_activity)
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(time_axis, average_activity, color='skyblue', linestyle='-', label='Average Activity')
        if display_sem:
            ax.fill_between(time_axis, average_activity - sem_activity, average_activity + sem_activity,
                            color='lightblue', alpha=0.5, label='SEM')
        if display_com:
            ax.axvline(x=com, color='red', linestyle='--', label=f'Center of Mass: {com:.2f} h')
        ax.set_xlabel('Time of Day (hours)')
        ax.set_ylabel('Average Activity Count')
        ax.set_title('Activity Profile')
        ax.set_xticks(np.arange(0, period_hours + 1, 2))
        ax.set_xlim(0, period_hours)
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        return fig

    def calculate_mean_on_off(self, onset_times, offset_times):
        """
        Calculate the average onset and offset times and return them as "HH:MM" strings.
        
        Parameters:
            onset_times (list[int]): List of onset bin indices.
            offset_times (list[int]): List of offset bin indices.
        
        Returns:
            (str, str): Average onset and offset times.
        """
        onset_hours = [(bin_val * self.bin_size_minutes) / 60 for bin_val in onset_times]
        onset_av = np.mean(onset_hours) if onset_hours else 0.0
        onset_time = f"{int(onset_av)%24:02d}:{int((onset_av-int(onset_av))*60):02d}"
        offset_hours = [(bin_val * self.bin_size_minutes) / 60 for bin_val in offset_times]
        offset_av = np.mean(offset_hours) if offset_hours else 0.0
        offset_time = f"{int(offset_av)%24:02d}:{int((offset_av-int(offset_av))*60):02d}"
        return onset_time, offset_time

    def results_writer(self, onset_res, offset_res, acro_res, batho_res, tau, onset_av, offset_av):
        """
        Combine analysis results into a DataFrame.
        
        Parameters:
            onset_res (list[str]): Onset times ("HH:MM").
            offset_res (list[str]): Offset times ("HH:MM").
            acro_res (list[str]): Acrophase times ("HH:MM").
            batho_res (list[str]): Bathyphase times ("HH:MM").
            tau (float): Estimated tau.
            onset_av (str): Average onset time.
            offset_av (str): Average offset time.
        
        Returns:
            pd.DataFrame: A DataFrame summarizing the results.
        """
        max_len = max(len(onset_res), len(offset_res), len(acro_res), len(batho_res))
        onset_res += [None] * (max_len - len(onset_res))
        offset_res += [None] * (max_len - len(offset_res))
        acro_res += [None] * (max_len - len(acro_res))
        batho_res += [None] * (max_len - len(batho_res))
        results = {
            'Onset_Times': onset_res,
            'Offset_Times': offset_res,
            'Acrophase_Times': acro_res,
            'Bathyphase_Times': batho_res,
            'Dominant_Period': [f'{tau:.2f}'] * max_len,
            'Onset_Average': [onset_av] * max_len,
            'Offset_Average': [offset_av] * max_len
        }
        df = pd.DataFrame(results)
        print(df)
        return df

    def csv_save(self, df, save_path):
        """
        Save the given DataFrame to a CSV file.
        
        Parameters:
            df (pd.DataFrame): The DataFrame to save.
            save_path (str): The target file path.
        """
        if not os.path.exists(save_path):
            df.to_csv(save_path, mode='w', header=True, index=False)
        else:
            df.to_csv(save_path, mode='a', header=False, index=False)

class DraggableMarker:
    """
    A helper class that makes matplotlib markers draggable.
    """
    def __init__(self, marker, times_list, index, bin_size_minutes, is_tuple=False,
                 on_drag=None, day=None, label_type=None, on_release_callback=None):
        """
        Initialize the DraggableMarker.

        Parameters:
            marker (matplotlib.lines.Line2D): The marker artist.
            times_list (list): List storing marker times.
            index (int): Index of the marker.
            bin_size_minutes (int): Bin size in minutes.
            is_tuple (bool): Whether times_list stores (hour, minute) tuples.
            on_drag (callable): Callback during dragging.
            day (int): The period index.
            label_type (str): Type of marker (e.g., 'onset').
            on_release_callback (callable): Callback after dragging.
        """
        self.marker = marker
        self.times_list = times_list
        self.index = index
        self.bin_size_minutes = bin_size_minutes
        self.is_tuple = is_tuple
        self.on_drag = on_drag
        self.day = day
        self.label_type = label_type
        self.on_release_callback = on_release_callback
        self.press = None
        self.background = None
        self.canvas = marker.figure.canvas
        self.connect()

    def connect(self):
        """Connect mouse events for dragging."""
        self.cidpress = self.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = self.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        """Handle mouse press event."""
        if event.inaxes != self.marker.axes:
            return
        contains, _ = self.marker.contains(event)
        if not contains:
            return
        self.press = (self.marker.get_xdata()[0], self.marker.get_ydata()[0]), event.xdata, event.ydata
        self.marker.set_animated(True)
        self.background = self.canvas.copy_from_bbox(self.marker.axes.bbox)
        self.marker.axes.draw_artist(self.marker)
        self.canvas.blit(self.marker.axes.bbox)

    def on_motion(self, event):
        """Handle mouse motion event."""
        if self.press is None or event.inaxes != self.marker.axes:
            return
        (x0, y0), xpress, ypress = self.press
        dx = event.xdata - xpress
        new_x = x0 + dx
        self.marker.set_data([new_x], [y0])
        if self.is_tuple:
            total_minutes = new_x * 60
            hours = int(total_minutes // 60)
            minutes = int(total_minutes % 60)
            self.times_list[self.index] = (hours % 24, minutes)
            new_time = hours + minutes / 60.0
        else:
            new_bin = int((new_x % 24) * 60 / self.bin_size_minutes)
            self.times_list[self.index] = new_bin
            new_time = new_x % 24
        self.canvas.restore_region(self.background)
        self.marker.axes.draw_artist(self.marker)
        self.canvas.blit(self.marker.axes.bbox)
        if self.on_drag is not None:
            self.on_drag(self.day, self.label_type, new_time)

    def on_release(self, event):
        """Handle mouse release event."""
        if self.press is None:
            return
        self.press = None
        self.marker.set_animated(False)
        self.background = None
        self.canvas.draw()
        if self.on_release_callback:
            self.on_release_callback()

    def disconnect(self):
        """Disconnect all events."""
        self.canvas.mpl_disconnect(self.cidpress)
        self.canvas.mpl_disconnect(self.cidrelease)
        self.canvas.mpl_disconnect(self.cidmotion)
