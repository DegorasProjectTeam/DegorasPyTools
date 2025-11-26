import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import csv
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.time import Time
import datetime
import time
import threading


class StarPlot:

    def __init__(self, window, star_predictor, stars_data=[],
                 default_time=datetime.datetime.utcnow().isoformat(timespec='seconds')):
        self.min_el = 20
        self.stars_points = None
        self.stars_predictions = []
        self.default_time = default_time
        self.star_predictor = star_predictor
        self.fig, self.ax = plt.subplots(subplot_kw={'projection': 'polar'})
        self.canvas = None
        self.info_label = None  # Label to display star info in the GUI
        self.time_entry = None  # Time input field
        self.selected_element = None  # Keep track of the selected star's plot element
        self.plot_frame = None  # Frame to hold the plot
        self.observed_positions = []  # Track the observed positions (name, theta, r, datetime, catalog details)
        self.auto_update_time = tk.BooleanVar(value=False)  # For enabling/disabling auto time update
        self.auto_time_thread = None  # Thread for automatic time update
        self.stop_auto_update = threading.Event()  # Event to stop auto-update
        self.init_plot(window)
        self.min_magnitude = 0
        self.max_magnitude = 0
        self.update_stars(stars_data)

    def init_plot(self, window):

        # Create a frame to hold the plot and allow it to resize
        self.plot_frame = tk.Frame(window)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

        # Set up the polar plot with azimuth and elevation limits
        self.ax.set_theta_zero_location('N')
        self.ax.set_theta_direction(-1)
        self.ax.set_ylim(0, 90)
        self.ax.set_yticks([13, 24, 35, 46, 57, 68, 79])
        self.ax.set_yticklabels(['79', '68°', '57°', '46°', '35°', '24°', '13°'])
        self.ax.set_xticks(np.deg2rad([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]))
        self.ax.set_xticklabels(['0', '30°', '60°', '90°', '120°', '150°', '180°', '210°', '240°', '270°', '300°', '330°'])

        # Add a label to display star info
        self.info_label = tk.Label(window, text="Click on a star to see details", font=("Arial", 12))
        self.info_label.pack(pady=10)
        self.info_label.config(text=f"Name: - | Catalog: - | Number: - | Magnitude: - \nAz: - | El: -")

        # Create a frame to align the time entry and update button horizontally
        time_frame = tk.Frame(window)
        time_frame.pack(pady=10)

        # Add an entry for time input inside the time frame
        self.time_entry = tk.Entry(time_frame)
        self.time_entry.pack(side=tk.LEFT)
        self.time_entry.insert(0, self.default_time)  # Default time format

        # Add a button to update the time inside the time frame
        update_button = tk.Button(time_frame, text="Update Time", command=self.update_time)
        update_button.pack(side=tk.LEFT, padx=5)

        # Add a checkbox to enable real-time updates
        auto_update_checkbox = tk.Checkbutton(time_frame, text="Auto-Update Time (UTC)", variable=self.auto_update_time, command=self.toggle_auto_time)
        auto_update_checkbox.pack(side=tk.LEFT, padx=5)

        # Create another frame for aligning the buttons horizontally
        button_frame = tk.Frame(window)
        button_frame.pack(pady=10)

        # Add a button to mark the current star as observed inside button frame
        mark_button = tk.Button(button_frame, text="Mark position as observed", command=self.mark_as_observed)
        mark_button.pack(side=tk.LEFT, padx=5)

        # Add a button to save observed positions to a CSV file inside button frame
        save_button = tk.Button(button_frame, text="Save Observed Positions", command=self.save_observed_positions)
        save_button.pack(side=tk.LEFT, padx=5)

        # Add a button to load observed positions from a CSV file inside button frame
        load_button = tk.Button(button_frame, text="Load Observed Positions", command=self.load_observed_positions)
        load_button.pack(side=tk.LEFT, padx=5)

        # Set up the canvas and Tkinter event handling
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)  # Make canvas expand and fill the frame
        self.canvas.mpl_connect('pick_event', self.on_pick)

        # Allow the window and the plot to resize
        self.plot_frame.pack_propagate(False)  # Prevent the frame from shrinking

    def update_stars(self, stars_data):

        filtered_stars = []
        for star in stars_data['stars']:
            if float(star.get('Vmag', 99)) <= 5:
                filtered_stars.append(star)

        self.stars_data = {'stars': filtered_stars}

        if len(stars_data) > 0:
            # First, calculate the minimum and maximum magnitudes in the data
            magnitudes = [star.get('Vmag', 1.0) for star in self.stars_data['stars']]
            self.min_magnitude = min(magnitudes)
            self.max_magnitude = max(magnitudes)
            self.update_time()

    def update_time(self):

        new_time_str = self.time_entry.get()

        try:
            # Parse the new time from the input
            new_time = Time(new_time_str)

            # Generate new predictions for the stars based on the new time
            new_predictions = self.star_predictor.predict(self.stars_data, [new_time])

            # Update the plot with the new predictions
            self.update_predictions(new_predictions)

        except ValueError:
            # If the time input is invalid, show an error in the info label
            self.info_label.config(text="Invalid time format. Use YYYY-MM-DDTHH:MM:SS")

    def update_predictions(self, new_predictions):

        # Store the currently selected star's name (if a star is selected)
        selected_star_name = self.selected_element[1]['name'] if self.selected_element else None

        # Update the internal star predictions with the new data
        self.stars_predictions = new_predictions

        # Clear the existing plot
        self.ax.clear()
        self.ax.set_theta_zero_location('N')
        self.ax.set_theta_direction(-1)
        self.ax.set_ylim(0, 90)
        self.ax.set_yticks([13, 24, 35, 46, 57, 68, 79])
        self.ax.set_yticklabels(['79', '68°', '57°', '46°', '35°', '24°', '13°'])
        self.ax.set_xticks(np.deg2rad([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]))
        self.ax.set_xticklabels(['0', '30°', '60°', '90°', '120°', '150°', '180°', '210°', '240°', '270°', '300°', '330°'])

        # Re-plot the stars with the new predictions
        self.stars_points = []

        # Re-plot all stars with black borders and sizes depending on magnitude
        for star_data in self.stars_predictions:

            predictions = star_data['predictions']

            if len(predictions) == 1:
                # Single prediction: plot a point
                azimuth = predictions[0]['azimuth']
                elevation = predictions[0]['elevation']

                if elevation > self.min_el:  # Only plot if the star is above the horizon

                    star = star_data['star']
                    magnitude = star.get('Vmag', 1.0)
                    size = int(100 - (magnitude - self.min_magnitude) * (
                            95 / (self.max_magnitude - self.min_magnitude + 1e-5)))
                    r = 90 - elevation  # Convert elevation to polar coordinates
                    theta = np.deg2rad(azimuth)  # Convert azimuth to radians

                    # Plot using scatter for custom size and black edge color
                    point = self.ax.scatter(theta, r, s=size, edgecolor='black', facecolor='blue', zorder=5, picker=6)
                    self.stars_points.append((point, star, r, theta))  # Store point and star for interaction

        # Re-plot observed positions
        for (star_name, theta, r, _, _, _) in self.observed_positions:
            self.ax.plot(theta, r, 'g*', markersize=12, zorder=10)  # Plot observed position with green asterisk

        # If a star was selected, find the updated position and keep it selected
        if selected_star_name:
            for plot_element, star, r, theta in self.stars_points:
                if star['name'] == selected_star_name:
                    # Update the selected element to reflect the new position
                    self.selected_element = (plot_element, star, r, theta)
                    plot_element.set_facecolor('red')  # Keep the selected star red
                    azimuth = np.rad2deg(theta)
                    elevation = 90 - r
                    # Update the star information label
                    self.info_label.config(
                        text=f"Name: {star['name']} | Catalog: {star['catalog_name']} | Number: {star['catalog_num']} | Magnitude: {round(star['Vmag'],4)}"
                             f"\nAz: {azimuth:.4f}° | El: {elevation:.4f}°"
                    )

        # Redraw the canvas to reflect the new data
        self.canvas.draw_idle()

    def close_plot(self):

        # Stop the auto-update thread before closing
        self.stop_auto_update.set()

        # Ensure any thread completes before closing the window
        if self.auto_time_thread:
            self.auto_time_thread.join()

    def mark_as_observed(self):
        """Mark the currently selected star's position as observed."""
        if self.selected_element:
            plot_element, star, r, theta = self.selected_element
            if r is not None and theta is not None:
                # Add the observed position to the list, along with the current datetime and catalog info
                observed_datetime = datetime.datetime.now(datetime.UTC).isoformat(timespec='milliseconds')
                catalog_name = star.get('catalog_name', 'Unknown')
                catalog_num = star.get('catalog_num', 'Unknown')

                self.observed_positions.append((star['name'], theta, r, observed_datetime, catalog_name, catalog_num))

                # Redraw the plot with the observed marker
                self.update_predictions(self.stars_predictions)

    def save_observed_positions(self):
        """Save the observed positions to a CSV file."""
        if not self.observed_positions:
            self.info_label.config(text="No observed positions to save.")
            return

        # Create the filename with the current datetime
        filename = f"observed_positions_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"

        with open(filename, 'w', newline='') as file:
            writer = csv.writer(file, delimiter=';')
            writer.writerow(['datetime', 'star_name', 'catalog_name', 'catalog_num', 'azimuth', 'elevation'])
            for star_name, theta, r, observed_datetime, catalog_name, catalog_num in self.observed_positions:
                azimuth_deg = np.rad2deg(theta)
                elevation_deg = 90 - r
                writer.writerow([observed_datetime, star_name, catalog_name, catalog_num, f"{azimuth_deg:.2f}", f"{elevation_deg:.2f}"])

        self.info_label.config(text=f"Observed positions saved to {filename}")

    def load_observed_positions(self):
        """Load observed positions from a CSV file."""
        filename = filedialog.askopenfilename(title="Select file", filetypes=[("CSV files", "*.csv")])
        if not filename:
            return

        with open(filename, 'r') as file:
            reader = csv.reader(file, delimiter=';')
            next(reader)  # Skip header row
            self.observed_positions = []
            for row in reader:
                observed_datetime, star_name, catalog_name, catalog_num, az, el = row
                theta = np.deg2rad(float(az))
                r = 90 - float(el)
                self.observed_positions.append((star_name, theta, r, observed_datetime, catalog_name, catalog_num))

        self.update_predictions(self.stars_predictions)
        self.info_label.config(text=f"Observed positions loaded from {filename}")

    def toggle_auto_time(self):

        if self.auto_update_time.get():
            self.stop_auto_update.clear()
            self.auto_time_thread = threading.Thread(target=self.auto_update_time_loop)
            self.auto_time_thread.start()
        else:
            self.stop_auto_update.set()

    def auto_update_time_loop(self):

        last_plot_update = time.time()  # Track the last time the plot was updated

        while not self.stop_auto_update.is_set():

            current_time = time.time()

            # Update the current UTC time in the time entry
            current_utc_time = datetime.datetime.utcnow().isoformat(timespec='seconds')
            self.time_entry.delete(0, tk.END)
            self.time_entry.insert(0, current_utc_time)

            # Update the selected star's position and label every 0.5 seconds
            if self.selected_element:
                plot_element, star, r, theta = self.selected_element
                new_time = Time(current_utc_time)  # Use the current time for prediction

                # Predict the new azimuth and elevation for the selected star
                star_coordinates = self.star_predictor.ra_dec_to_sky_coord_j2000(
                    star['ra_h'], star['ra_m'], star['ra_s'], star['dec_d'], star['dec_m'], star['dec_s']
                )
                azimuth, elevation = self.star_predictor.calculate_az_el(star_coordinates, new_time)

                # Update the star's position (theta, r) based on new predictions
                theta = np.deg2rad(azimuth)
                r = 90 - elevation

                # Update the label to show the current azimuth, elevation, and other info
                magnitude = round(star['Vmag'], 4)
                catalog_name = star['catalog_name']
                catalog_num = star['catalog_num']

                # Update the star information label with the new predicted values
                self.info_label.config(
                    text=f"Name: {star['name']} | Catalog: {catalog_name} | Number: {catalog_num} | Magnitude: {magnitude}"
                         f"\nAz: {azimuth:.4f}° | El: {elevation:.4f}°")

            # Update the entire plot every 10 seconds
            if current_time - last_plot_update >= 120:
                self.update_time()  # Refresh the entire plot
                current_utc_time = datetime.datetime.utcnow().isoformat(timespec='seconds')
                last_plot_update = current_time  # Reset the last plot update time

            time.sleep(10)  # Sleep for 0.5 seconds before the next loop iteration

    # Event handler for clicking on stars
    def on_pick(self, event):

        # Restore the color of the previously selected element
        if self.selected_element:
            self.selected_element[0].set_facecolor('blue')  # Reset the color to blue for previous selection

        # Find the star that was clicked on
        for plot_element, star, r, theta in self.stars_points:

            if event.artist == plot_element:
                star_name = star['name']
                catalog_num = star['catalog_num']
                catalog_name = star['catalog_name']
                magnitude = round(star['Vmag'], 4)

                # Check if only one prediction exists for this star
                for star_data in self.stars_predictions:
                    if star_data['star']['name'] == star_name:
                        predictions = star_data['predictions']

                        if len(predictions) == 1:
                            # Extract azimuth and elevation from the single prediction
                            azimuth = predictions[0]['azimuth']
                            elevation = predictions[0]['elevation']
                            self.info_label.config(
                                text=f"Name: {star_name} | Catalog: {catalog_name} | Number: {catalog_num} | Magnitude: {magnitude}"
                                     f"\nAz: {azimuth:.4f}° | El: {elevation:.4f}°"
                            )
                        else:
                            # If more than one prediction, show just the star information
                            self.info_label.config(
                                text=f"Name: {star_name} | Name: {catalog_name} | Number: {catalog_num}\nAz: - | El: -"
                            )

                # Change the color of the clicked element to red
                plot_element.set_facecolor('red')  # Set the color to red for new selection

                # Keep track of the selected element
                self.selected_element = (plot_element, star, r, theta)

                # Redraw the canvas to show the updated plot
                self.canvas.draw()

                break
