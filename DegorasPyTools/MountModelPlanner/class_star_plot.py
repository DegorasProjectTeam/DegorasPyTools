import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import csv
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.time import Time
import datetime
import time
from matplotlib import colors as mcolors
import matplotlib.patheffects as patheffects
import matplotlib as mpl

plt.rcParams['font.family'] = 'Segoe UI'


class StarPlot:

    def __init__(self, window, star_predictor, stars_data=[],
                 default_time=datetime.datetime.utcnow().isoformat(timespec='seconds')):

        # Color configuration
        self.color_bg_hex = "#000885"
        self.color_star_edge_hex = "#F39C12"  # Orange border
        self.color_star_fill_hex = "#F4D03F"  # Warm yellow fill
        self.color_star_selected_hex = "#E74C3C"  # Reddish for selected star
        self.color_star_fill_rgba = np.array([0.957, 0.816, 0.247, 1.0])  # #F4D03F
        self.color_star_selected_rgba = np.array([0.784, 0.133, 0.117, 1.0])  # ~#C9221E

        # Init members.
        self.min_el = 15
        self.max_el = 85
        self.default_time = default_time
        self.star_predictor = star_predictor
        self.fig, self.ax = plt.subplots(subplot_kw={'projection': 'polar'})
        self.canvas = None
        self.info_label = None
        self.time_entry = None
        self.selected_element = None
        self.plot_frame = None
        self.observed_positions = []
        self.auto_update_time = tk.BooleanVar(value=False)
        self.init_plot(window)
        self.min_magnitude = 0
        self.max_magnitude = 0
        self.stars_predictions = []
        self.scatter = None
        self.thetas = None
        self.rs = None
        self.visible_stars = []
        self.visible_predictions = []
        self.selected_index = None
        self.selected_star = None
        self.auto_update_job = None



        # Update the stars.
        self.update_stars(stars_data)

    def _configure_axes(self):
        """Configure polar axes (background, grid, azimuth & elevation ticks/labels)."""

        # --- Base polar config & background ---
        self.ax.set_facecolor(self.color_bg_hex)
        self.ax.set_theta_zero_location('N')
        self.ax.set_theta_direction(-1)
        self.ax.set_ylim(0, 80)

        # --- Elevation (El) major/minor ticks ---
        # Major elevation rings
        el_major_deg = [15, 25, 35, 45, 55, 65, 75, 85]
        r_major = np.array([90 - el for el in el_major_deg], dtype=float)

        # Minor elevation rings: one division between each pair of major rings
        el_minor_deg = 0.5 * (np.array(el_major_deg[:-1]) + np.array(el_major_deg[1:]))
        r_minor = np.array([90 - el for el in el_minor_deg], dtype=float)

        # Set major/minor radial ticks (no labels here, los pondremos a mano)
        self.ax.set_yticks(r_major)
        self.ax.set_yticks(r_minor, minor=True)

        # --- Azimuth (Az) major/minor ticks ---
        azi_major_deg = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
        azi_minor_deg = np.arange(0, 360, 10)  # 0,10,20,...,350

        self.ax.set_xticks(np.deg2rad(azi_major_deg))
        self.ax.set_xticklabels([f"{azi}" for azi in azi_major_deg])

        self.ax.set_xticks(np.deg2rad(azi_minor_deg), minor=True)

        # --- Shaded zone between 10° and 15° elevation (r from 80 to 75) ---
        theta = np.linspace(0, 2 * np.pi, 360)
        r_min = 90 - 15  # 75
        r_max = 90 - 10  # 80
        self.ax.fill_between(theta, r_min, r_max, color='gray', alpha=0.3, zorder=0)

        # --- Grid: major + minor (same style for Az and El) ---
        # Major grid: stronger
        self.ax.grid(
            which='major',
            color="#EAF4FF",
            alpha=0.5,
            linewidth=1.1
        )

        # Minor grid: lighter & dashed
        self.ax.grid(
            which='minor',
            color="#EAF4FF",
            alpha=0.3,
            linewidth=0.5,
            linestyle='--'
        )

        # --- Custom elevation labels (inside plot, near Az = 90°) ---
        # Remove automatic radial labels (major & minor)
        self.ax.set_yticklabels([])
        self.ax.set_yticklabels([], minor=True)

        theta_pos = np.deg2rad(90 + 1)  # Slight offset from pure 90° for readability

        for el, r in zip(el_major_deg, r_major):
            # Small radial offset inward so they do not sit exactly on the grid circle
            r_text = r - 1

            self.ax.text(
                theta_pos,
                r_text,
                f"{el}",
                ha='right',
                va='top',
                color="#EAF4FF",
                alpha=0.7,
                fontsize=7.5,
                bbox=dict(
                    boxstyle='round,pad=0.0',
                    facecolor=self.color_bg_hex,  # same as background
                    edgecolor='none',
                    alpha=1.0
                ),
                zorder=5
            )

        # --- Azimuth labels styling ---
        for label in self.ax.get_xticklabels():
            label.set_fontsize(8.5)
            label.set_color("#EAF4FF")
            label.set_fontweight('normal')

    def _apply_scatter_colors(self, selected_index: int | None):
        """
        Apply face colors & edges to scatter points:
          - Normal stars: warm yellow + orange border.
          - Selected star: reddish + bright border (thicker).
        """
        if self.scatter is None:
            return

        n = len(self.visible_stars)
        if n == 0:
            return

        facecolors = np.tile(self.color_star_fill_rgba, (n, 1))
        edgecolors = np.tile(np.array(mcolors.to_rgba(self.color_star_edge_hex)), (n, 1))
        linewidths = np.full(n, 1.0)

        # If something is selected, mark it more strongly
        if selected_index is not None and 0 <= selected_index < n:
            edgecolors[selected_index] = np.array(self.color_star_selected_rgba)
            linewidths[selected_index] = 3

        # Apply colors
        self.scatter.set_facecolors(facecolors)
        self.scatter.set_edgecolors(edgecolors)
        self.scatter.set_linewidths(linewidths)

    def init_plot(self, window):

        self.plot_frame = tk.Frame(window)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

        self._configure_axes()

        # Minimiza márgenes dentro de la figura
        self.fig.tight_layout(pad=0.2)
        self.fig.subplots_adjust(top=0.97, bottom=0.05, left=0.03, right=0.97)

        # Label info
        self.info_label = tk.Label(window, text="Click on a star to see details", font=("Arial", 12))
        self.info_label.pack(pady=3)

        time_frame = tk.Frame(window)
        time_frame.pack(pady=3)

        self.time_entry = tk.Entry(time_frame)
        self.time_entry.pack(side=tk.LEFT)
        self.time_entry.insert(0, self.default_time)

        update_button = tk.Button(time_frame, text="Update Time", command=self.update_time)
        update_button.pack(side=tk.LEFT, padx=5)

        auto_update_checkbox = tk.Checkbutton(time_frame, text="Auto-Update Time (UTC)",
                                              variable=self.auto_update_time, command=self.toggle_auto_time)
        auto_update_checkbox.pack(side=tk.LEFT, padx=5)

        button_frame = tk.Frame(window)
        button_frame.pack(pady=2)

        mark_button = tk.Button(button_frame, text="Mark position as observed", command=self.mark_as_observed)
        mark_button.pack(side=tk.LEFT, padx=5)

        save_button = tk.Button(button_frame, text="Save Observed Positions", command=self.save_observed_positions)
        save_button.pack(side=tk.LEFT, padx=5)

        load_button = tk.Button(button_frame, text="Load Observed Positions", command=self.load_observed_positions)
        load_button.pack(side=tk.LEFT, padx=5)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect('pick_event', self.on_pick)

        self.plot_frame.pack_propagate(False)

    def update_stars(self, stars_data):

        filtered_stars = []
        for star in stars_data['stars']:
            if float(star.get('Vmag', 99)) <= 8:
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
            t0 = time.perf_counter()
            new_predictions = self.star_predictor.predict(self.stars_data, [new_time])
            t1 = time.perf_counter()
            print(f"predict() took {1000 * (t1 - t0):.3f} ms")

            # Update the plot with the new predictions
            self.update_predictions(new_predictions)

        except ValueError:
            # If the time input is invalid, show an error in the info label
            self.info_label.config(text="Invalid time format. Use YYYY-MM-DDTHH:MM:SS")

    def _mag_to_size(self, mag: float) -> float:
        """
        Map visual magnitude to marker size (scatter 's' parameter).

        Brighter stars (lower/negative mag) get larger sizes.
        Fainter stars get smaller but still visible markers.
        """

        if mag <= 0.0:
            return 50.0   # very bright stars (Sirius-like)
        elif mag <= 1.0:
            return 40.0
        elif mag <= 2.0:
            return 30.0
        elif mag <= 3.0:
            return 20.0
        elif mag <= 4.0:
            return 15.0
        elif mag <= 5.0:
            return 10.0
        elif mag <= 6.0:
            return 5.0
        elif mag <= 7.0:
            return 5.0
        else:
            return 9.0    # faintest stars you are plotting

    def update_predictions(self, new_predictions):

        # Preserve previously selected star by name
        selected_star_name = self.selected_star['name'] if self.selected_star else None
        self.stars_predictions = new_predictions

        # Reset state structures for this update
        self.visible_stars = []
        self.visible_predictions = []
        self.thetas = []
        self.rs = []
        self.selected_index = None

        # Build arrays of visible stars
        for star_data in self.stars_predictions:
            predictions = star_data['predictions']
            if len(predictions) != 1:
                continue

            azimuth = predictions[0]['azimuth']
            elevation = predictions[0]['elevation']

            if elevation <= self.min_el or elevation >= self.max_el:
                continue

            star = star_data['star']
            self.visible_stars.append(star)
            self.visible_predictions.append(predictions[0])
            self.thetas.append(np.deg2rad(azimuth))  # az -> radians
            self.rs.append(90 - elevation)  # el -> polar radius

        # If there are no visible stars, remove scatter if it exists and exit
        if not self.thetas:
            if self.scatter is not None:
                self.scatter.remove()
                self.scatter = None
            self.canvas.draw_idle()
            return

        # Create or update the single scatter
        thetas_arr = np.array(self.thetas)
        rs_arr = np.array(self.rs)

        '''
        sizes_arr = np.array([
            int(4 * np.sqrt(max(0.1, self.max_magnitude - s.get('Vmag', 1.0))))
            for s in self.visible_stars
        ], dtype=float)
        '''

        sizes_arr = np.array([
            self._mag_to_size(s.get('Vmag', 1.0))
            for s in self.visible_stars
        ], dtype=float)

        if self.scatter is None:
            # Create unified scatter (first time)
            self.scatter = self.ax.scatter(
                thetas_arr,
                rs_arr,
                s=sizes_arr,
                linewidths=1,
                zorder=5,
                picker=6
            )
        else:
            # Update existing scatter
            offsets = np.column_stack([thetas_arr, rs_arr])
            self.scatter.set_offsets(offsets)
            self.scatter.set_sizes(sizes_arr)

        # Restore selected star if present
        restored_index = None
        if selected_star_name:
            for idx, star in enumerate(self.visible_stars):
                if star['name'] == selected_star_name:
                    self.selected_index = idx
                    self.selected_star = star
                    restored_index = idx

                    pred = self.visible_predictions[idx]
                    azimuth = pred['azimuth']
                    elevation = pred['elevation']
                    self.info_label.config(
                        text=f"Name: {star['name']} | Catalog: {star['catalog_name']} | Number: {star['catalog_num']} | Magnitude: {round(star['Vmag'], 4)}"
                             f"\nAz: {azimuth:.4f}° | El: {elevation:.4f}°"
                    )
                    break

        self._apply_scatter_colors(restored_index)

        # Re-plot observed positions (will add markers on top of scatter)
        for (star_name, theta, r, _, _, _) in self.observed_positions:
            self.ax.plot(theta, r, 'g*', markersize=12, zorder=10)

        self.canvas.draw_idle()

    def close_plot(self):
        """Cleanly stop auto-update and close the plot."""
        if self.auto_update_job is not None:
            self.plot_frame.after_cancel(self.auto_update_job)
            self.auto_update_job = None

    def mark_as_observed(self):
        """Mark the currently selected star's position as observed."""
        if self.selected_index is None or self.selected_star is None:
            self.info_label.config(text="No star selected to mark as observed.")
            return

        idx = self.selected_index
        star = self.selected_star

        if self.thetas is None or self.rs is None:
            self.info_label.config(text="Internal error: no positions available for selected star.")
            return

        if idx < 0 or idx >= len(self.thetas) or idx >= len(self.rs):
            self.info_label.config(text="Internal error: invalid selected index.")
            return

        theta = self.thetas[idx]
        r = self.rs[idx]

        if theta is None or r is None:
            self.info_label.config(text="Selected star has no valid position to mark as observed.")
            return

        observed_datetime = datetime.datetime.now(datetime.UTC).isoformat(timespec='milliseconds')
        catalog_name = star.get('catalog_name', 'Unknown')
        catalog_num = star.get('catalog_num', 'Unknown')

        self.observed_positions.append(
            (star['name'], theta, r, observed_datetime, catalog_name, catalog_num)
        )

        # Draw only the new observed marker
        self.ax.plot(
            theta, r,
            marker='X',  # Better symbol
            markersize=10,  # Larger size
            markerfacecolor='#32CD32',
            markeredgecolor='red',
            markeredgewidth=1.5,
            zorder=12
        )

        self.canvas.draw_idle()

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
        """Enable or disable automatic time updates using Tk's event loop."""
        if self.auto_update_time.get():
            if self.auto_update_job is None:
                self.auto_update_time_step()
        else:
            if self.auto_update_job is not None:
                self.plot_frame.after_cancel(self.auto_update_job)
                self.auto_update_job = None

    def auto_update_time_step(self):
        """
        Single auto-update step.

        This function:
          - Updates the time entry with current UTC time.
          - Updates the info label for the selected star (if any).
          - Optionally refreshes the full plot at a slower cadence.
          - Reschedules itself using Tk's 'after' if auto-update is still enabled.
        """

        # If auto-update was turned off meanwhile, do nothing
        if not self.auto_update_time.get():
            self.auto_update_job = None
            return

        # 1) Update current UTC time in the entry
        current_utc_time = datetime.datetime.utcnow().isoformat(timespec='seconds')
        self.time_entry.delete(0, tk.END)
        self.time_entry.insert(0, current_utc_time)

        # 2) Update selected star label (only label, not full plot)
        if self.selected_star is not None:
            star = self.selected_star
            new_time = Time(current_utc_time)

            # Compute new az/el for the currently selected star
            star_coordinates = self.star_predictor.ra_dec_to_sky_coord_j2000(
                star['ra_h'], star['ra_m'], star['ra_s'],
                star['dec_d'], star['dec_m'], star['dec_s']
            )
            azimuth, elevation = self.star_predictor.calculate_az_el(star_coordinates, new_time)

            magnitude = round(star['Vmag'], 4)
            catalog_name = star['catalog_name']
            catalog_num = star['catalog_num']

            self.info_label.config(
                text=f"Name: {star['name']} | Catalog: {catalog_name} | Number: {catalog_num} | Magnitude: {magnitude}"
                     f"\nAz: {azimuth:.4f}° | El: {elevation:.4f}°"
            )

        # 3) Optionally refresh the entire plot every N steps (optional)
        #    For example, refresh every 2 minutes.
        #    You can count calls using an attribute counter if needed.
        #    For now, you can keep manual refresh with the "Update Time" button
        #    or call self.update_time() less frequently from here if you want.

        # 4) Reschedule next auto-update (e.g. every 10 seconds)
        #    Adjust interval_ms as needed (1000 = 1 second, 10000 = 10 seconds)
        interval_ms = 1000
        self.auto_update_job = self.plot_frame.after(interval_ms, self.auto_update_time_step)

    # Event handler for clicking on stars
    def on_pick(self, event):

        # If there is no scatter, there is nothing to pick
        if self.scatter is None:
            return

        # Ensure the event comes from our scatter
        if event.artist is not self.scatter:
            return

        # Indices of the picked points (can be multiple)
        inds = event.ind
        if not len(inds):
            return

        # Use the first picked index
        idx = int(inds[0])

        # Safety check for index range
        if idx < 0 or idx >= len(self.visible_stars):
            return

        # Get star and prediction for this index
        star = self.visible_stars[idx]
        pred = self.visible_predictions[idx]

        # Update info label
        star_name = star['name']
        cat_name = star['catalog_name']
        cat_num = star['catalog_num']
        mag = round(star['Vmag'], 4)
        azimuth = pred['azimuth']
        elevation = pred['elevation']

        self.info_label.config(
            text=f"Name: {star_name} | Catalog: {cat_name} | Number: {cat_num} | Magnitude: {mag}"
                 f"\nAz: {azimuth:.4f}° | El: {elevation:.4f}°"
        )

        # Store current selection
        self.selected_index = idx
        self.selected_star = star

        # Apply colors: base + selected
        self._apply_scatter_colors(idx)

        # Non-blocking redraw
        self.canvas.draw_idle()
