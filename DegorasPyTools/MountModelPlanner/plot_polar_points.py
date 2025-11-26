import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import numpy as np


# Function to add a polar point
def add_point():
    try:
        azimuth = float(azimuth_entry.get()) * np.pi / 180  # Convert degrees to radians
        elevation = float(elevation_entry.get())

        # Convert elevation to polar coordinates
        r = 90 - elevation  # r represents the distance from the center in the plot (90°-elevation)
        ax.plot(azimuth, r, 'ro')  # Plot the point on the polar plot
        canvas.draw()
    except ValueError:
        print("Please enter valid numeric values for azimuth and elevation.")


# Function to create the polar plot
def create_polar_plot():
    global ax, canvas
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, polar=True)

    # Set the range for azimuth (0 to 360 degrees)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)

    # Set the range for elevation (0 to 90 degrees)
    ax.set_ylim(0, 90)
    ax.set_yticks([0, 30, 60, 90])  # Customize the elevation rings
    ax.set_yticklabels(['90°', '60°', '30°', '0°'])

    # Initialize the polar plot in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.get_tk_widget().pack()


# Main interface setup
def main():
    global window, azimuth_entry, elevation_entry
    window = tk.Tk()
    window.title("Azimuth and Elevation Polar Plot")

    # Create the plot area
    create_polar_plot()

    # Create input fields for azimuth and elevation
    tk.Label(window, text="Azimuth (0 to 360°):").pack()
    azimuth_entry = tk.Entry(window)
    azimuth_entry.pack()

    tk.Label(window, text="Elevation (0 to 90°):").pack()
    elevation_entry = tk.Entry(window)
    elevation_entry.pack()

    # Button to add the point to the plot
    add_point_button = tk.Button(window, text="Add Point", command=add_point)
    add_point_button.pack()

    # Start the Tkinter loop
    window.mainloop()


# Ensure the script runs only if executed directly
if __name__ == "__main__":
    main()