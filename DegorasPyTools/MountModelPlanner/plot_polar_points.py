import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import numpy as np


# Configuración de bandas y sectores (ejemplo)
# el_min / el_max en grados; az_edges en grados
bands = [
    # Elevaciones bajas: 3 sectores en azimut
    {"el_min": 0.0, "el_max": 24.0, "az_edges": [0.0, 10.0, 20.0, 30.0]},
    # Elevaciones altas: 1 sector en azimut
    {"el_min": 24.0, "el_max": 90.0, "az_edges": [0.0, 30.0]},
]


def draw_sectors(ax):
    """
    Dibuja las divisiones de sectores en el diagrama polar, usando la
    convención r = 90 - elevación.
    """
    for band in bands:
        el_min = band["el_min"]
        el_max = band["el_max"]
        az_edges = band["az_edges"]

        # Conversión elevación -> radio
        # r = 90 - elevación
        r_outer = 90.0 - el_min  # límite exterior (más lejos del centro)
        r_inner = 90.0 - el_max  # límite interior (más cerca del centro)

        az_min = az_edges[0]
        az_max = az_edges[-1]

        # Grid fino de azimut para dibujar las líneas de elevación (arcos)
        theta = np.linspace(np.deg2rad(az_min), np.deg2rad(az_max), 200)

        # Líneas de elevación (bordes de la banda)
        ax.plot(theta, np.full_like(theta, r_inner), linestyle="-", linewidth=1)
        ax.plot(theta, np.full_like(theta, r_outer), linestyle="-", linewidth=1)

        # Líneas de azimut (bordes de sectores) dentro de esta banda
        for az in az_edges:
            t = np.deg2rad(az)
            r_line = np.linspace(r_inner, r_outer, 2)
            ax.plot([t, t], [r_inner, r_outer], linestyle="-", linewidth=1)

        # Etiquetas de sectores (opcional)
        for i in range(len(az_edges) - 1):
            az_c = 0.5 * (az_edges[i] + az_edges[i + 1])
            el_c = 0.5 * (el_min + el_max)
            t_c = np.deg2rad(az_c)
            r_c = 90.0 - el_c
            ax.text(
                t_c,
                r_c,
                f"S{i}",
                ha="center",
                va="center",
                fontsize=8,
            )


# Function to add a polar point
def add_point():
    try:
        azimuth = float(azimuth_entry.get()) * np.pi / 180.0  # Convert degrees to radians
        elevation = float(elevation_entry.get())

        # Convert elevation to polar coordinates
        r = 90.0 - elevation  # r representa la distancia desde el centro en el plot
        ax.plot(azimuth, r, "ro")  # Plot the point on the polar plot
        canvas.draw()
    except ValueError:
        print("Please enter valid numeric values for azimuth and elevation.")


# Function to create the polar plot
def create_polar_plot():
    global ax, canvas
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, polar=True)

    # Set the range for azimuth (0 to 360 degrees)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    # Set the range for elevation (0 to 90 degrees)
    ax.set_ylim(0, 90)
    ax.set_yticks([0, 30, 60, 90])  # Customize the elevation rings
    ax.set_yticklabels(["90°", "60°", "30°", "0°"])

    # Dibujar los sectores definidos
    draw_sectors(ax)

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
