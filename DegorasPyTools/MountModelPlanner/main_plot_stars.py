import tkinter as tk
from astropy.time import Time
from class_star_file_manager import StarDataManager
from class_star_predictor import StarPredictor
from class_star_plot import StarPlot

# ----------------------------------------------------------------------------------------------------------------------
# WARNING
# THIS IS AN AUXILIARY SOFTWARE FOR TESTING PURPOSES ONLY.
# Angel Vera Herrera
# ----------------------------------------------------------------------------------------------------------------------

# Basic configuration
json_file = "star_raw_parameters_example.json"
json_file = "star_raw_parameters_all.json"
lat = 36.465257734376407939
lon = -6.20530535896
alt = 98.2496715541929
min_alt = 20

# Close function.

# MAIN
if __name__ == "__main__":

    # Prepare the location.
    data_manager = StarDataManager(json_file)
    stars_data = data_manager.read_json()

    # Prepare the predictor.
    star_predictor = StarPredictor(lat, lon, alt)

    # Prepare the GUI.
    window = tk.Tk()
    window.title("Star Polar Plot")
    window.geometry("800x800")

    # Create an instance of StarPlot
    plot = StarPlot(window, star_predictor, stars_data)

    # Close procedure.
    def on_closing():
        if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
            plot.close_plot()
            window.quit()
            window.destroy()
    window.protocol("WM_DELETE_WINDOW", on_closing)

    # Run the Tkinter event loop
    window.mainloop()