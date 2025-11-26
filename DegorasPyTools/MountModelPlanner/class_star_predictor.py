from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u


class StarPredictor:

    def __init__(self, lat, lon, alt):

        location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=alt * u.m)
        self.location = location

    def ra_dec_to_sky_coord_j2000(self, ra_h, ra_m, ra_s, dec_h, dec_m, dec_s):

        ra = f"{ra_h}h{ra_m}m{ra_s}s"
        dec = f"{dec_h}d{dec_m}m{dec_s}s"
        coordinates = SkyCoord(ra, dec, frame='icrs')
        return coordinates

    def calculate_az_el(self, star_coordinates, time):

        altaz_frame = AltAz(obstime=time, location=self.location)
        altaz_coordinates = star_coordinates.transform_to(altaz_frame)
        elevation = altaz_coordinates.alt.deg  # Elevation in degrees
        azimuth = altaz_coordinates.az.deg  # Azimuth in degrees
        return azimuth, elevation

    def predict(self, stars_data, times):

        stars_predictions = []

        # Iterate through each star in the stars_data
        for star in stars_data['stars']:
            star_predictions = {
                'star': star,
                'predictions': []
            }

            # Get RA/DEC coordinates in J2000 (ICRS) frame
            star_coordinates = self.ra_dec_to_sky_coord_j2000(
                star['ra_h'], star['ra_m'], star['ra_s'],
                star['dec_d'], star['dec_m'], star['dec_s']
            )

            # Calculate azimuth and elevation for each time
            for time in times:
                azimuth, elevation = self.calculate_az_el(star_coordinates, time)
                # Append the time, azimuth, and elevation to the star's prediction
                star_predictions['predictions'].append({
                    'time': time,
                    'azimuth': azimuth,
                    'elevation': elevation
                })

            # Store the predictions for this star
            stars_predictions.append(star_predictions)

        return stars_predictions