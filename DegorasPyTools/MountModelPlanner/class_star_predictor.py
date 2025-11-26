from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
from astropy.time import Time
import astropy.units as u


class StarPredictor:

    def __init__(self, lat, lon, alt):

        self.location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=alt * u.m)

    def ra_dec_to_sky_coord_j2000(self, ra_h, ra_m, ra_s, dec_h, dec_m, dec_s):

        ra = f"{ra_h}h{ra_m}m{ra_s}s"
        dec = f"{dec_h}d{dec_m}m{dec_s}s"
        coordinates = SkyCoord(ra, dec, frame='icrs')
        return coordinates

    def calculate_az_el(self, star_coordinates, time):

        altaz_frame = AltAz(obstime=time, location=self.location)
        altaz_coordinates = star_coordinates.transform_to(altaz_frame)
        elevation = altaz_coordinates.alt.deg
        azimuth = altaz_coordinates.az.deg
        return azimuth, elevation

    def predict(self, stars_data, times):

        stars = stars_data.get('stars', [])
        if not stars:
            return []

        # --- Normalize 'times' to a Time array ---
        if isinstance(times, Time):
            times_arr = times
        else:
            times_arr = Time(times)

        # Ensure we have at least one instant
        if np.size(times_arr) == 0:
            return []

        # --- Build RA/Dec arrays in degrees for all stars ---
        n_stars = len(stars)

        ra_h = np.array([float(s['ra_h']) for s in stars], dtype=float)
        ra_m = np.array([float(s['ra_m']) for s in stars], dtype=float)
        ra_s = np.array([float(s['ra_s']) for s in stars], dtype=float)
        dec_d = np.array([float(s['dec_d']) for s in stars], dtype=float)
        dec_m = np.array([float(s['dec_m']) for s in stars], dtype=float)
        dec_s = np.array([float(s['dec_s']) for s in stars], dtype=float)

        # RA: h/m/s → deg (1h = 15 deg)
        ra_deg = (ra_h + ra_m / 60.0 + ra_s / 3600.0) * 15.0

        # Dec: d/m/s → deg, preserving sign in dec_d
        dec_abs = np.abs(dec_d) + dec_m / 60.0 + dec_s / 3600.0
        dec_deg = np.sign(dec_d) * dec_abs

        # --- Vector SkyCoord for all stars ---
        coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')

        # --- AltAz frame for all times ---
        altaz_frame = AltAz(obstime=times_arr, location=self.location)

        # Broadcasting: coords (N,) with times (T,) → altaz with shape (N, T) or (T, N)
        altaz = coords.transform_to(altaz_frame)

        elev = np.asarray(altaz.alt.deg)
        az = np.asarray(altaz.az.deg)

        # Normalize shapes to (n_stars, n_times)
        if elev.ndim == 0:
            # Single star, single time (very degenerate case)
            elev = elev.reshape(1, 1)
            az = az.reshape(1, 1)
            n_times = 1
        elif elev.ndim == 1:
            # Single time → (n_stars,) → make it (n_stars, 1)
            elev = elev.reshape(n_stars, 1)
            az = az.reshape(n_stars, 1)
            n_times = 1
        else:
            # 2D → could be (n_stars, n_times) or (n_times, n_stars)
            n_times_candidate0 = elev.shape[1]
            n_times_candidate1 = elev.shape[0]

            # We know how many stars and how many times we expect:
            # expected n_times from times_arr:
            n_times_expected = np.size(times_arr)

            if elev.shape[0] == n_stars and elev.shape[1] == n_times_expected:
                # shape (n_stars, n_times)
                n_times = elev.shape[1]
            elif elev.shape[0] == n_times_expected and elev.shape[1] == n_stars:
                # shape (n_times, n_stars) → transpose
                elev = elev.T
                az = az.T
                n_times = elev.shape[1]
            else:
                # Unexpected shape, but we try to fail loudly
                raise RuntimeError(
                    f"Unexpected alt/az shape {elev.shape}, "
                    f"expected (n_stars={n_stars}, n_times={n_times_expected}) or transposed."
                )

        # --- Pack results back into your original structure ---
        stars_predictions = []
        for i, star in enumerate(stars):
            preds = []

            if n_times == 1:
                # Single time case
                t = times_arr[0] if isinstance(times_arr, Time) and times_arr.isscalar is False else times_arr
                preds.append({
                    'time': t,
                    'azimuth': float(az[i, 0]),
                    'elevation': float(elev[i, 0])
                })
            else:
                # Multiple times
                for j in range(n_times):
                    t_j = times_arr[j]
                    preds.append({
                        'time': t_j,
                        'azimuth': float(az[i, j]),
                        'elevation': float(elev[i, j])
                    })

            stars_predictions.append({
                'star': star,
                'predictions': preds
            })

        return stars_predictions