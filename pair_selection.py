import os
import json
import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from tqdm import tqdm

c_km_s = 299792.458
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Input and output
base_catalog = "/home/josephine/Sample/Catalogs/"
base_output = "/home/josephine/Sample/"

# Calculates angular distance using spherical law of cosines
def angular_distance(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = np.radians([ra1, dec1, ra2, dec2])
    cos_angle = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    return np.arccos(np.clip(cos_angle, -1, 1))
    
# Function to append pair to the respective FITS file
def append_to_fits_file(filename, best_pair, ra, dec, z):
    try:
        col1 = fits.Column(name='galaxy_1', format='K', array=[best_pair[0]])
        col2 = fits.Column(name='galaxy_2', format='K', array=[best_pair[1]])
        col3 = fits.Column(name='ra_1', format='E', array=[ra[best_pair[0]]])
        col4 = fits.Column(name='dec_1', format='E', array=[dec[best_pair[0]]])
        col5 = fits.Column(name='z_1', format='E', array=[z[best_pair[0]]])
        col6 = fits.Column(name='ra_2', format='E', array=[ra[best_pair[1]]])
        col7 = fits.Column(name='dec_2', format='E', array=[dec[best_pair[1]]])
        col8 = fits.Column(name='z_2', format='E', array=[z[best_pair[1]]])
        col9 = fits.Column(name='rp', format='E', array=[best_pair[2]])
        col10 = fits.Column(name='dv', format='E', array=[best_pair[3]])
        col11 = fits.Column(name='angular_distance', format='E', array=[best_pair[4]])
        col12 = fits.Column(name='angular_diameter_distance', format='E', array=[best_pair[5]])

        cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12])

        if os.path.exists(filename):
            with fits.open(filename, mode='update') as hdul:
                new_hdu = fits.BinTableHDU.from_columns(cols)
                hdul[1].data = np.hstack([hdul[1].data, new_hdu.data])
        else:
            fits.BinTableHDU.from_columns(cols).writeto(filename, overwrite=True)
    except Exception as e:
        print(f"Error writing to {filename}: {e}")

for subset_num in range(1, 10):
    print(f"\n===== Processing subset {subset_num} =====")

    input_file = os.path.join(base_catalog, f"subset{subset_num}.fits")
    output_dir = os.path.join(base_output, f"subset_{subset_num}")
    os.makedirs(output_dir, exist_ok=True)

    progress_file = f"{output_dir}/progress.json"

    try:
        with fits.open(input_file) as hdul:
            data = hdul[1].data
    except Exception as e:
        print(f"Error loading {input_file}: {e}")
        continue

    ra = data['ra_1']
    dec = data['dec_1']
    z = data['z']
    N = len(z)

    if os.path.exists(progress_file):
        with open(progress_file, 'r') as f:
            start_i = json.load(f).get("last_processed_index", 0)
    else:
        start_i = 0
 # Pair every galaxy i, with every other galaxy j

    for i in tqdm(range(start_i, N), desc=f'Subset {subset_num}'):
        best_pair = None
        best_file_name = None

        for j in range(N):
            if i == j:
                continue

            z1, z2 = z[i], z[j]
            ang_dist = angular_distance(ra[i], dec[i], ra[j], dec[j])
            ang_diam_dist = cosmo.angular_diameter_distance(z1).value * 1000 # converts to kpc
            rp = ang_dist * ang_diam_dist  # calculates projected distance separation in kpc
            dv = np.abs(c_km_s * (z1 - z2)) # velocity differences in kpc
        # Conditions to categorise the pairs
            file_name = None
            if rp < 20:
                file_name = 'A1' if dv < 500 else 'A2' if dv < 1000 else None
            elif rp < 50:
                file_name = 'B1' if dv < 500 else 'B2' if dv < 1000 else None
            elif rp < 100:
                file_name = 'C1' if dv < 500 else 'C2' if dv < 1000 else None
            elif rp < 250:
                file_name = 'D1' if dv < 500 else 'D2' if dv < 1000 else None
            elif rp < 500:
                file_name = 'E1' if dv < 500 else 'E2' if dv < 1000 else None
            elif rp < 1000:
                file_name = 'F1' if dv < 500 else 'F2' if dv < 1000 else None
            
          # Determine if this is the best pair for galaxy i
            if file_name and (best_pair is None or rp < best_pair[2] or (rp == best_pair[2] and dv < best_pair[3])):
                best_pair = (i, j, rp, dv, ang_dist, ang_diam_dist)
                best_file_name = file_name
                
    # Save the best pair for galaxy i if a valid pair was found
        if best_pair:
            append_to_fits_file(os.path.join(output_dir, f"{best_file_name}.fits"),
                               best_pair, ra, dec, z)
    # Update progress file
        with open(progress_file, 'w') as f:
            json.dump({"last_processed_index": i}, f)

    print(f"Subset {subset_num} complete!")

print(" All 9 subsets processed successfully!")

