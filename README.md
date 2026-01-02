# Pair_selection
This repository contains Python code used to select and classify galaxy pairs from volume-limited samples. The script selects the **best companion** for each galaxy based on projected separation and velocity difference, and saves outputs into categorised FITS tables. Note that each galaxy is only paired once.

## Working of the code

For each galaxy in the input dataset:

1. Computes angular separation using the spherical law of cosines
2. Computes the angular-diameter distance from a Flat ΛCDM cosmology
3. Calculates **projected distance** (rp in kpc) using the angular distance and angular diameter distance
4.Calculates **velocity difference** (dv in km/s) using redshifts z1 and z2
5. Categorises galaxy pairs into six projected separation bins and two velocity differences.
   
Only the **closest valid companion** (minimum rp, then dv) is stored for each galaxy.
# Contact
For inquiries, contact:
Josephine Chishala
josephinechishala2@gmail.com
