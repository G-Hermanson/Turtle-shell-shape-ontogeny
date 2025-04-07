## This folder includes:

- `Script 1 (GMM, cluster).R` : code to read landmark data of Stayton et al. (2018) as well as new specimens into R environment. This script also runs allometric regressions of species-specific subsets of the data and cluster analysis to aggregate specimens into 'small', 'intermediate' and 'large' groups
- `Script 2 (disparity).R` : code to perform disparity analysis of turtle shell shape based on the pre-defined ontogenetic stages (see `Script 1` above)
- `Script 3 (final figures).R` : code to generate `png` and/or `pdf` files correspoding to figures in the Main and Supplementary texts
- `Stayton_landmarks.txt` : text file (`txt`) that is pasted into the Avizo console to create empty 3D landmark objects
- `new_info.csv` : simple `csv` containing the same metadata for the newly added specimens as the one used in Stayton et al. (2018), for matching purposes in `R`
- `sea_turtle_growth_Kordikova.csv` : data on carapace growth rate (cm/year) retrieved from Kordikova (2002) and from Bjorndal & Bolten (1988)
- `turtle_landm-1.csv` : raw 3D landmark coordinates of turtle shells retrieved from Stayton et al. (2018)
