# TopoLAB
Topography-based routines for generating terrain parameters (TopoPAR) as well as downscaling (TopoSCALE) and clustering (TopoSUB) meteorological forcing data. These routines have been ported and adapted to Matlab based on initial work by Joel Fiddes and routines in R and Python. 

The TopoLAB routines have been tested in several sites including: Svalbard, Finnmark, Hokkaido, the Sierra Nevada, Finse, Mongolia, and the Pyrenees.

## General recipe for running TopoSCALE

1. Download ERA5 data using python scripts, then merge the .nc files using the Matlab script get_ERA5.m, resulting in a sturcutre ERA5_merged.mat
2. Download the required SRTM (or ArcticDEM, or ASTER GDEM, or local) DEM tiles for your area of interest. Then run the matlab script get_SRTM.m to reproject (to UTM x,y coordinates), resample, and crop out a square padded domain (set pad accordingly) centered on your area of interest. You create your area of interest using the "mask" array (this can be a square area, or any arbitrary shape, as well as scattered points in the domain). TopoPAR and TopoSCALE will only run for points where mask=1. You can also define an arbitrary target spatial resolution when you resample the DEM, through first nearest neighbor resampling then spatial aggregation (averaging). 
3. Run TopoPAR (see the function TopoPAR.m) using get_TopoPAR.m, this calculates the terrain parameters (elevation, sky view factor, aspect, slope, and horizon angles in 36 bins) for the points in your domain. 
4. To ease the compuational burden, especially for large (>1000s of grid cells) domains consider running the TopoSUB routine. This will perform a dimensionality reduction on your topograhic parameters (from TopoPAR) by assigning them to clusters. Within each cluster, the different parameters should be relatively similar compared to across the clusters. The number of clusters is (by design) several orders of magnitude lower than the number of grid cells. Due to the within cluster similarity running TopoSCALE for the clusters rather than the full DEM will lead to only a minor reduction in performance, but with huge gains in run time and storage. 
4. Run TopoSCALE (see the function TopoSCALE.m) using get\_TopoSCALE.m. Here you can set the time that you are interested in downscaling for. Note that this time vector (t) must be the completely contained within your ERA5 data (ERA5_merged.mat) and also be at the same resolution (typically hourly). 

**References**

You can read more about the routines making up TopoLAB (TopoSCALE, TopoSUB, TopoPAR) are described in the following publications:

- Fiddes, J. and Gruber, S.: TopoSCALE v.1.0: downscaling gridded climate data in complex terrain, Geosci. Model Dev., 7, 387–405, https://doi.org/10.5194/gmd-7-387-2014, 2014.
- Fiddes, J. and Gruber, S.: TopoSUB: a tool for efficient large area numerical modelling in complex topography at sub-grid scales, Geosci. Model Dev., 5, 1245–1257, https://doi.org/10.5194/gmd-5-1245-2012, 2012.
- Dozier, J. and Frew, J.: Rapid calculation of terrain parameters for radiation modeling from digital elevation data, IEEE TGRS, https://doi.org/10.1109/36.58986, 1990

The equations and workflow in TopoLAB is also outlined more extensively in Section 2.1 (Forcing) of my thesis:

- Aalstad, K.: Ensemble-based retrospective analysis of the seasonal snowpack, PhD Thesis at the University of Oslo, http://dx.doi.org/10.13140/RG.2.2.11709.20961, 2019. 






