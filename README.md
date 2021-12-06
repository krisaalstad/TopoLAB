# TopoLAB
Topography-based routines for generating terrain parameters (TopoPAR) as well as downscaling (TopoSCALE) and clustering (TopoSUB) meteorological forcing data. These routines have been ported and adapted to Matlab based on initial work by Joel Fiddes and routines in R and Python. 

To test the TopoLAB routines simply run the main.m script which goes through the steps of calculating terrain parameters (TopoPAR), splitting the domain into clusters (TopoSUB), and running TopoSCALE. Note that this requires downloading some input data (roughly 70 MB) in the form of processed DEMs and ERA5 data. The two example sites are the Brøgger peninsula (the area around Ny-Ålesund) on Svalbard and the area aroudn Terelj in central Mongolia.  

These TopoLAB routines have been tested in several sites including: Svalbard, Finnmark, Hokkaido, the Sierra Nevada, Finse, Mongolia, and the Pyrenees. The original TopoSUB and TopoSCALE code that TopoLAB is based on has been tested extensively in the Swiss Alps as well as in central Asia. 

Contributors (both science & code): Kristoffer Aalstad, Joel Fiddes, Simon Filhol, Sebastian Westermann.

## General recipe for running TopoSCALE

1. Download ERA5 data using python scripts, then merge the .nc files using the Matlab script get_ERA5.m, resulting in a sturcutre ERA5_merged.mat
2. Download the required SRTM (or ArcticDEM, or ASTER GDEM, or local) DEM tiles for your area of interest. Then run the matlab script get_SRTM.m to reproject (to UTM x,y coordinates), resample, and crop out a square padded domain (set pad accordingly) centered on your area of interest. You create your area of interest using the "mask" array (this can be a square area, or any arbitrary shape, as well as scattered points in the domain). TopoPAR and TopoSCALE will only run for points where mask=1. You can also define an arbitrary target spatial resolution when you resample the DEM, through first nearest neighbor resampling then spatial aggregation (averaging). 
3. Run TopoPAR (see the function TopoPAR.m) using get_TopoPAR.m, this calculates the terrain parameters (elevation, sky view factor, aspect, slope, and horizon angles in 36 bins) for the points in your domain. 
4. To ease the compuational burden, especially for large (>1000s of grid cells) domains consider running the TopoSUB routine. This will perform a dimensionality reduction on your topograhic parameters (from TopoPAR) by assigning them to clusters. Within each cluster, the different parameters should be relatively similar compared to across the clusters. The number of clusters is (by design) several orders of magnitude lower than the number of grid cells. Due to the within cluster similarity running TopoSCALE for the clusters rather than the full DEM will lead to only a minor reduction in performance, but with huge gains in run time and storage. 
4. Run TopoSCALE (see the function TopoSCALE.m) using get\_TopoSCALE.m. Here you can set the time that you are interested in downscaling for. Note that this time vector (t) must be the completely contained within your ERA5 data (ERA5_merged.mat) and also be at the same resolution (typically hourly). 

## References

You can read more about the routines making up TopoLAB (TopoSCALE, TopoSUB, TopoPAR) in the following publications:

- Fiddes, J. and Gruber, S.: TopoSCALE v.1.0: downscaling gridded climate data in complex terrain, Geosci. Model Dev., 7, 387–405, https://doi.org/10.5194/gmd-7-387-2014, 2014.
- Fiddes, J. and Gruber, S.: TopoSUB: a tool for efficient large area numerical modelling in complex topography at sub-grid scales, Geosci. Model Dev., 5, 1245–1257, https://doi.org/10.5194/gmd-5-1245-2012, 2012.
- Dozier, J. and Frew, J.: Rapid calculation of terrain parameters for radiation modeling from digital elevation data, IEEE TGRS, https://doi.org/10.1109/36.58986, 1990

In addition, most of the equations in TopoLAB are outlined in Section 2.1 of my thesis:

- Aalstad, K.: Ensemble-based retrospective analysis of the seasonal snowpack, PhD thesis at the University of Oslo, http://dx.doi.org/10.13140/RG.2.2.11709.20961, 2019. 

Some examples of recent studies applying TopoSCALE and TopoSUB include:

- Fiddes, J., Aalstad, K., and Westermann, S.: Hyper-resolution ensemble-based snow reanalysis in mountain regions using clustering, Hydrol. Earth Syst. Sci., 23, 4717–4736, https://doi.org/10.5194/hess-23-4717-2019, 2019.

- Fiddes, J., Aalstad, K., and Lehning, M.: TopoCLIM: Rapid topography-based downscaling of regional climate model output in complex terrain v.1.0, Geosci. Model Dev. Discuss. [preprint], https://doi.org/10.5194/gmd-2021-60, in review, 2021.

with more on the way.

## Related repositories

- https://github.com/joelfiddes/topoMAPP : Modeling framework combining TopoSUB and TopoSCALE with ensemble-based data assimilation methods (Fiddes et al., 2019).

- https://github.com/joelfiddes/toposcale : The original TopoSCALE repository (Fiddes and Gruber, 2014). 

- https://github.com/joelfiddes/toposub : The original TopoSUB directory (Fiddes and Gruber, 2012). 

with more in the pipeline.

## Acknowledgements

A big thanks to ECMWF and the Copernicus Climate Change Service (C3S) for providing the ERA5 reanalysis data (Hersbach et al. 2020). Thansk as well to the agencies that provide access to global or regional DEMs such as SRTM (USGS), ASTER GDEM (NASA), the ArcticDEM (Polar Geospatial Center), and DEMs of Svalbard (Norwegian Polar Institute). Thanks also to Crameri et al. (2020) for providing some nice colormaps.  

- Crameri, F., G.E. Shephard, and P.J. Heron (2020), The misuse of colour in science communication, Nature Communications, 11, 5444. doi:10.1038/s41467-020-19160-7

- Hersbach, H, Bell, B, Berrisford, P, et al. The ERA5 global reanalysis. Q J R Meteorol Soc. 2020; 146: 1999– 2049. https://doi.org/10.1002/qj.3803





