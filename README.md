# TopoLAB
Topography-based routines for generating terrain parameters (TopoPAR) as well as downscaling (TopoSCALE) and clustering (TopoSUB) meteorological forcing data. These routines have been ported and adapted to Matlab based on initial work by Joel Fiddes and routines in R and Python. 

The TopoLAB routines have been tested in several sites including: Svalbard, Finnmark, Hokkaido, the Sierra Nevada, Finse, Mongolia, and the Pyrenees.



## General recipe for running TopoSCALE

1. Download ERA5 data using python scripts, then merge the .nc files using the Matlab script get_ERA5.m, resulting in a sturcutre ERA5_merged.mat
2. Download the required SRTM (or ArcticDEM, or ASTER GDEM, or local) DEM tiles for your area of interest. Then run the matlab script get_SRTM.m to reproject (to UTM x,y coordinates), resample, and crop out a square padded domain (set pad accordingly) centered on your area of interest. You create your area of interest using the "mask" array (this can be a square area, or any arbitrary shape, as well as scattered points in the domain). TopoPAR and TopoSCALE will only run for points where mask=1. You can also define an arbitrary target spatial resolution when you resample the DEM, through first nearest neighbor resampling then spatial aggregation (averaging). 
3. Run TopoPAR (see the function TopoPAR.m) using get_TopoPAR.m, this calculates the terrain parameters (elevation, sky view factor, aspect, slope, and horizon angles in 36 bins) for the points in your domain. 
4. To ease the compuational burden, especially for large (>1000s of grid cells) domains consider running the TopoSUB routine. This will perform a dimensionality reduction on your topograhic parameters (from TopoPAR) by assigning them to clusters. Within each cluster, the different parameters should be relatively similar compared to across the clusters. The number of clusters is (by design) several orders of magnitude lower than the number of grid cells. Due to the within cluster similarity running TopoSCALE for the clusters rather than the full DEM will lead to only a minor reduction in performance, but with huge gains in run time and storage. 
4. Run TopoSCALE (see the function TopoSCALE.m) using get_TopoSCALE. Here you can set the time that you are interested in downscaling for. Not that this time vector (t) must be the completely contained within your ERA5 data (ERA5_merged.mat) and also be at the same resolution (typically hourly). The results are stored in "TopoSCALE.mat". There is an example of how to visualize the outputs. All the forcing variables have dimensions Number of spatial points (grid cells) x Number of time steps. Not that these are stored in a compressed format in the "f" structure. To decompress, for air temperature simply do 
T=double(f.T).*f.T_sf 
which converts to double and multiplies by the relevant scale. Shortwave radiation would be
SW=double(f.SW).*f.rad_sf
and so on for the other variables (just use the corresponding factor). 
Units: T [K], SW [Wm^-2], LW [Wm^-2], q [kg/kg], U [ms^-1], P [mm/hour], p [Pa].
To convert from a list of points to a grid (raster) for visualization, you can do the following for mean air temperature for example:
I=nan(size(f.mask));
T=double(f.T).*f.T_sf;
Tmean=mean(T,2);
I(f.mask)=Tmean;
figure(1); clf;
imagesc(f.xg,f.yg,I,'AlphaData',f.mask);
axis xy;
