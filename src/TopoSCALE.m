function [ fout ] = TopoSCALE( t, tp, fin, nosurf, scaleout)
%% TopoSCALE.m
% [ fout ] = TopoSCALE( t, tp, fin, nosurf, scaleout)
% A MATLAB implementation of the TopoSCALE routine for topography-based
% downscaling of gridded meteorological data (reanalyses, climate model
% outputs, weather forecasts, etc...).
% Modified to interpolate surface fields as well.
% This version is built around downscaling ERA5 (or ERA-Interim) data, but
% is readily adapted to other products as well.
% Dimensions: Nt = Number of timesteps, Nf = Number of fine (target) grid
% cells, Nc = Number of coarse horizontal (lon,lat) grid cells, Np = Number
% of pressure levels.
% Inputs:
%           t = (Nt x 1) time array. Time is in UTC (Zulu time) as the
%               number of days since 01-Jan-0000 + 1. This is the standard
%               MATLAB datenum format.
%           tp = A topographic parameter structure. This must contain a 'utmz' field
%           specifiying the UTM zone, as well as the following arrays
%               tp.lat = Latitudes (Nf x 1)*
%               tp.lon = Longitudes (Nf x 1)*
%               tp.utmz = Utm zone string (e.g. '32V')*
%               tp.x = Easting coordinates [m] in UTM zone (Nf x 1)*
%               tp.y = Northing coordinates (Nf x 1)
%               tp.asp = Aspect [radians*] (Nf x 1)
%               tp.slp = Slope [radians] (Nf x 1)
%               tp.svf = Sky view factor [-] (Nf x 1)
%               tp.h = Horizon angle [radians] (Nf x Nh)
%               tp.hbins = Horizon angle bins [radians**] (Nh x 1)
%                   *: Coordinates can be specified either in lat,lon or as
%                   utm coordinates.
%                   *: Azimuthal quantities are defined in radians,
%                   clockwise from south i.e.: [N,E,S,W]=[pi,pi/2,0,-pi/2]
%           fin = A coarse scale input forcing data structure. Contains:
%               lat = 1D latitude array.
%               lon = 1D longitude array.
%               p = 1D array containing pressure levels [Pa].
%               Zs =  2D surface elevation array [m].
%               t = 1D timestamp array (datenum format).
%               SW = 3D downwelling shortwave array [W/m^2].
%               LW = 3D downwelling longwave array [W/m^2].
%               P = 3D precipitation array [mm/hour].
%               Td2 = 3D dewpoint temperature array [Kelvin].
%               T2 = 3D 2-meter temperature array [Kelvin].
%               u10 = 3D 10-meter zonal wind [m/s].
%               v10 = 3D 10-meter meridional wind [m/s].
%               ps = 3D surface pressure array [Pa].
%               u = 4D u-wind component array [m/s]
%               v = 4D v-wind component array [m/s]
%               T = 4D temperature array [Kelvin].
%               q = 4D specific humidity array [-].
%               Z = 4D geopotential height of pressure levels array [m] (geopotential/g).
%               Array dimensions: 2D = lon x lat, 3D = lon x lat x time, 4D = lon x lat x pressure x time
%           nosurf = Binary switch to not include (1=don't include,
%                    0=include) the surface level in vertical
%                    interpolations. Recommended value=1 (don't include).
%           scaleout = Binary switch (1=scale, 0=don't scale) to scale
%                    outputs and compress them to lower precision than
%                    double. This will reduce the size of the fout struct
%                    considerably without affecting accuracy considerably.
%                    Recommended value=1 if you are saving the output,
%                    otherwise 0 is fine if you will work wtih the
%                    variables directly in memory.
%
%
% Outputs:
%           fout = A downscaled output forcing data structure.
% Requires:
%           The tp structure can be obtained following a call to the
%           "TopoPAR.m" routine using a DEM. The input forcing structure
%            can, e.g., be obstained following a call to the "get_ERA.m"
%            routine which processes "raw" ERA (5 or Interim) data.
%
% References:
% Fiddes and Gruber (2014), doi:10.5194/gmd-7-387-2014
% Fiddes, Aalstad, and Westermann (2019), doi:10.5194/hess-23-4717-2019
% Aalstad (2019), PhD Thesis, doi:10.13140/RG.2.2.11709.20961
% Dozier and Frew (1990), doi:10.1109/36.58986
%
% Code by Kristoffer Aalstad (March, 2018). Last revision March 2021.
%
% Output variables units are the same as the input, unless compressed where:
%       u,v = [cm/s] as int16. So to the nearest 0.01 ms^-1 and with
%           a maximum absolute value of 327 [ms^-1] [unphysical].
%       q   = [mg/kg] as uint16. So to the nearest 0.001 gkg^-1 and with a
%           maximum value of  65 gkg^-1 [unphysical].
%       ps  = [hPa] as uint16. So to the nearest hPa with a maximum value
%           of  65000 hPa (unphysical).
%     SW,LW = [dWm^-2] as uint16. So to the nearest 0.1 Wm^-2 and with a
%           maximum value of 6500 Wm^-2 [unphysical].
%       T   = [c\deg C] as int16. So to the nearest 0.01 \deg C and with
%           a maximum absolute value of 328 \deg C [unphysical].
%       P   = [cmm/hour] as uint16 so to the nearest 0.01 mm/hour and with
%           a maximum avalue of 655 mm/hour which is twice the world
%           record hourly rainfall rate.
%     All scale factors are saved to allow easy conversion to the typical
%           units: m/s, kg/kg, hPa, Wm^-2, \deg C, mm/hour.
%     So, e.g. to convert temperature to double and \deg C:
%           T=double(fout.T).*fout.T_sf
%
% Gridding:
%       The variables are saved as a 2d (space,time) array. To recover
%       the spatial positions of the variables use a mask array. E.g.:
%               Tscaled=double(fout.T(:,n)).*fout.T_sf;
%               T=nan(size(mask));
%               T(mask)=Tscaled;
%       The mask array is a 2D logical array where mask=1 means a grid cell
%       is in the model domain.
% NB this version is modified to take a list (1D array) of target points
% (fine grid cells) in the "tp" structure.

% Scale factors:
if scaleout==1
    fout.wind_sf=1e-2;
    fout.q_sf=1e-6;
    fout.ps_sf=1e2;
    fout.rad_sf=1e-1;
    fout.T_sf=1e-2;
    fout.P_sf=1e-2;
end


% Prepare the spatiotemporal referencing of the output struct "fout".
if size(t,2)>1 % Make sure the target time is a column vector.
    t=t';
end
fout.t=t;
utmstruct=defaultm('utm');
utmstruct.geoid=wgs84Ellipsoid;
% Calculate target utm coordinates if these don't exist in the tp struct.
if ~isfield(tp,'x')
    tp.utmz=utmzone(mean(tp.lat),mean(tp.lon));
    utmstruct.zone=tp.utmz;
    utmstruct=defaultm(utmstruct);
    [tp.x,tp.y]=minvtran(utmstruct,tp.lat,tp.lon);
else
    utmstruct.zone=tp.utmz;
    utmstruct=defaultm(utmstruct);
    % Calculate target lat,lon coordinates if these don't exist.
    if ~isfield(tp,'lat')
        [tp.lat,tp.lon]=minvtran(utmstruct,tp.x,tp.y);
    end
end
fout.x=tp.x;
fout.y=tp.y;
fout.utmz=tp.utmz;
fout.lat=tp.lat;
fout.lon=tp.lon;

% Get the solar geometry, which is assumed to be approximately invariant
% for the entire area of interest. This is relatively unproblematic
% even for a domain with an extent of 0.5 deg (500 km) or so.
[saz,sze]=solargeom(t,mean(tp.lat),mean(tp.lon));

% Calculate input utm coords based on the input lat lon grid.
% N.B. the ERA grid is (lon,lat) not (lat,lon)
[LAT,LON]=meshgrid(double(fin.lat),double(fin.lon));
[Xc,Yc]=mfwdtran(utmstruct,LAT,LON);

% Constants.
R=287.05;  % Gas constant for dry air [JK^-1kg^-1]
g=9.81; % Acceleration of gravity [ms^-1]
eps0=0.622; % Ratio of molecular weight of water and dry air [-]
S0=1370; % Solar constat (total TOA solar irradiance) [Wm^-2] used in ECMWF's IFS

% Fine "grid" (the target)
% This is a vector of (x,y,z,t) points with a certain spatial extent given
% by the resolution of the terrain parameters (z,aspect,slope).
Nt=numel(t);
Z=tp.z;
X=tp.x;
Nf=numel(X);
Y=tp.y;
asp=tp.asp;
slp=tp.slp;
svf=tp.svf;
h=tp.h;
hbins=tp.hbins;

% Repeat copies of the coarse grid for each pressure level (3D arrays).
ppl=fin.p;
Np=numel(ppl);
ppl=repmat(ppl,Nf,1); % Pressure levels for each fine grid cell (Nf x Np)

% Output structure:
temp=zeros(Nf,Nt,'double'); 
fout.T=temp; fout.U=temp;  fout.q=temp;
fout.LW=temp; fout.SW=temp; fout.P=temp; fout.ps=temp; 

% Calculate inverse distance weights (IDWs).
Nc=numel(Xc(:)); % Number of coarse grid cells (in the horizontal plane).
% Number of nearest neighbors used in inverse distance weighted interpolation
kNN=min(9,Nc);
%kNN=4;
%NN=knnsearch([Xc(:) Yc(:)],[X(:)  Y(:)],'k',kNN);
Xcr=repmat(Xc(:)',Nf,1); Ycr=repmat(Yc(:)',Nf,1);
Xfr=repmat(X,1,Nc); Yfr=repmat(Y,1,Nc);
d=sqrt((Xcr-Xfr).^2+(Ycr-Yfr).^2); % Nf x Nc
idw=1./(d.^2);
[~,rankis]=sort(d,2,'ascend');
idw(rankis>kNN)=0; % Set the weights other than the kNN to zero.
% Set non NN weights to 0.
idw=idw./(repmat(sum(idw,2),1,Nc)); % Nf x Nc
% Can replace IDW with Natural Neighbor at a later stage, would work in the
% same way.


% Height of the "grid" (coarse scale) surface (not pressure levels).
Zsl=double(fin.Zs); % Coarse-scale surface elevation.
Zsl=idw*Zsl(:); % Interpolated coarse-scale surface elevation to the fine grid.

% Beware that including surface levels (nosurf==1) can lead to unphysical
% extrapolations.

% Useful in-line functions.
K2C=@(Tk) Tk-273.15;
q2w=@(q) 0.5.*(1-sqrt(1-4.*q)); % Mixing ratio from specific humidity based on 2nd order Taylor series expansion.
wp2e=@(w,p) 0.5.*p.*(-1+sqrt(1+4.*w./eps0)); % Vapor pressure from mixing ratio based on 2nd order Taylor series expansion.
% This wp2e was bugged

% AERK from Alduchov and Eskridge (1996).
A1=17.625; B1=243.04; C1=610.94;
Magnus=@(tc) C1.*exp(A1.*tc./(B1+tc)); % A version of the Magnus formula with the AERK parameters.
% Note, e=Magnus(tdc) and es=Magnus(tc)

% Sync the time series (assumed to have the same time step).
these=(fin.t>=min(t))&(fin.t<=max(t));
fin.t=fin.t(these);
if numel(fin.t)~=Nt
    error('Missing input data for this time period');
end
fin.P=double(fin.P(:,:,these));
fin.SW=double(fin.SW(:,:,these));
fin.LW=double(fin.LW(:,:,these));
fin.ps=double(fin.ps(:,:,these));
fin.T2=double(fin.T2(:,:,these));
fin.u10=double(fin.u10(:,:,these));
fin.v10=double(fin.v10(:,:,these));
fin.Td2=double(fin.Td2(:,:,these));
fin.Z=double(fin.Z(:,:,:,these));
fin.q=double(fin.q(:,:,:,these));
fin.T=double(fin.T(:,:,:,these));
fin.u=double(fin.u(:,:,:,these));
fin.v=double(fin.v(:,:,:,these));

%fout.LW=uint16(fout.LW);

outstatus=0:10:100;
for n=1:Nt
    rpc=round(1e2*(n/Nt));
    here=outstatus==rpc;
    if any(here)
        outstatus(here)=NaN;
        fprintf('\n %d percent of toposcaling complete \n',rpc);
    end
    
    
    %% Identify the nearest levels for each grid point.
    
    % Create a column matrix (Nc x Np) containing the elevations of each
    % coarse grid cell (rows) for each pressure level (columns).
    Zpl=reshape(fin.Z(:,:,:,n),Nc,Np); % Nc x Np = (I x J) x Np
    %Znm=[Zn Zsc(:)]; % Merge surface elevations as final column.
    %Znm=idw*Znm; % Estimate surface elevation of each level at query points
    Zpl=idw*Zpl;
    % Define the effective coarse surface for vertical interpolation,
    % potentially ignoring it. 
    if nosurf==1
        Zsl=nan(size(Zsl)); % Ignore surface in vertical interpolation.
    end 
    Zm=[Zpl Zsl]; % Merge pressure and surface levels Nc x (Np + 1)
    dZ=Zm-repmat(Z,1,Np+1); % Difference between level elevations and query elevations
    adZ=abs(dZ); % Absolute difference.
    above=(dZ>=0);
    below=~above;
    adZabove=adZ; adZabove(below)=NaN;
    closestabove=min(adZabove,[],2); % Should always exist and be unique for each query point (row).
    overlying=(adZabove==closestabove);
    closest=overlying; % Store closest overlying level pointer, useful for later pressure interpolations.
    adZbelow=adZ; adZbelow(above)=NaN;
    closestbelow=min(adZbelow,[],2);
    underlying=(adZbelow==closestbelow);
    doextrap=~any(underlying,2); % Underlying may not exist, in that case extrapolate.
    if any(doextrap) % Check if extrapolation is necessary for some query points.
        % In that case, pick the two nearest overlying levels instead so
        % that we extrapolate down from the layer above.
        underlying(doextrap,:)=overlying(doextrap,:); % Set underlying to closest overlying level.
        adZabovex=adZabove(doextrap,:); % Only pick those query points you need to extrapolate.
        adZabovex(overlying(doextrap,:))=NaN; % Set old minima to NaN.
        nextabove=min(adZabovex,[],2);
        overlying(doextrap,:)=(adZabovex==nextabove); % Set overlying to next overlying level.
    end
    % Now we can pick out the elevations at the relevant ("underlying" and
    % "overlying") levels.
    %Z1=Zm;
    %Z1(~underlying)=0; % Would be nice to speed this up!
    %Z1=sum(Z1,2);
    Zmt=Zm';
    underlying=underlying';
    Z1=Zmt(underlying);
    
    %Z2=Zm;
    %Z2(~overlying)=0;
    %Z2=sum(Z2,2);
    overlying=overlying';
    Z2=Zmt(overlying);
    
    dZ21r=(Z2-Z1).^(-1); % Useful for vertical linear interpolation: 1/(Z2-Z1)
    dZ1=Z-Z1;
    vsc=dZ1.*dZ21r; % Use in the vertical interpolation, avoiding repeated calculation.
    
    % Seem to not be picking the right overlying and underlying levels.
    % Most likely because the logical indexing doesn't work as intended,
    % instead it follows column major order so that e.g. underlying which
    % is Ng x (Np+1) [Number of fine grid cells x (Number of plevs + 1)]
    % will when used in Z1=Zm(underlying) pick out elevations in the wrong
    % order, traversing rows first (down each column) across different fine
    % grid cells, rather than within fine grid cells. 
    % Solved easily by instead working with the transposed matrices.
    
    % For each variable, pick out the relevant levels for each query point.
    
    %% Temperature
    Tpl=reshape(fin.T(:,:,:,n),Nc,Np); Tpl=idw*Tpl;
    Tsl=reshape(fin.T2(:,:,n),Nc,1); Tsl=idw*Tsl;
    Tm=[Tpl Tsl];
    Tout=idw_vinterp(Tm,underlying,overlying,vsc);
    fout.T(:,n)=Tout;
    
    % Rather than just linearly interpolating, it would be much more elegant to
    % use MOST flux profile relationships with changing land surface
    % conditions; this is what is done to interpolate to the measurement
    % level anyways in most models, so we inherit the bias of the models'
    % surface conditions.
    % In this way, we would could just consider the large scale ("free
    % atmosphere") forcing for T, u, v, q
    
    % Huge extrapolations can happen due to surface effects. Want to avoid
    % this, so maybe don't use surface (2m) values?
    
    %% Wind speed
    % Zonal wind
    upl=reshape(fin.u(:,:,:,n),Nc,Np); upl=idw*upl;
    usl=reshape(fin.u10(:,:,n),Nc,1); usl=idw*usl; % May not be needed.
    um=[upl usl];
    uout=idw_vinterp(um,underlying,overlying,vsc);
    
    % Meridional wind
    vpl=reshape(fin.v(:,:,:,n),Nc,Np); vpl=idw*vpl;
    vsl=reshape(fin.u10(:,:,n),Nc,1); vsl=idw*vsl; % May not be needed.
    vm=[vpl vsl];
    vout=idw_vinterp(vm,underlying,overlying,vsc);
    
    % Compute horizontal wind speed
    U=sqrt(uout.^2+vout.^2);
    fout.U(:,n)=U;
    
    %% Humidity
    psl=reshape(fin.ps(:,:,n),Nc,1);
    Td2sl=reshape(fin.Td2(:,:,n),Nc,1);
    vpsl=Magnus(K2C(Td2sl));
    wsl=eps0.*vpsl./(psl-vpsl);
    qsl=wsl./(1+wsl); qsl=idw*qsl;
    qpl=reshape(fin.q(:,:,:,n),Nc,Np); qpl=idw*qpl;
    qm=[qpl qsl];
    qout=idw_vinterp(qm,underlying,overlying,vsc);
    fout.q(:,n)=qout;
    
    %% Pressure
    psl=idw*psl;
    pm=[ppl psl];
    %pclosest=pm;
    %pclosest(~closest)=0;
    %pclosest=sum(pclosest,2);
    %Zclosest=Zm;
    %Zclosest(~closest)=0;
    %Zclosest=sum(Zclosest,2);
    %Tclosest=Tm;
    %Tclosest(~closest)=0;
    %Tclosest=sum(Tclosest,2);
    %pclosest=pm(closest);
    %Zclosest=Zm(closest);
    %Tclosest=Tm(closest);
    closest=closest';
    pmt=pm';
    pclosest=pmt(closest);
    Zclosest=Zmt(closest);
    Tmt=Tm';
    Tclosest=Tmt(closest);
    
    Tbar=0.5.*(Tout+Tclosest);
    H=(R/g).*Tbar;
    dZ=Z-Zclosest;
    pout=pclosest.*exp(-dZ./H);
    fout.ps(:,n)=pout;
    
    %% Longwave routine.
    LWc=reshape(fin.LW(:,:,n),Nc,1);
    LWc=max(idw*LWc,0); % Interpolate to fine grid (but Zsl elevations).
    vpc=idw*vpsl; % Can also use qsl directly and convert this, but (maybe) more costly?
    wf=q2w(qout); % Convert to mixing ratio at fine grid.
    vpf=wp2e(wf,pout);
    
    % Use the vapor pressure and temperature to calculate clear sky
    % emssivity at grid and subgrid. [also function]
    x1=0.43; x2=5.7;
    cef=0.23+x1.*(vpf./Tout).^(1/x2); % Obs! got negative vapor pressure-> imaginary number in LW calc
    cec=0.23+x1.*(vpc./Tsl).^(1/x2);
    
    % Diagnose the all sky emissivity at grid.
    sbc=5.67e-8;
    aec=LWc./(sbc.*Tsl.^4);
    
    % Calculate the "cloud" emissivity at grid, assume this is the same at
    % subgrid.
    deltae=aec-cec;
    
    % Use the former cloud emissivity to compute the all sky emissivity at subgrid.
    aef=cef+deltae;
    LWf=aef.*sbc.*Tout.^4;
    
    % Scale LW with terrain configuration, considering both occlusion by
    % and emissions from the surrounding terrain. From Dozier & Frew 1990.
    LWf=svf.*LWf;
    %if rpc>50
    %    disp('stop!');
    %end
    if any(imag(LWf))
        error('something went wrong');
    end
    fout.LW(:,n)=LWf; % This step is really slow, no idea why (checked if conflicting data type).
    
    %tmp=svf.*LWf;
    %fout.LW(:,n)=LWf; % Why is this step slow?
    %fout.LW(:,n)=LWf;%+0.5.*(1+cos(slp)).*(1-svf).*0.99.*5.67e-8.*(273.15.^4);
    % Removed emissions from the surrounding terrain for now.
    
    
    %% Shortwave radiation
    % Estimate the average solar zenith and azimuth over the last time step
    % (if available). This is important, otherwise the (average) shortwave
    % radiation for a timestep (which is what we are estimating) can be quite biased.
    if n~=1
        sazn=0.5*(saz(n-1)+saz(n));
        szen=0.5*(sze(n-1)+sze(n));
    else
        sazn=saz(n);
        szen=sze(n);
    end
    
    % Compute downwelling TOA SW irradiance (i.e. the incoming shortwave
    % incident on a horizontal plane at the TOA), by accounting for the
    % solar zentih angle.
    mu0=max(cos(szen),0); % Trunacte negative values.
    % May also want to consider treating values mu0<0 for prominent topography
    % when the horizon  angles are less than 0.
    sunset=mu0<cosd(89);%(mu0==0); % Sunset switch.
    % Note, it would be better to use the true average ((1/tau) integral_t^(t+tau) mu0 dt)
    % but this approximation should be ok.
    SWtoa=S0.*mu0;
    
    % Get the surface shortwave from ERA5 and interpolate to fine grid.
    SWc=reshape(fin.SW(:,:,n),Nc,1);
    SWc=max(idw*SWc,0);
    
    if sunset % If the sun has "surely" set (mu0=cos of average zenith<=0).
        SWout=zeros(Nf,1); % For now, we ignore the small effect of topographic prominence.
    else
        kt=SWc./SWtoa; % Clearness index.
        kd=0.952-1.041.*exp(-1.*exp(2.3-4.702.*kt)); % Diffuse fraction.
        kd=max(kd,0);
        
        % Diffuse component.
        SWcdiff=kd.*SWc; % Diffuse shortwave radiation.
        SWfdiff=svf.*SWcdiff; % Scale with the sky-view factor to get the
        % Note that even if it is physically correct, the above sky view scaling
        % can bias models (by removing a compensating error source) that
        % don't account for reflections from the surrounding terrain.
        
        % Direct component
        SWcdir=SWc-SWcdiff;
        
       
        
        
        % Scale direct shortwave using Beer's law (see Aalstad 2019, Appendix A)
        ka=(g.*mu0./(psl)).*log(SWtoa./SWcdir); % Note, we don't get log(0) due to "if sunset" condition.
        SWfdir=SWtoa.*exp(-ka.*pout./(g*mu0));
        %SWfdir2=SWtoa.*((SWcdir./SWtoa).^(pout./psl)); % Combined but
        %slower expression.
        
        % Illumination angles.
        cosif=mu0.*cos(slp)+sin(szen).*sin(slp).*cos(sazn-asp); % Cosine of illumination angle at fine (target) grid.
        selfshadow=cosif<0; % Self shadowing, occurs when |saz-asp|>90
        cosif(selfshadow)=0;
        cosic=mu0; % Cosine of illumination angle at coarse grid (ERA5 assumes slope=0 deg)
        
        % Binary shadow masks.
        absdiff=abs(hbins-sazn); % To identify the relevant (minimum difference) horizon azimuth bin.
        [~,binis]=min(absdiff); % Index of the relevant bin, returns first optimum index if there are multiple.
        horizons=h(:,binis); % Identify the relevant horizon angles for each target grid cell.
        seln=max(pi/2-szen,0); % Solar elevation angle (in degrees for consistency with horizon angles).
        shade=horizons>seln; % True if a grid cell is in shadow.
        
        % Terrain corrected direct shortwave.
        SWfdir=SWfdir.*(cosif./cosic).*(1-shade); % Terrain corrected direct shortwave at subgrid.
        
        % Combine diff and dir shortwave.
        SWout=SWfdir+SWfdiff;
    end
    fout.SW(:,n)=SWout;
    
    %% Precipitation.
    Psl=reshape(fin.P(:,:,n),Nc,1);
    Pout=idw*Psl;
    fout.P(:,n)=Pout; % mm/h
end

fout.mask=tp.mask;
fout.xg=tp.xg;
fout.yg=tp.yg;

% Compress outputs if desired.
if scaleout==1
    fout.U=int16(fout.U./fout.wind_sf);
    fout.q=uint16(fout.q./fout.q_sf);
    fout.ps=uint16(fout.ps./fout.ps_sf);
    fout.SW=uint16(fout.SW./fout.rad_sf);
    fout.LW=uint16(fout.LW./fout.rad_sf); % Do this in the time loop.
    fout.T=int16((fout.T-273.15)./fout.T_sf);
    fout.P=uint16(fout.P./fout.P_sf);
end







    function fq=idw_vinterp(fm,these1,these2,vscale)
        % This function performs vertical interpolation between pressure levels.
        % It's just here to keep the code cleaner, since these operations
        % are called in an identical manner for many of the fields.
        %fm=[fp fs(:)]; % Merge pressure and surface fields.
        %fm=invdw*fm; % Perform inverse distance weighting, converting from Nc x Np to Nf x Np.
        %fm=[fp fs];
        %f1=fm(these1); % lower pressure level
        %f2=fm(these2); % upper pressure level
        %tmp=fm;
        %tmp(~these1)=0;
        %f1=sum(tmp,2);
        %tmp=fm;
        %tmp(~these2)=0;
        %f2=sum(tmp,2);
        fm=fm';
        f1=fm(these1);
        f2=fm(these2);
        
        fq=f1+(f2-f1).*vscale; % vertical (elevation based) interpolation between pressure levels
    end




    function [solar_azimuth,solar_zenith]=solargeom(Time,Latitude,Longitude)
        %% [saz,szen]=solargeom(time,latitude,longitude)
        % Adopted from the Sandia National Labs PVL Toolbox ephemeris routine.
        % Inputs:
        %   time = Time stamp vector (matlab datenum format) assumed to be in UTC
        %   latitude = Latitude
        %   longitude = Longitude
        % Outputs:
        %   saz = Solar azimuth angle [radians, anticlockwise from south]
        %   szen = Solar zentih angle [radians].
        % Link to the original toolbox:
        % https://pvpmc.sandia.gov/applications/pv_lib-toolbox/
        % References:
        % Stein et al. (2012), doi:10.1109/PVSC.2012.6318225 [MATLAB version]
        % Holmgren et al. (2018), doi:10.21105/joss.00884 [Python version]
        
        
        TZone=0;
        Longitude=-Longitude;
        tv=datevec(Time);
        Year=tv(:,1);
        v0=zeros(size(Year)); v1=ones(size(Year));
        DayOfYear=Time-datenum([Year v1 v1 v0 v0 v0])+1;
        DecHours=(Time - floor(Time)) .* 24;
        RadtoDeg=180/pi;
        DegtoRad=pi/180;
        Abber = 20/3600;
        LatR = Latitude * DegtoRad;
        UnivDate = DayOfYear + floor((DecHours + TZone)/24);
        UnivHr = mod((DecHours + TZone), 24);
        Yr = Year-1900;
        YrBegin = 365 * Yr + floor((Yr-1)/4)-0.5;
        Ezero = YrBegin + UnivDate;
        T = Ezero / 36525;
        GMST0 = 6/24 +38/1440 + (45.836 + 8640184.542 * T + 0.0929 * T.^2)/86400;
        GMST0 = 360 * (GMST0 - floor(GMST0));
        GMSTi = mod(GMST0 + 360*(1.0027379093 * UnivHr / 24),360);
        LocAST = mod((360 + GMSTi - Longitude), 360);
        EpochDate = Ezero + UnivHr / 24;
        T1 = EpochDate / 36525;
        ObliquityR = DegtoRad * (23.452294 - 0.0130125 * T1 - 0.00000164 * T1.^2 ...
            + 0.000000503 * T1.^3);
        MlPerigee = 281.22083 + 0.0000470684 * EpochDate + 0.000453 * T1 .^ 2 + ...
            0.000003 * T1 .^ 3;
        MeanAnom = mod((358.47583 + 0.985600267 * EpochDate - 0.00015 * T1 .^ 2 - ...
            0.000003 * T1 .^ 3), 360);
        Eccen = 0.01675104 - 0.0000418 * T1 - 0.000000126 * T1 .^ 2;
        EccenAnom = MeanAnom;
        E=0;
        while max(abs(EccenAnom - E)) > 0.0001
            E = EccenAnom;
            EccenAnom = MeanAnom + RadtoDeg .* Eccen .* sin(DegtoRad .* E);
        end
        TrueAnom = 2 * mod(RadtoDeg * atan2(((1 + Eccen) ./ (1 - Eccen)).^ 0.5 .* tan(DegtoRad * EccenAnom / 2), 1), 360) ;
        EcLon = mod(MlPerigee + TrueAnom, 360) - Abber ;
        EcLonR = DegtoRad * EcLon;
        DecR = asin(sin(ObliquityR) .* sin(EcLonR));
        %Dec = RadtoDeg * DecR;
        RtAscen = RadtoDeg * atan2(cos(ObliquityR).*(sin(EcLonR)),cos(EcLonR));
        HrAngle = LocAST - RtAscen ;
        HrAngleR = DegtoRad .* HrAngle ;
        %HrAngle = HrAngle - (360 .* sign(HrAngle) .* (abs(HrAngle) > 180));
        SunAz = RadtoDeg .* atan2(-1 * sin(HrAngleR), cos(LatR) .* tan(DecR) - sin(LatR) .* cos(HrAngleR));
        SunAz = SunAz + (SunAz < 0) * 360; %shift from range of [-180,180] to [0,360]
        SunEl = asind(cos(LatR) .* cos(DecR) .* cos(HrAngleR) + sin(LatR) .* sin(DecR));
        
        % Convert solar azimuth angle from [N,E,S,W]=[0,90,180,270] to [180, 90, 0
        % -90], i.e. the same as the aspect and horizon angle system.
        solar_azimuth=deg2rad(SunAz);
        solar_azimuth=(5*pi/2)-solar_azimuth;
        solar_azimuth=solar_azimuth-2.*pi.*(solar_azimuth>2.*pi);
        solar_azimuth=solar_azimuth+pi/2;
        solar_azimuth=solar_azimuth-2.*pi.*(solar_azimuth>pi);
        
        % Calculate solar zenith angle from solar elevation angle
        SunEl=deg2rad(SunEl);
        solar_zenith=(pi/2)-SunEl;
    end


end

