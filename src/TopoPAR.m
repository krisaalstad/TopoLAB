function [ tp ] = TopoPAR_par( DEM )
%% TopoPAR_par.m (parallel version of TopoPAR)
% Calculates topographic parameters, including solar geometry, for a
% given region and time period as specified by an input time array and
% an input digital elevation model (DEM) structure.
% Inputs: 
%           t = [Nt x 1] time array. Time is in UTC (Zulu time) as the
%           number of days since 01-Jan-0000 + 1. This is the standard
%           MATLAB datenum format.
%           DEM = A DEM structure. This must contain a 'utmz' field
%           specifiying the UTM zone, a 'domain' field that contains the UTM
%           coordinate arrays  (x,y) and elevation (z) for the raster
%           domain that the topographic parameters will be calculated over,
%           and a 'pad' field that comtains the same arrays as the 'domain'
%           field but for a larger padded domain that encompasses the true
%           domain.
% Outputs:
%           tp = A topographic parameter structure that contains the
%           aspect, slope, elevation, solar zenith angle (time series),
%           solar azimuth angle (time series), horizon angle and sky view
%           factor for the raster domain.
% References:
% Zaksek et al. (2011), doi:10.3390/rs3020398 
% Dozier and Frew (1990), doi:10.1109/36.58986
%
% Code by Kristoffer Aalstad (March 7, 2018).
% Parallel version October 24, 2021.

%% Related to parallelization
ncores=2;%20;
Nc=20;%sum(DEM.mask(:));%100; % Number of chunks (*must be a multiple of ncores). 
% Chunk size of Np=sum(DEM.mask(:)) is the same as looping over all cells,
% but with parfor overhead. So best is some order of magnitude less than
% Np, but greater than ncores.


%% Extract DEM data from the DEM structure.

mask=DEM.mask;
x=DEM.x;
y=DEM.y;%
z=DEM.Z;
[X,Y]=meshgrid(x,y);
del=x(2)-x(1); % dx=dy


%% Elevation derivatives. N.B. not computed for the boundary of the DEM.
dzdy=NaN.*z;  dzdx=dzdy;
dzdy(2:end-1,2:end-1)=-1.*(z(3:end,2:end-1)-z(1:end-2,2:end-1))./(2.*del); % N.B. y decreases with i.
dzdx(2:end-1,2:end-1)=(z(2:end-1,3:end)-z(2:end-1,1:end-2))./(2.*del);

%% Slope and aspect. 
% Aspect definition: N,E,S,W = pi, pi/2, 0, -pi/2
asp=atan2(-1.*dzdx,dzdy);
slp=atan(sqrt(dzdx.^2+dzdy.^2));
asp(~mask)=NaN; slp(~mask)=NaN;

%% Horizon angles.

% Quite slow due to detailed computation, but you only have to do this
% once per site.
bins=0:10:350; nbins=numel(bins);
ii=numel(y); jj=numel(x);
h=255.*uint8(ones(ii,jj,nbins));
pcts=1:1:1e2;

% Only loop over valid points (mask=1).
is=1:size(mask,1); is=is';
js=1:size(mask,2); js=js';
[Js,Is]=meshgrid(js,is);
Np=sum(mask(:));
Ism=Is(mask);
Jsm=Js(mask);
nsm=sub2ind(size(mask),Ism,Jsm);
xm=X(nsm); ym=Y(nsm);
hv=zeros(Np,nbins);
svfv=zeros(Np,1);
c=parcluster('local');
par=gcp('nocreate');
if isempty(par)
    delete(par);
    parpool('local',ncores);
else
    if ~par.Connected||par.NumWorkers~=ncores
        delete(par);
        parpool('local',ncores);
    end
end
%tic;

chunksize=ceil(Np/Nc);
disp('Starting horizon computations');
parfor_progress(Nc);
parfor n=1:Nc
    these=((n-1)*chunksize+1):min(n*chunksize,Np);
    [s(n).h,s(n).svf]=horizon_chunk(xm(these),ym(these),x,y,z);
    parfor_progress;
end
% Merge chunks
for n=1:Nc
    these=((n-1)*chunksize+1):min(n*chunksize,Np);
    hv(these,:)=s(n).h;
    svfv(these,:)=s(n).svf;
end
clear s;
parfor_progress(0);
SVF=zeros(size(mask));
SVF(mask)=svfv;
for nb=1:nbins
    tmp=nan(size(mask));
    tmp(mask)=hv(:,nb);
    h(:,:,nb)=tmp;
end



% Save the center of the coarse bins. Convert these to the same coordinate
% system as the aspect. 
hbins=bins+5; % This is maybe bad, hard coded increment?
hbins=deg2rad(hbins);
hbins=(5*pi/2)-hbins;
hbins=hbins-2.*pi.*(hbins>2.*pi);
hbins=hbins+pi/2;
hbins=hbins-2.*pi.*(hbins>pi);



%% Create a structure of topographic parameters as output. 
tp.z=z(mask);
Np=numel(tp.z);
tp.asp=asp(mask);
tp.slp=slp(mask);
tp.h=zeros(Np,numel(hbins));
for j=1:numel(hbins)
    tmp=h(:,:,j); 
    tmp=tmp(mask);
    tp.h(:,j)=tmp;
end
tp.h=deg2rad(tp.h);
tp.svf=SVF(mask);
[X,Y]=meshgrid(x,y);
tp.x=X(mask); % Masked coordinates.
tp.y=Y(mask);
tp.xg=x; % To reconstruct the whole grid
tp.yg=y;
tp.mask=mask;
tp.utmz=DEM.utmz;
tp.hbins=hbins';


end

