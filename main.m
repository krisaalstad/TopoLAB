%% Main script demonstrating how to use the different components of TopoLAB.
clearvars;
paths;
load('batlow.mat'); % Nice "scientific colormap"

%% Pick a study site for this demonstration
% Options: 'Brøgger' (on Svalbard) or 'Terelj' (in Mongolia)
site='Brøgger';
% Note, Terelj is a smaller area so it will run faster.

%% Pre-processing fetching input

% Fetch the input data if you haven't yet.
if ~exist('input','dir')
    disp('Fetching input data');
    urlis='https://www.dropbox.com/s/9qckrp6tp7v6a81/input.zip?dl=1'; % Changed dl=0 to dl=1
    tarf='input.zip';
    websave(tarf,urlis);
    unzip(tarf,''); % Unzip
    paths; % Update path
end

if any(strcmp(site,'Brøgger'))
    demfile='NP_DEM_Brøgger.mat';
elseif any(strcmp(site,'Terelj'))
    demfile='SRTM_Terelj.mat';
else
    error('Please specify either Brøgger or Terelj as the study site');
end
tardir=[site '/DEM/'];
tarfile=[tardir demfile];
load(tarfile);
if isa(DEM.Z,'uint16') % Convert Z to double if necessary
    DEM.Z=double(DEM.Z).*DEM.Z_sf;
end
% See the "get_DEM" scripts for examples of how to read in and pre-process
% DEMs for TopoLAB.


%% Get the TopoPARs 
doTopoPAR=false; % Switched off by default since it's slow
% The topographic parameters have been pre-computed for these sites.

if doTopoPAR
    ncores=2; % Set the number of cores (assumes you have at least 2).
    tp=TopoPAR(DEM,ncores);
else
    tarfile=[tardir 'TopoPARs.mat'];
    load(tarfile);
end


%% Run TopoSUB
Nc=2e2; % Number of clusters
maxit=1e4; % Maximum number of iterations for k-means
doxy=1; % Include x,y coordinates as predictors (1=true, 0=false)
disp('Running TopoSUB');
tic;
[cn,indsc]=TopoSUB(tp,Nc,maxit,doxy);
toc;


%% Visualize the TopoSub clustering
figure(1); clf;
xl=[min(tp.x) max(tp.x)]; dxl=xl(2)-xl(1);
yl=[min(tp.y) max(tp.y)]; dyl=yl(2)-yl(1);

ax1=subplot(1,3,1);
img=nan(size(tp.mask));
img(tp.mask)=tp.z;
imagesc(tp.xg,tp.yg,img,'AlphaData',tp.mask);
axis xy;
xlim(xl); ylim(yl);
pbaspect([1 dyl/dxl 1]);
title('DEM elevation $z$ [m asl]','Interpreter','Latex','FontSize',16);
set(gca,'XTick',[],'YTick',[]);
colormap(ax1,batlow);

ax2=subplot(1,3,2);
img=nan(size(tp.mask));
img(tp.mask)=cn;
imagesc(tp.xg,tp.yg,img,'AlphaData',tp.mask);
axis xy;
xlim(xl); ylim(yl);
pbaspect([1 dyl/dxl 1]);
title('TopoSUB clusters','Interpreter','Latex','FontSize',16);
set(gca,'XTick',[],'YTick',[]);
colormap(ax2,'hsv');


ax3=subplot(1,3,3);
img=nan(size(tp.mask));
zc=tp.z(indsc); % Elevation of the sample centroids for each cluster
zcs=zc(cn); % Assign these to each point in the respective clusters.
img(tp.mask)=zcs;
imagesc(tp.xg,tp.yg,img,'AlphaData',tp.mask);
axis xy;
xlim(xl); ylim(yl);
pbaspect([1 dyl/dxl 1]);
title('TopoSUB elevation $z$ [m asl]','Interpreter','Latex','FontSize',16);
set(gca,'XTick',[],'YTick',[]);
colormap(ax3,batlow);

%% Make a new topographic parameter structure using the TopoSUB output
tps=tp;
tps.cn=cn;
tps.indsc=indsc;
tps.z=tp.z(indsc);
tps.asp=tp.asp(indsc);
tps.slp=tp.slp(indsc);
tps.svf=tp.svf(indsc);
tps.h=tp.h(indsc,:);
tps.x=tp.x(indsc);
tps.y=tp.y(indsc);
tps.xg=tp.xg;
tps.yg=tp.yg;
tps.mask=tp.mask;
tps.utmz=tp.utmz;
tps.hbinds=tp.hbins;

%% Read in the forcing data.

% The ERA5 data has already been processed, check out the
% input/site/Forcing/get_ERA5.m script for details.
tarfile=[site '/Forcing/ERA5_merged.mat'];
load(tarfile); % Load the "era" input forcing structure.

%% Run TopoSCALE

years=2019;
% Specify the time vector (only 2019 data is available in this example)
t=datenum(sprintf('01-01-%d',min(years))):(1/24):datenum(sprintf('01-01-%d',max(years)+1)); % Target time steps.
t=t(1:(end-1)); % Cut out the last timestep which is in 2020
t=t';

% Switches:
nosurf=1; % Exclude ERA5 surface levels from interpolation (recommended)
scaleout=1; % Scale outputs for compression (recommended)

% Call TopoSCALE
disp('Running TopoSCALE');
tic;
f=TopoSCALE(t,tps,era,nosurf,scaleout);
toc;

%% Visualize some of the outputs.
load('batlow.mat'); % Load a nice colormap.

figure(2); clf;

ax1=subplot(2,3,1);
img=nan(size(tp.mask));
% Extract shortwave radiation (SW), remembering to decompress with the scale 
% factor which in the case of radiation varibles is f.rad_sf 
SW=double(f.SW).*f.rad_sf;
SWi=sum(SW,2).*(60^2/1e6); % Convert to mega Joules (MJ) of annual insolation.
SWis=SWi(cn); % Spatialize (assign clusters value to grid cells).
img(tp.mask)=SWis;
imagesc(tp.xg,tp.yg,img,'AlphaData',tp.mask);
axis xy;
xlim(xl); ylim(yl);
pbaspect([1 dyl/dxl 1]);
%title('Annual incoming shortwave','Interpreter','Latex','FontSize',16);
set(gca,'XTick',[],'YTick',[],'LineWidth',1.5);
colormap(ax1,batlow);
c=colorbar;
c.TickLabelInterpreter='latex';
c.TickDirection='out';
c.FontSize=10;
c.LineWidth=1.5;
caxis([min(img(:)) max(img(:))]);
xlabel(c,'Shortwave, $\int S^\downarrow \,  dt$ [MJm$^{-2}$]',...
    'Interpreter','Latex','FontSize',16);

ax2=subplot(2,3,2);
img=nan(size(tp.mask));
% Extract longwave radiation (LW), remembering to decompress with the scale 
% factor which in the case of radiation varibles is f.rad_sf 
LW=double(f.LW).*f.rad_sf;
LWi=sum(LW,2).*(60^2/1e6); % Convert to mega Joules (MJ) of annual insolation.
LWis=LWi(cn); % Spatialize (assign clusters value to grid cells).
img(tp.mask)=LWis;
imagesc(tp.xg,tp.yg,img,'AlphaData',tp.mask);
axis xy;
xlim(xl); ylim(yl);
pbaspect([1 dyl/dxl 1]);
%title('Annual incoming longwave','Interpreter','Latex','FontSize',16);
set(gca,'XTick',[],'YTick',[],'LineWidth',1.5);
colormap(ax2,batlow);
c=colorbar;
c.TickLabelInterpreter='latex';
c.TickDirection='out';
c.FontSize=10;
c.LineWidth=1.5;
caxis([min(img(:)) max(img(:))]);
xlabel(c,'Longwave, $\int L^\downarrow \,  dt$ [MJm$^{-2}$]',...
    'Interpreter','Latex','FontSize',16);

ax2=subplot(2,3,3);
img=nan(size(tp.mask));
% Extract shortwave radiation (SW), remembering to decompress with the scale 
% factor which in the case of radiation varibles is f.rad_sf 
T=double(f.T).*f.T_sf;
Tm=mean(T,2); % Calculate the annual mean air temperature.
Tms=Tm(cn); % Spatialize (assign clusters value to grid cells).
img(tp.mask)=Tms;
imagesc(tp.xg,tp.yg,img,'AlphaData',tp.mask);
axis xy;
xlim(xl); ylim(yl);
pbaspect([1 dyl/dxl 1]);
%title('Annual mean air temperature','Interpreter','Latex','FontSize',16);
set(gca,'XTick',[],'YTick',[],'LineWidth',1.5);
colormap(ax2,batlow);
c=colorbar;
c.TickLabelInterpreter='latex';
c.TickDirection='out';
c.FontSize=10;
c.LineWidth=1.5;
caxis([min(img(:)) max(img(:))]);
xlabel(c,'Mean air temperature, $\overline{T}$ [$^\circ$C]',...
    'Interpreter','Latex','FontSize',16);



subplot(2,3,4:6);
% Plot the air temperature for min and max elevation
here=tps.z==max(tps.z)|tps.z==min(tps.z);
zex=round(tps.z(here));
T=double(f.T(here,:)).*f.T_sf;
yl=[min(T(:)) max(T(:))];
leglab=cell(numel(zex),1);
plot(f.t,zeros(size(f.t)),'Color',[0 0 0]); % 0 line
hold on;
cs=[0.8 0 0; 0 0 0.8];
for j=1:numel(zex)
    pt(j)=plot(f.t,T(j,:),'LineWidth',2,'Color',cs(j,:)); 
    leglab{j}=sprintf('z=%d [m asl]',zex(j));
end
xlim([min(t) max(t)]);
ylim(yl);
leg=legend(pt,leglab,'Interpreter','Latex','FontSize',16);
set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.005, 0.005]);
set(groot, 'defaultAxesTickLabelInterpreter','Latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',18);
ylabel('Air temperature, $T$ [$^\circ$C]','Interpreter','Latex','FontSize',18);
datetick('x','keeplimits','keepticks');
if strcmp(site,'Brøgger')
    sites='Br\o gger (Svalbard)'; % Make Brøgger latex compatible
else
    sites='Terelj (Mongolia)';
end
title(sprintf('ERA5 data downscaled with TopoSUB and TopoSCALE for %s in 2019',sites),...
    'Interpreter','Latex','FontSize',18);


print('-djpeg','tmp','-r300','-opengl');
