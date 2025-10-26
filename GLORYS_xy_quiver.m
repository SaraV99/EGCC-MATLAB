% GLORYS Velocity Quiver and Salinity Analysis
    % Evaluating the time-dependent flow at various depths in GLORYS to
    % determine when the EGCC splits from the EGC and hops onto the shelf
% Created by: Sara Vianco
% Masters Thesis
% 04 August 2023

%% Preliminaries
close all 
clear all

%% Setting working directory 
cd('/Users/saravianco/Documents/Research/MATLAB/GLORYS/GLORYS quiver/')
addpath('/Users/saravianco/Documents/Research/MATLAB/GLORYS/')
addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/') % for cmocean
addpath('/Users/saravianco/Documents/Research/MATLAB/m_map/');
addpath('/Users/saravianco/Documents/MATLAB//haversine/') % this is for haversine
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/html/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/library/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/thermodynamics_from_t/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/pdf/");
addpath('/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_quiver/');

% this was not preliminarily added until later when I went to compare the
% volumetric transport between RARE and GLORYS
load('date_num.mat')

%% Establishing variables
addpath('/Volumes/One Touch/Research Datasets/GLORYS/larger_data/')

% one-time variables
lat = ncread('GLORYS_091503_071507.nc','latitude');
lon = ncread('GLORYS_091503_071507.nc','longitude');
z = ncread('GLORYS_091503_071507.nc','depth');

% varying variables 09/15/2003 - 07/15/2007
time_1 = ncread("GLORYS_091503_071507.nc",'time');
v_1 = ncread("GLORYS_091503_071507.nc",'vo');
u_1 = ncread("GLORYS_091503_071507.nc",'uo');
tempp_1 = ncread("GLORYS_091503_071507.nc",'thetao');
psal_1 = ncread("GLORYS_091503_071507.nc",'so');

% 07/16/2007 - 07/16/2010
time_2 = ncread("GLORYS_071607_071610.nc",'time');
v_2 = ncread("GLORYS_071607_071610.nc",'vo');
u_2 = ncread("GLORYS_071607_071610.nc",'uo');
tempp_2 = ncread("GLORYS_071607_071610.nc",'thetao');
psal_2 = ncread("GLORYS_071607_071610.nc",'so');

% 07/17/2010 - 07/16/2013
time_3 = ncread("GLORYS_071710_071613.nc",'time');
v_3 = ncread("GLORYS_071710_071613.nc",'vo');
u_3 = ncread("GLORYS_071710_071613.nc",'uo');
tempp_3 = ncread("GLORYS_071710_071613.nc",'thetao');
psal_3 = ncread("GLORYS_071710_071613.nc",'so');

% 07/17/2013 - 07/16/2016
time_4 = ncread("GLORYS_071713_071616.nc",'time');
v_4 = ncread("GLORYS_071713_071616.nc",'vo');
u_4 = ncread("GLORYS_071713_071616.nc",'uo');
tempp_4 = ncread("GLORYS_071713_071616.nc",'thetao');
psal_4 = ncread("GLORYS_071713_071616.nc",'so');

% 07/17/2016 - 08/15/2019
time_5 = ncread("GLORYS_071716_081519.nc",'time');
v_5 = ncread("GLORYS_071716_081519.nc",'vo');
u_5 = ncread("GLORYS_071716_081519.nc",'uo');
tempp_5 = ncread("GLORYS_071716_081519.nc",'thetao');
psal_5 = ncread("GLORYS_071716_081519.nc",'so');

%% Concatenating into Single Variables
time = cat(1,time_1, time_2, time_3, time_4, time_5);
v = cat(4,v_1,v_2,v_3,v_4,v_5);
u = cat(4,u_1,u_2,u_3,u_4,u_5);
tempp = cat(4,tempp_1,tempp_2,tempp_3,tempp_4,tempp_5);
psal = cat(4,psal_1,psal_2,psal_3,psal_4,psal_5);

%% Adjusting the time
time_pivot = datenum(1950,1,1,0,0,0);
time_matlab = (time/24) + time_pivot;
time_matlab = floor(time_matlab);
date = datevec(time_matlab);

%% Establishing 3D Variables
% determining lat and lon limits for figure
latS_map = 73;
latN_map = 85;
lonW_map = -35;
lonE_map = 5;

lat2D_i = find(lat>latS_map & lat<latN_map);
lon2D_i = find(lon<lonE_map & lon>lonW_map);

lon2D = lon(lon2D_i);
lat2D = lat(lat2D_i);

%% Loading in bathymetry data
% Load 1-minute ETOPO data

addpath('/Users/saravianco/Documents/Research/MATLAB/ETOPO/')
info = ncinfo('ETOPO_2022_v1_60s_N90W180_bed.nc');

lonbath = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','lon');
latbath = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','lat');

idxSy = find(abs(latbath-latS_map)==nanmin(abs(latbath-latS_map)),1,'first');
idxNy = find(abs(latbath-latN_map)==nanmin(abs(latbath-latN_map)),1,'first');
idxWx = find(abs(lonbath-lonW_map)==nanmin(abs(lonbath-lonW_map)),1,'first');
idxEx = find(abs(lonbath-lonE_map)==nanmin(abs(lonbath-lonE_map)),1,'first');

% setting the bounds of the bathymetry data to the area that we're
% interested in
lonbath = lonbath(idxWx:idxEx);
latbath = latbath(idxSy:idxNy);
nlons = length(lonbath);
nlats = length(latbath);
bath = -ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','z',[idxWx idxSy],[nlons nlats]);

%% Separating V into Monthly Means

% January
jani = find(date(:,2) == 1);
jan_v_total = v(lon2D_i,lat2D_i,:,jani);
jan_v_mean = squeeze(mean(jan_v_total,4,'omitnan'));

% February
febi = find(date(:,2) == 2);
feb_v_total = v(lon2D_i,lat2D_i,:,febi);
feb_v_mean = squeeze(mean(feb_v_total,4,'omitnan'));

% March
mari = find(date(:,2) == 3);
mar_v_total = v(lon2D_i,lat2D_i,:,mari);
mar_v_mean = squeeze(mean(mar_v_total,4,'omitnan'));

% April
apri = find(date(:,2) == 4);
apr_v_total = v(lon2D_i,lat2D_i,:,apri);
apr_v_mean = squeeze(mean(apr_v_total,4,'omitnan'));

% May
mayi = find(date(:,2) == 5);
may_v_total = v(lon2D_i,lat2D_i,:,mayi);
may_v_mean = squeeze(mean(may_v_total,4,'omitnan'));

% June
juni = find(date(:,2) == 6);
jun_v_total = v(lon2D_i,lat2D_i,:,juni);
jun_v_mean = squeeze(mean(jun_v_total,4,'omitnan'));

% July
juli = find(date(:,2) == 7);
jul_v_total = v(lon2D_i,lat2D_i,:,juli);
jul_v_mean = squeeze(mean(jul_v_total,4,'omitnan'));

% August
augi = find(date(:,2) == 8);
aug_v_total = v(lon2D_i,lat2D_i,:,augi);
aug_v_mean = squeeze(mean(aug_v_total,4,'omitnan'));

% September
sepi = find(date(:,2) == 9);
sep_v_total = v(lon2D_i,lat2D_i,:,sepi);
sep_v_mean = squeeze(mean(sep_v_total,4,'omitnan'));

% October
octi = find(date(:,2) == 10);
oct_v_total = v(lon2D_i,lat2D_i,:,octi);
oct_v_mean = squeeze(mean(oct_v_total,4,'omitnan'));

% November
novi = find(date(:,2) == 11);
nov_v_total = v(lon2D_i,lat2D_i,:,novi);
nov_v_mean = squeeze(mean(nov_v_total,4,'omitnan'));

% December
deci = find(date(:,2) == 12);
dec_v_total = v(lon2D_i,lat2D_i,:,deci);
dec_v_mean = squeeze(mean(dec_v_total,4,'omitnan'));

%% Separating U into monthly means
% January
jani = find(date(:,2) == 1);
jan_u_total = u(lon2D_i,lat2D_i,:,jani);
jan_u_mean = squeeze(mean(jan_u_total,4,'omitnan'));

% February
febi = find(date(:,2) == 2);
feb_u_total = u(lon2D_i,lat2D_i,:,febi);
feb_u_mean = squeeze(mean(feb_u_total,4,'omitnan'));

% March
mari = find(date(:,2) == 3);
mar_u_total = u(lon2D_i,lat2D_i,:,mari);
mar_u_mean = squeeze(mean(mar_u_total,4,'omitnan'));

% April
apri = find(date(:,2) == 4);
apr_u_total = u(lon2D_i,lat2D_i,:,apri);
apr_u_mean = squeeze(mean(apr_u_total,4,'omitnan'));

% May
mayi = find(date(:,2) == 5);
may_u_total = u(lon2D_i,lat2D_i,:,mayi);
may_u_mean = squeeze(mean(may_u_total,4,'omitnan'));

% June
juni = find(date(:,2) == 6);
jun_u_total = u(lon2D_i,lat2D_i,:,juni);
jun_u_mean = squeeze(mean(jun_u_total,4,'omitnan'));

% July
juli = find(date(:,2) == 7);
jul_u_total = u(lon2D_i,lat2D_i,:,juli);
jul_u_mean = squeeze(mean(jul_u_total,4,'omitnan'));

% August
augi = find(date(:,2) == 8);
aug_u_total = u(lon2D_i,lat2D_i,:,augi);
aug_u_mean = squeeze(mean(aug_u_total,4,'omitnan'));

% September
sepi = find(date(:,2) == 9);
sep_u_total = u(lon2D_i,lat2D_i,:,sepi);
sep_u_mean = squeeze(mean(sep_u_total,4,'omitnan'));

% October
octi = find(date(:,2) == 10);
oct_u_total = u(lon2D_i,lat2D_i,:,octi);
oct_u_mean = squeeze(mean(oct_u_total,4,'omitnan'));

% November
novi = find(date(:,2) == 11);
nov_u_total = u(lon2D_i,lat2D_i,:,novi);
nov_u_mean = squeeze(mean(nov_u_total,4,'omitnan'));

% December
deci = find(date(:,2) == 12);
dec_u_total = u(lon2D_i,lat2D_i,:,deci);
dec_u_mean = squeeze(mean(dec_u_total,4,'omitnan'));

%% Separating Practical Salinity Variable into monthly means
% January
jani = find(date(:,2) == 1);
jan_s_total = psal(lon2D_i,lat2D_i,:,jani);
jan_s_mean = squeeze(mean(jan_s_total,4,'omitnan'));

% February
febi = find(date(:,2) == 1);
feb_s_total = psal(lon2D_i,lat2D_i,:,febi);
feb_s_mean = squeeze(mean(feb_s_total,4,'omitnan'));

% March
mari = find(date(:,2) == 3);
mar_s_total = psal(lon2D_i,lat2D_i,:,mari);
mar_s_mean = squeeze(mean(mar_s_total,4,'omitnan'));

% April
apri = find(date(:,2) == 4);
apr_s_total = psal(lon2D_i,lat2D_i,:,apri);
apr_s_mean = squeeze(mean(apr_s_total,4,'omitnan'));

% May
mayi = find(date(:,2) == 5);
may_s_total = psal(lon2D_i,lat2D_i,:,mayi);
may_s_mean = squeeze(mean(may_s_total,4,'omitnan'));

% June
juni = find(date(:,2) == 6);
jun_s_total = psal(lon2D_i,lat2D_i,:,juni);
jun_s_mean = squeeze(mean(jun_s_total,4,'omitnan'));

% July
juli = find(date(:,2) == 7);
jul_s_total = psal(lon2D_i,lat2D_i,:,juli);
jul_s_mean = squeeze(mean(jul_s_total,4,'omitnan'));

% August
augi = find(date(:,2) == 8);
aug_s_total = psal(lon2D_i,lat2D_i,:,augi);
aug_s_mean = squeeze(mean(aug_s_total,4,'omitnan'));

% September
sepi = find(date(:,2) == 9);
sep_s_total = psal(lon2D_i,lat2D_i,:,sepi);
sep_s_mean = squeeze(mean(sep_s_total,4,'omitnan'));

% October
octi = find(date(:,2) == 10);
oct_s_total = psal(lon2D_i,lat2D_i,:,octi);
oct_s_mean = squeeze(mean(oct_s_total,4,'omitnan'));

% November
novi = find(date(:,2) == 11);
nov_s_total = psal(lon2D_i,lat2D_i,:,novi);
nov_s_mean = squeeze(mean(nov_s_total,4,'omitnan'));

% December
deci = find(date(:,2) == 12);
dec_s_total = psal(lon2D_i,lat2D_i,:,deci);
dec_s_mean = squeeze(mean(dec_s_total,4,'omitnan'));

%% Converting SP to SA
z_a = repmat(z,1,143,479); % use repmat to create a 3x10x10 copy
z_3D = permute(z_a,[3 2 1]);

for j = 1:36
p(:,:,j) = gsw_p_from_z(-z_3D(:,:,j),lat2D);
end

% January
for j = 1:36
    [jan_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(jan_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% February
for j = 1:36
    [feb_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(feb_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% March
for j = 1:36
    [mar_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(mar_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% April
for j = 1:36
    [apr_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(apr_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% May
for j = 1:36
    [may_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(may_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% June
for j = 1:36
    [jun_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(jun_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% July
for j = 1:36
    [jul_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(jul_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% August
for j = 1:36
    [aug_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(aug_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% September
for j = 1:36
    [sep_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(sep_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% October
for j = 1:36
    [oct_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(oct_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% November
for j = 1:36
    [nov_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(nov_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

% December
for j = 1:36
    [dec_sa_mean(:,:,j), in_ocean] = gsw_SA_from_SP(dec_s_mean(:,:,j),p(:,:,j),lon2D,lat2D);
end

%% Flow and Salinity top 600m

sa_600 = cat(3,mean(jan_sa_mean(:,:,1:33),3,'omitnan'),mean(feb_sa_mean(:,:,1:33),3,'omitnan'),mean(mar_sa_mean(:,:,1:33),3,'omitnan'),...
    mean(apr_sa_mean(:,:,1:33),3,'omitnan'),mean(may_sa_mean(:,:,1:33),3,'omitnan'),mean(jun_sa_mean(:,:,1:33),3,'omitnan'),...
    mean(jul_sa_mean(:,:,1:33),3,'omitnan'),mean(aug_sa_mean(:,:,1:33),3,'omitnan'),mean(sep_sa_mean(:,:,1:33),3,'omitnan'),...
    mean(oct_sa_mean(:,:,1:33),3,'omitnan'),mean(nov_sa_mean(:,:,1:33),3,'omitnan'),mean(dec_sa_mean(:,:,1:33),3,'omitnan'));
u_600 = cat(3,mean(jan_u_mean(:,:,1:33),3,'omitnan'),mean(feb_u_mean(:,:,1:33),3,'omitnan'),mean(mar_u_mean(:,:,1:33),3,'omitnan'),...
    mean(apr_u_mean(:,:,1:33),3,'omitnan'),mean(may_u_mean(:,:,1:33),3,'omitnan'),mean(jun_u_mean(:,:,1:33),3,'omitnan'),...
    mean(jul_u_mean(:,:,1:33),3,'omitnan'),mean(aug_u_mean(:,:,1:33),3,'omitnan'),mean(sep_u_mean(:,:,1:33),3,'omitnan'),...
    mean(oct_u_mean(:,:,1:33),3,'omitnan'),mean(nov_u_mean(:,:,1:33),3,'omitnan'),mean(dec_u_mean(:,:,1:33),3,'omitnan'));
v_600 = cat(3,mean(jan_v_mean(:,:,1:33),3,'omitnan'),mean(feb_v_mean(:,:,1:33),3,'omitnan'),mean(mar_v_mean(:,:,1:33),3,'omitnan'),...
    mean(apr_v_mean(:,:,1:33),3,'omitnan'),mean(may_v_mean(:,:,1:33),3,'omitnan'),mean(jun_v_mean(:,:,1:33),3,'omitnan'),...
    mean(jul_v_mean(:,:,1:33),3,'omitnan'),mean(aug_v_mean(:,:,1:33),3,'omitnan'),mean(sep_v_mean(:,:,1:33),3,'omitnan'),...
    mean(oct_v_mean(:,:,1:33),3,'omitnan'),mean(nov_v_mean(:,:,1:33),3,'omitnan'),mean(dec_v_mean(:,:,1:33),3,'omitnan'));

[LON,LAT] = meshgrid(lon2D,lat2D);

%% Creating an animation: Reusing code for all depth levelsa    

% cd('/Users/saravianco/Documents/Research/MATLAB/GLORYS/GLORYS quiver/');
% writerObj = VideoWriter('vs_300_quiver_GLORYS.avi','Uncompressed AVI');
% writerObj.FrameRate = 1;
% open(writerObj);

u_ref = NaN(size(LON));
v_ref = NaN(size(LON));
u_ref(41,50) = 0.25;
v_ref(41,50) = 0;

% 
% for ii = 1:12
% fig1 = figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
% clf
% hold on
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
% m_pcolor(LON,LAT,sa_300(:,:,ii)');
% cmocean('haline')
% x = colorbar;
% m_quiver(LON(1:3:end,1:5:end),LAT(1:3:end,1:5:end),u_300(1:5:end,1:3:end,ii)',...
%     v_300(1:5:end,1:3:end,ii)',2.5,'k','LineWidth',0.5,'Clipping','on');
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r'); % this plots the FSAOO mooring locations
% 
% 
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','w');
% ylabel(x,'SA [g/kg]')
% set(x,'fontsize',10,'fontname','helvetica','position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',29:1:35);
% clim([29.5 35.5]);
% % m_quiver(LON,LAT,u_ref,v_ref,2.5,'k','LineWidth',0.7);
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
%     'ytick',60:5:80,'xtick',-35:15:5);
% m_quiver(LON,LAT,u_ref,v_ref,2.5,'k','LineWidth',0.7);
% 
% dim = [0.38 0.35 .045 .06];
% str = '0.25 m/s';
% annotation('textbox',dim,'String',str,'FontSize',13,'FitBoxToText','off','Rotation',-14.5);
% 
% formatSpec = '%s NEGS Velocity and Salinity in the top 300m (GLORYS)';
% months = ["January","February","March","April","May","June","July","August",...
%     "September","October","November","December"];
% title(sprintf(formatSpec,months(ii)),'FontSize',16)
% 
% 
% F = getframe(fig1);
% writeVideo(writerObj, F);
% 
% end
% close(writerObj);


%% Time mean Quiver Plot for top 600 m
v_600_avg = mean(v_600,3,'omitnan');
u_600_avg = mean(u_600,3,'omitnan');
sa_600_avg = mean(sa_600,3,'omitnan');
speed_600_avg = sqrt(u_600_avg.^2+v_600_avg.^2);

%% fixing reference quiver
u_600_avg(41,50) = 0.25;
v_600_avg(41,50) = 0;
% 
% Salinity and Velocity
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
% hold on
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
% m_pcolor(LON,LAT,sa_300_avg');
% cmocean('haline')
% x = colorbar;
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r'); % this plots the FSAOO mooring locations
% 
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
%     'ytick',60:5:80,'xtick',-35:15:5);
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','w');
% ylabel(x,'SA [g/kg]')
% set(x,'fontsize',10,'fontname','helvetica','position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',29:1:35);
% clim([29.5 35.5]);
% % m_quiver(LON,LAT,u_ref,v_ref,2.5,'k','LineWidth',0.7);
% m_quiver(LON(1:3:end,1:5:end),LAT(1:3:end,1:5:end),u_300(1:5:end,1:3:end,ii)',...
%   v_300(1:5:end,1:3:end,ii)',2.5,'k','LineWidth',0.5,'Clipping','on');
% dim = [0.36 0.35 .045 .06];
% str = '0.25 m/s';
% annotation('textbox',dim,'String',str,'FontSize',13,'FitBoxToText','off','Rotation',-14.5);
% 
% title('Time Mean Velocity and Salinity: GLORYS 2003-2019','FontSize',16)
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/GLORYS quiver/GLORYS_tmean_v_sal.jpg']);


%% Determining the angle of the bathymetry (Adding this here becasue I use the 
% 600m contour in the following plot)

%% Determining the angle of the bathymetry
% C gives me the x and y coordinates of the 600m contour line
bath_600 = bath;
% manually reworking the bathymetry matrix to yield only the outer shelf
% 600m contour
bath_600(1:1057,1:422) = 0;
bath_600(1501:1801,1:62) = 1000;
bath_600(241:313,680:end) = 1000;
[C,h] = contour(lonbath,latbath,bath_600',[600 600]);
x_600 = C(1,2:end);
y_600 = C(2,2:end);

% smoothing the x and y data (still in x and y coordinates
% to lat and lon...)
lat_600_sm = smooth(y_600,100)';
lon_600_sm = smooth(x_600,100)';


% this figure tests the original contour against the smoothed one
figure
hold on
plot(x_600,y_600,'-b')
plot(lon_600_sm,lat_600_sm,'-r');

%% Velocity and Speed
addpath('/Users/saravianco/Documents/Research/MATLAB/Share w Sara Vianco/')

figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(LON,LAT,speed_600_avg');
cmocean('speed')
x = colorbar;
m_quiver(LON(1:2:end,1:10:end),LAT(1:2:end,1:10:end),u_600(1:10:end,1:2:end,12)',...
  v_600(1:10:end,1:2:end,12)',2.5,'k','LineWidth',0.5,'Clipping','on','AutoScale','on');
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','w');
% b6 = m_plot(lon_600_sm,lat_600_sm,'-m','LineWidth',1);
m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');

% this part is plotting the k mean clusters
% m_plot([lon_600_poly(idx2_1_ii)-2 lon_600_poly(idx2_1_ii)+2],[lat_600_poly(idx2_1_ii) lat_600_poly(idx2_1_ii)],'linewi',2,'color','g');
% m_plot([lon_600_poly(idx2_1_if)-2 lon_600_poly(idx2_1_if)+2],[lat_600_poly(idx2_1_if) lat_600_poly(idx2_1_if)],'linewi',2,'color','g');
% m_plot([lon_600_poly(idx2_1_fi)-2 lon_600_poly(idx2_1_fi)+2],[lat_600_poly(idx2_1_fi) lat_600_poly(idx2_1_fi)],'linewi',2,'color','g');
% m_plot([lon_600_poly(idx2_1_ff)-2 lon_600_poly(idx2_1_ff)+2],[lat_600_poly(idx2_1_ff) lat_600_poly(idx2_1_ff)],'linewi',2,'color','g');
% m_plot([lon_600_poly(idx2_2_i)-2 lon_600_poly(idx2_2_i)+2],[lat_600_poly(idx2_2_i) lat_600_poly(idx2_2_i)],'linewi',2,'color','b');
% m_plot([lon_600_poly(idx2_2_f)-2 lon_600_poly(idx2_2_f)+2],[lat_600_poly(idx2_2_f) lat_600_poly(idx2_2_f)],'linewi',2,'color','b');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);

set(x,'fontsize',10,'fontname','helvetica','position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',0:0.05:0.2);
clim([0 0.2]);
m_quiver(LON,LAT,u_ref,v_ref,1,'k','LineWidth',0.7,'MaxHeadSize',2);

dim = [0.36 0.35 .045 .06];
% dim = [0.38 0.35 .045 .06];
str = '0.25 m/s';
annotation('textbox',dim,'String',str,'FontSize',13,'FitBoxToText','off','Rotation',-14.5);

ylabel(x,'Speed [m/s]')
title('Average December Velocity and Speed (GLORYS)','FontSize',16)
% legend(b6,'600 m isobath')
% this specific file accounts for adding in the smoothed bathymetry (600m
% contour)
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/velcity_speed_monthly/v_spd_12.jpg']);

%% Identifying Polygon coordinates: Nordostrundingen latitude
% [x,y] = ginput(1);
load('nord_xy.mat');
[lon_nord,lat_nord] = m_xy2ll(x,y);

% attempting to determine the intersection between the northern boundary
% and the 600m isobath
int_0 = lat_600_sm-lat_nord;

nord_600 = min(abs(int_0));
[nord_600_r,nord_600_c] = find(abs(int_0)==nord_600);

polyNE_lon = x_600(nord_600_r,nord_600_c);
polyNE_lat = y_600(nord_600_r,nord_600_c);

%% ID-ing cordinates: 74N latitude
% attempting to determine the intersection between the northern boundary
% and the 600m isobath
int_0 = lat_600_sm-74;

seven4_600 = min(abs(int_0));
[seven4_600_r,seven4_600_c] = find(abs(int_0)==seven4_600);

polySE_lon = x_600(seven4_600_r,seven4_600_c);
polySE_lat = y_600(seven4_600_r,seven4_600_c);

%%
% [x,y] = ginput(1); determining the coordinate of where the 74N lat
% encounters the coast

load("74_xy.mat")
[lon_74,lat_74] = m_xy2ll(x,y);

%% Curtailed 600m isobath

% cutting down the 600m isobath contour to only contain the section from
% 74N to Nord lat
lat_600_poly = lat_600_sm(1,nord_600_c:seven4_600_c);
lon_600_poly = lon_600_sm(1,nord_600_c:seven4_600_c);
% cutting down the slope angle variable to match our curtailed contour
% slope_ang_poly = slope_ang(1,nord_600_c:seven4_600_c);


%% a testing figure for the polygon limits

figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
cmocean('deep')

b6 = m_plot(lon_600_poly,lat_600_poly,'-r','LineWidth',2);
m_plot([lon_nord polyNE_lon],[lat_nord lat_nord],'-r','LineWidth',2)
m_plot([lon_74 polySE_lon],[polySE_lat polySE_lat],'-r','LineWidth',2)
m_plot(polyNE_lon,polyNE_lat,'*w')
m_plot(polySE_lon,polySE_lat,'*w')

% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
title('Shelf Polygon','FontSize',16)

% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/GLORYS quiver/shelf_polygon.jpg']);


%% Finding the coordinates of the quivers along the 600m isobath

% for the original coordinate system (lon2D and lat2D) I need to identify
% the coordinates of the 600m isobath, roughly

% first to constrain the original latitudes by the polygon limits
idxlatN = find(abs(lat2D-polyNE_lat)==nanmin(abs(lat2D-polyNE_lat)),1,'first');
idxlatS = find(abs(lat2D-polySE_lat)==nanmin(abs(lat2D-polySE_lat)),1,'first');
idxlonW = find(abs(lon2D-polySE_lon)==nanmin(abs(lon2D-polySE_lon)),1,'first');
idxlonE = find(abs(lon2D-polyNE_lon)==nanmin(abs(lon2D-polyNE_lon)),1,'first');

for i = 1:12
v_600m_res(i,:) = interp2(lat2D,lon2D,v_600(:,:,i),lat_600_poly,lon_600_poly,'nearest');
u_600m_res(i,:) = interp2(lat2D,lon2D,u_600(:,:,i),lat_600_poly,lon_600_poly,'nearest');
end
% it worked!
v_600m_res = v_600m_res';
u_600m_res = u_600m_res';

% finding the high-res shelf angle
clear x
clear y

[x,y] = m_ll2xy(lon_600_poly, lat_600_poly);

% calculating a simple slope 
for m = 1:(length(lon_600_poly)-1)
dy(m) = (y(m+1)-y(m));
dx(m) = (x(m+1)-x(m));
end

ang_rough = atan2d(dy,dx);
vec_ang = -ang_rough;
vec_ang(1,1606) = vec_ang(1,1605);

v_600_cont = reshape(v_600m_res,[1,(1606*12)]);
u_600_cont = reshape(u_600m_res,[1,(1606*12)]);

U = cat(1,u_600_cont,v_600_cont);

R = NaN(2,2,(length(u_600_cont)));

vec_ang_yr = repmat(vec_ang,1,12);

for jj = 1:length(vec_ang_yr)
R(:,:,jj)=[cosd(vec_ang_yr(jj)), -sind(vec_ang_yr(jj)) ; sind(vec_ang_yr(jj)), cosd(vec_ang_yr(jj))];
U_r(:,jj) = R(:,:,jj)*U(:,jj);
end

u_600r = U_r(1,:);
v_600r = U_r(2,:);

u_600r = reshape(u_600r,[1606,12]);
v_600r = reshape(v_600r,[1606,12]);

u_600r = u_600r(2:end,:);
v_600r = v_600r(2:end,:);
% %% testing the coordinates in GLORYS coordinate system
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
% hold on
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
% m_pcolor(lonbath,latbath,bath')
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
% cmocean('deep')
% m_plot(lon_600,lat_600,'-r','LineWidth',2)
% % b6 = m_plot(lon_600_poly,lat_600_poly,'-r','LineWidth',2);
% % m_plot([lon_nord polyNE_lon],[lat_nord lat_nord],'-r','LineWidth',2)
% % m_plot([lon_74 polySE_lon],[polySE_lat polySE_lat],'-r','LineWidth',2)
% % m_plot(polyNE_lon,polyNE_lat,'*w')
% % m_plot(polySE_lon,polySE_lat,'*w')
% 
% % m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
%     'ytick',60:5:80,'xtick',-35:15:5);

%% now to test in a plot
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
cmocean('deep')
m_plot(lon_600_poly,lat_600_poly,'-r','LineWidth',2)   
% m_quiver(ll_600(1:5:end,1),ll_600(1:5:end,2),u_600_con(1:5:end),...
%   v_600_con(1:5:end),1,'k','LineWidth',0.5,'Clipping','on');
% m_quiver(lon_600_poly(1:5:end),lat_600_poly(1:5:end),u_600_con(1:5:end),v_600_con(1:5:end),1,'k','LineWidth',0.5,'Clipping','on');
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
m_quiver(lon_600_poly(1:100:end), lat_600_poly(1:100:end), dx(1:100:end), dy(1:100:end), '-r')                          % Draw lquivero Arrows
hold off
axis('equal')

angstr = sprintfc('\\angle %.0f \\circ',ang_rough);
text(x(1:100:end), y(1:100:end), angstr(1:100:end), 'HorizontalAlignment','left')
% 
% m_quiver(LON,LAT,u_ref,v_ref,2.5,'k','LineWidth',0.7);
title('Tangent Angles along 600m isobath')
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/shelf_angles.jpg']);

%%
% speed_600c = sqrt(u_600_con.^2+v_600_con.^2)';
% 
% for jj = 1:length(idx_lon600)
% speed_600_check(jj) = speed_600_avg(idx_lon600(jj),idx_lat600(jj));
% end
% 
% test = speed_600c-speed_600_check;
% % these should come out to zero, and they do! so maybe the difference is in
% % the scaling

%% to diagnose onshore and offshore transport, -v is onshore and +v is offshore

lon_km = NaN(length(lon_600_poly)-1,1);

for jj = 1:(length(lon_600_poly)-1)
lon_km(jj,1) = haversine([lat_600_poly(jj),lon_600_poly(jj)],[lat_600_poly(jj+1),lon_600_poly(jj+1)]); 
% haversign calculates the geometric distance in km between two coordinates
end

lon_m = lon_km*1000;

% to calculate onshelf transport, we're going to find where v is negative and
% multiply that against the size of our grid, as well as by 600m, since our
% velocities are averaged over the top 600m. That way our final units will
% be in Sv
vol_tx_600G = NaN(size(v_600r));

onshelf_idx = find(v_600r<0);
v_onshelf = zeros(size(v_600r));
v_onshelf(onshelf_idx) = v_600r(onshelf_idx);

onshelf_vx = 600*v_onshelf.*lon_m;
vol_tx_600G(onshelf_idx) = onshelf_vx((onshelf_idx));

for i = 1:12
onshelf_vx_tot(i) = (nansum(onshelf_vx(:,i),1))/10^6;
end
% this is the total onshelf volumetric transport (average) for each month,
% summed along the 600m polygon

% same thing but off the shelf
offshelf_idx = find(v_600r>0);
v_offshelf = zeros(size(v_600r));
v_offshelf(offshelf_idx) = v_600r(offshelf_idx);

offshelf_vx = 600*v_offshelf.*lon_m;
vol_tx_600G(offshelf_idx) = offshelf_vx((offshelf_idx));

for i = 1:12
offshelf_vx_tot(i) = (nansum(offshelf_vx(:,i),1))/10^6;
end

% converting to sverdrups
vol_tx_600G = vol_tx_600G/10^6;

% averaging across all years
vx_yr_avg = mean(vol_tx_600G,2);

% %% testing out hist3
% markerSize = 100;
% figure
% geoscatter(lat_600_poly, lon_600_poly, markerSize, vol_tx_600,'Filled');
% colorbar

%% displaying transport on/off the shelf
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
% cmocean('deep')
% m_plot(lon_600,lat_600,'-g','LineWidth',1)
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
m_scatter(lon_600_poly(1:end-1),lat_600_poly(1:end-1),100,(1000*vx_yr_avg),'Filled');
clim([-30 5])
cmocean('balance','pivot',0)
h = colorbar;
% m_quiver(lon_600_poly,lat_600_poly,u_600_con,v_600_con,1,'k','LineWidth',0.5,'Clipping','on');
title('Average Cross-Shelf Volumetric Transport in the top 600m (GLORYS)')
ylabel(h,'Volumetric Transport [mSv]')

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/Vol_tx_600.jpg']);

%%
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
% hold on
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
% %m_pcolor(lonbath,latbath,bath')
% % m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
% m_plot(lon_600,lat_600,'-k','LineWidth',1)
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
%     'ytick',60:5:80,'xtick',-35:15:5);
% % m_plot(ll_600_poly(:,1),ll_600_poly(:,2),'-g','LineWidth',1)
% m_quiver(lon_600_poly(1:2:end), lat_600_poly(1:2:end), u_600_con(1:2:end),...
%     v_600_con(1:2:end),'b','LineWidth',0.5,'Clipping','on');                       
% % m_quiver(lon_600_poly(1:2:end-1,1), lat_600_poly(1:2:end-1,1), dx(1,1:2:end), dy(1,1:2:end),...
% %     '-r','ShowArrowHead','off')                         
% hold off
% % axis('equal')
% 
% angstr = sprintfc('\\angle %.0f \\circ',ang_rough);
% text(x(1:2:end-1), y(1:2:end-1), angstr(1:2:end), 'HorizontalAlignment','left')
% title('600m Isobath Angles')
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/GLORYS quiver/600m_angles.jpg']);



%% displaying transport on/off the shelf 
% (same code as above, but running for different months)
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
% cmocean('deep')
% m_plot(lon_600,lat_600,'-g','LineWidth',1)
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
m_scatter(lon_600_poly(1:end-1),lat_600_poly(1:end-1),100,(1000*vol_tx_600G(:,12)),'Filled');
h = colorbar;
clim([-30 5])
colormap(cmocean('balance','pivot',0))
title('Average December Cross-Shelf Volumetric Transport (GLORYS)')
ylabel(h,'Volumetric Transport [mSv]')

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/vol tx monthly/vx_600_12.jpg']);

%% visualizing the volumetric transport by plotting against distance

months = {'Jan','Feb','Mar','Apr','May','Jul','Jul','Aug','Sep','Oct','Nov','Dec'};
% establishing a variable that's distance from the northern tip of our
% polygon

dist_nord = cumsum(lon_m(2:end));

figure
hold on
imagesc([1:12],dist_nord/1000,1000*vol_tx_600G(2:end,:));
set(gca,'Ydir','reverse');
x = colorbar;
clim([-40 40])
colormap(cmocean('balance','pivot',0))

xlim([0.5 12.5])
xticks(1:12);
xticklabels(months)
ylabel(x,'Volumetric Transport [mSv]','FontSize',10)
title('Volumetric Transport along the NEGS (GLORYS)')
subtitle('Top 600m')
set(gcf,'color','w');
ylabel('Distance from Nordostrundingen [km]','FontSize',10)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/vol tx monthly/monthly_vx_03_19.jpg']);

%% using a line plot to image the temporal (not super helpful)

c = cmocean('thermal',1605);
shelf_mean_vx = mean(vol_tx_600G,1);
sm_vx_smooth = movmean(shelf_mean_vx,6);

p = polyfit(1:12,shelf_mean_vx,4);
y1 = polyval(p,1:12);

y1_sm = smooth(y1,2);

figure
hold on
% for ii = 1:100:1606
% plot([1:12],1000*vol_tx_600(ii,:),'linewi',1,'Color',c(ii,:));
% end
plot(1:12, 1000*shelf_mean_vx,'k','linewi',2)
% plot(1:12, 1000*sm_vx_smooth,'b','linewi',2);
plot(1:12,1000*y1_sm,'r','linewi',2)

xlim([1 12])
xticks(1:12);
xticklabels(months)
legend('Mean VT','Polynomial Fit')
% ylim([-0.18 0])
ylabel('VT [mSv]')

title('Seasonal Cycle of Cross-shelf Volumetric Transport')
set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/GLORYS quiver/vol tx monthly/sea_cyc_line.jpg']);

% Across the shelf, average monthly volumetric transport is all on-shelf.
% Meaning that there must be some dynamic on the shelf that's exporting
% this water that's coming on

%% Looking at the Volumetric Transport for different seasons

% to plot this value, we can look at different seasons, and further down,
% we ca look at indvidual years if we wanted to 
vol_tx_600G = vol_tx_600G(2:end,:);

VT_winter = 1000*mean(vol_tx_600G(:,[12,1,2]),2);
VT_spring = 1000*mean(vol_tx_600G(:,[3,4,5]),2);
VT_summer = 1000*mean(vol_tx_600G(:,[6,7,8]),2);
VT_fall = 1000*mean(vol_tx_600G(:,[9,10,11]),2);

figure
hold on
plot(dist_nord/1000,VT_winter,'b')
plot(dist_nord/1000,VT_spring,'g')
plot(dist_nord/1000,VT_summer,'m')
plot(dist_nord/1000,VT_fall,'r')
yline(0,'--k','linewidth',1)

ylim([-50 20])
legend('Winter (DJF)','Spring (MAM)','Summer (JJA)','Autumn (SON)','location','southeast')
xlabel('Distance along 600m Isobath [km]')
ylabel('VT [mSv]')
title('Seasonal Volumetric Transport (GLORYS)')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/vol tx monthly/sea_vt_line.jpg']);

%% looking at the seasonal cycle across all years and months

u_all_z = squeeze(mean(u(:,:,1:33,:),3,'omitnan'));
v_all_z = squeeze(mean(v(:,:,1:33,:),3,'omitnan'));

for i = 1:191
v_600_mo(:,i) = interp2(lat2D,lon2D,v_all_z(lon2D_i,lat2D_i,i),lat_600_poly(2:end),lon_600_poly(2:end),'nearest');
u_600_mo(:,i) = interp2(lat2D,lon2D,u_all_z(lon2D_i,lat2D_i,i),lat_600_poly(2:end),lon_600_poly(2:end),'nearest');
end

v_600_long = reshape(v_600_mo,[1,(1605*191)]);
u_600_long = reshape(u_600_mo,[1,(1605*191)]);

U_mo = cat(1,u_600_long,v_600_long);

R_mo = NaN(2,2,(length(u_600_long)));

vec_ang_mo = repmat(vec_ang(2:end),1,191);

for jj = 1:length(vec_ang_mo)
R_mo(:,:,jj)=[cosd(vec_ang_mo(jj)), -sind(vec_ang_mo(jj)) ; sind(vec_ang_mo(jj)), cosd(vec_ang_mo(jj))];
U_mo_r(:,jj) = R_mo(:,:,jj)*U_mo(:,jj);
end

u_600r_mo = U_mo_r(1,:);
v_600r_mo = U_mo_r(2,:);

u_600r_mo = reshape(u_600r_mo,[1605,191]);
v_600r_mo = reshape(v_600r_mo,[1605,191]);

%% determining on/off-shore transport
vx_600_GLORYS = NaN(size(v_600r_mo));

onshelf_idx = find(v_600r_mo<0);
v_onshelf = zeros(size(v_600r_mo));
v_onshelf(onshelf_idx) = v_600r_mo(onshelf_idx);

onshelf_vx = 600*v_onshelf.*lon_m;
vx_600_GLORYS(onshelf_idx) = onshelf_vx((onshelf_idx));

for i = 1:191
onshelf_vx_tot_mo(i) = (nansum(onshelf_vx(:,i),1))/10^6;
end
% this is the total onshelf volumetric transport (average) for each month,
% summed along the 600m polygon

% same thing but off the shelf
offshelf_idx = find(v_600r_mo>0);
v_offshelf = zeros(size(v_600r_mo));
v_offshelf(offshelf_idx) = v_600r_mo(offshelf_idx);

offshelf_vx = 600*v_offshelf.*lon_m;
vx_600_GLORYS(offshelf_idx) = offshelf_vx((offshelf_idx));

for i = 1:191
offshelf_vx_tot_mo(i) = (nansum(offshelf_vx(:,i),1))/10^6;
end

% converting to sverdrups
vx_600_GLORYS = vx_600_GLORYS/10^6;

% averaging across all months/years
vx_mo_avg = mean(vx_600_GLORYS,2);

%% Imaging this Yearly Volumetric Transport

% calculating different means
% Mean volumetric transport across our entire section for each month of
% each year

vx_yr_SectionAvg = squeeze(mean(vx_600_GLORYS,2,'omitnan'));
vx_year_tavg = nanmean(vx_yr_SectionAvg,2);

c = cmocean('thermal',17);

% figure
% hold on
% for K = 1:17
%     plot([1:12],vx_yr_SectionAvg(:,K),'Color',c(K,:))
% end
%     plot([1:12],vx_year_tavg,'-k','LineWidth',2)
% legend('2003','2004','2005','2006','2007','2008','2009','2010','2011',...
%     '2012','2013','2014','2015','2016','2017','2018','2019','location','eastoutside')
% 
% xlim([1 12])
% xticks(1:12);
% xticklabels(months)
% ylim([-0.4 0.15])
% ylabel('VT [Sv]')
% title('Cross-shelf Volumetric Transport for 2003-2019 (GLORYS)')
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/GLORYS quiver/vol tx monthly/year_seas.jpg']);

%% Volumetric Transport Anomalies

vx_yr_anom = vx_yr_SectionAvg-vx_year_tavg;

% figure
% hold on
% for K = 1:17
%     plot([1:12],vx_yr_anom(:,K),'Color',c(K,:))
% end
% legend('2003','2004','2005','2006','2007','2008','2009','2010','2011',...
%     '2012','2013','2014','2015','2016','2017','2018','2019','location','eastoutside')
% 
% xlim([1 12])
% xticks(1:12);
% xticklabels(months)
% ylim([-0.2 0.2])
% ylabel('VT [Sv]')
% title('Cross-shelf Volumetric Transport Anomalies for 2003-2019 (GLORYS)')
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/GLORYS quiver/vol tx monthly/year_seas_anom.jpg']);

%% Creating a colormap across all months/years
% need to adjust the time and the volumetric transport arrays to both
% resemble 1606x204

% [r,c,l] = size(vol_tx_GLORYS);
% vx_all = reshape(vol_tx_GLORYS,[r,c*l]);

date_mat = cell2mat(date_num);
[r,c] = size(date_mat);
date_num_all = reshape(date_mat,[1,r*c]);
date_hyp = datetime(date_num_all, 'ConvertFrom', 'datenum', 'Format', 'MM-dd-yyyy');  

%%
vx_fit_GLORYS = NaN(1605,204);
vx_fit_GLORYS(:,9:199) = vx_600_GLORYS;

figure
hold on
imagesc(date_num_all,dist_nord/1000,1000*vx_fit_GLORYS);
set(gca,'Ydir','reverse');
dateFormat = 11;
datetick('x',dateFormat)
xlim([date_num_all(1,1) date_num_all(end)])
x = colorbar;
clim([-50 50])
colormap(cmocean('balance','pivot',0))
ylabel('Distance from Nordostrundingen [km]','FontSize',10)
ylabel(x,'VT [mSv]')
title('Cross-shelf Volumetric Transport for 2003-2019 (GLORYS)')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/vol tx monthly/vt_all.jpg']);

%% Volumetric Transport Anomalies (x shelf)

% first looking at the anomalies across all months/years based on the
% average at each location
vx_test = nanmean(vx_fit_GLORYS,2);
vx_anom_shelf = vx_fit_GLORYS - vx_test;

figure
hold on
imagesc(date_num_all,dist_nord/1000,1000*vx_anom_shelf);
set(gca,'Ydir','reverse');
dateFormat = 11;
datetick('x',dateFormat)
xlim([date_num_all(1,1) date_num_all(end)])
x = colorbar;
clim([-40 40])
colormap(cmocean('balance','pivot',0))
ylabel('Distance from Nordostrundingen [km]','FontSize',10)
ylabel(x,"VT' [mSv]")
title('Cross-shelf Volumetric Transport Anomalies for 2003-2019 (GLORYS)')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/vol tx monthly/vt_all_anom.jpg']);

%% Cluster Analysis: Creating a Dendrogram 
% the purpose of this is to try and diagnose how many clusters there are
% within the volumetric transport data

tree = linkage(vx_600_GLORYS,'average');

figure()
dendrogram(tree,25,'ColorThreshold','default')

%% K-means clustering

% 1 cluster
[idx1,C] = kmeans(vx_600_GLORYS,1);
[silh1,h] = silhouette(vx_600_GLORYS,idx1);

cluster1 = mean(silh1); % NaN

% 2 cluster
[idx2,C] = kmeans(vx_600_GLORYS,2);
[silh2,h] = silhouette(vx_600_GLORYS,idx2);

cluster2 = mean(silh2); 

% 3 cluster
[idx3,C] = kmeans(vx_600_GLORYS,3);
[silh3,h] = silhouette(vx_600_GLORYS,idx3);

cluster3 = mean(silh3); 

% 4 cluster
[idx4,C] = kmeans(vx_600_GLORYS,4);
[silh4,h] = silhouette(vx_600_GLORYS,idx4);

cluster4 = mean(silh4); 

% 5 cluster
[idx5,C] = kmeans(vx_600_GLORYS,5);
[silh5,h] = silhouette(vx_600_GLORYS,idx5);

cluster5 = mean(silh5); 

% 6 cluster
[idx6,C] = kmeans(vx_600_GLORYS,6);
[silh6,h] = silhouette(vx_600_GLORYS,idx6);

cluster6 = mean(silh6); 

% 7 cluster
[idx7,C] = kmeans(vx_600_GLORYS,7);
[silh7,h] = silhouette(vx_600_GLORYS,idx7);

cluster7 = mean(silh7); 

% 8 cluster
[idx8,C] = kmeans(vx_600_GLORYS,8);
[silh8,h] = silhouette(vx_600_GLORYS,idx8);

cluster8 = mean(silh8); 

% 9 cluster
[idx9,C] = kmeans(vx_600_GLORYS,9);
[silh9,h] = silhouette(vx_600_GLORYS,idx9);

cluster9 = mean(silh9); 

% 10 cluster
[idx10,C] = kmeans(vx_600_GLORYS,10);
[silh10,h] = silhouette(vx_600_GLORYS,idx10);

cluster10 = mean(silh10); 


%%
figure
scatter(idx4(2:end),dist_nord/1000,'k','Filled');
set(gca,'Ydir','reverse');

title('K-Means Clustering Results (4 clusters, GLORYS)')
xlabel('Cluster Number','FontSize',10)
ylabel('Distance from Nordostrundingen [km]','FontSize',10)

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/vol tx monthly/idx4.jpg']);


% this appears to explain a lot of the variation in the section in GLORYS.
% There's two clusters but 3 regimes that seem to be well mimicked on the
% volumetric transport heatmaps. However, the top and bottom are different
% areas dynamically, so there shouldn't be that much of a correlation...

% 2 clusters also has the highest silhouette value out of the first 8
% situations, meaning that on average, 2 clusters yields the shortest
% distance between the different individual points. AKA 2 clusters most
% accurately represents the grouping of points wihin this population of
% obserations. 

%% Plotting 100/200/300...km along shelfbreak (taken from RARE code to get var)

dist_nord_km = dist_nord/1000;
xq = [100:100:900];
dist_nord_km_highres = interp1(dist_nord_km,dist_nord_km,xq);

for dd = 1:length(xq)
dist_100s_idx(dd) = find(abs(dist_nord_km-xq(dd))==nanmin(abs(dist_nord_km-xq(dd))));
end

xq_string = string(xq);

%%
% idx2_2_i = 43;
% idx2_2_f = find(idx2 == 2,1,'last');
% idx2_1_ii = find(idx2 == 1,1,'first');
% idx2_1_if = 42;
% idx2_1_fi = idx2_2_f+1;
% idx2_1_ff = find(idx2 == 1,1,'last');
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
% hold on
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
% m_pcolor(lonbath,latbath,bath')
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
% cmocean('deep')
% 
% b6 = m_plot(lon_600_poly,lat_600_poly,'-r','LineWidth',2);
% m_plot([lon_nord polyNE_lon],[lat_nord lat_nord],'-r','LineWidth',2)
% m_plot([lon_74 polySE_lon],[polySE_lat polySE_lat],'-r','LineWidth',2)
% % m_plot(polyNE_lon,polyNE_lat,'*w')
% % m_plot(polySE_lon,polySE_lat,'*w')
% m_scatter(lon_600_poly(dist_100s_idx),lat_600_poly(dist_100s_idx),100,lat_600_poly(dist_100s_idx)','*w');
% m_plot([lon_600_poly(idx2_1_ii)-2 lon_600_poly(idx2_1_ii)+2],[lat_600_poly(idx2_1_ii) lat_600_poly(idx2_1_ii)],'linewi',2,'color','g');
% m_plot([lon_600_poly(idx2_1_if)-2 lon_600_poly(idx2_1_if)+2],[lat_600_poly(idx2_1_if) lat_600_poly(idx2_1_if)],'linewi',2,'color','g');
% m_plot([lon_600_poly(idx2_1_fi)-2 lon_600_poly(idx2_1_fi)+2],[lat_600_poly(idx2_1_fi) lat_600_poly(idx2_1_fi)],'linewi',2,'color','g');
% m_plot([lon_600_poly(idx2_1_ff)-2 lon_600_poly(idx2_1_ff)+2],[lat_600_poly(idx2_1_ff) lat_600_poly(idx2_1_ff)],'linewi',2,'color','g');
% m_plot([lon_600_poly(idx2_2_i)-2 lon_600_poly(idx2_2_i)+2],[lat_600_poly(idx2_2_i) lat_600_poly(idx2_2_i)],'linewi',2,'color','b');
% m_plot([lon_600_poly(idx2_2_f)-2 lon_600_poly(idx2_2_f)+2],[lat_600_poly(idx2_2_f) lat_600_poly(idx2_2_f)],'linewi',2,'color','b');
% 
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
%     'ytick',60:5:80,'xtick',-35:15:5);
% 
% title('Cluster Locations along 600m (GLORYS)','FontSize',16)
% 
% % imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% % 'GLORYS/GLORYS quiver/vol tx monthly/geo_kmean_sections.jpg']);
% 
% %% sections on heatmap
% figure
% hold on
% imagesc(date_num_all,dist_nord_km,vx_all);
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_ii) dist_nord_km(idx2_1_ii)],'linewi',2,'color','g');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_if) dist_nord_km(idx2_1_if)],'linewi',2,'color','g');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_fi) dist_nord_km(idx2_1_fi)],'linewi',2,'color','g');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_ff) dist_nord_km(idx2_1_ff)],'linewi',2,'color','g');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_2_i) dist_nord_km(idx2_2_i)],'linewi',2,'color','b');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_2_f) dist_nord_km(idx2_2_f)],'linewi',2,'color','b');
% 
% set(gca,'Ydir','reverse');
% dateFormat = 11;
% datetick('x',dateFormat)
% xlim([date_num_all(1,1) date_num_all(end)])
% x = colorbar;
% clim([-0.8 0.8])
% colormap(cmocean('balance','pivot',0))
% ylabel('Distance from Nordostrundingen [km]','FontSize',10)
% ylabel(x,'VT [Sv]')
% title('Cross-shelf Volumetric Transport for 2003-2019 (GLORYS)')
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/GLORYS quiver/vol tx monthly/heatmap_kmean_sections.jpg']);
% 
% %% sections on anomaly heatmap
% 
% figure
% hold on
% imagesc(date_num_all,dist_nord/1000,vx_anom_shelf);
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_ii) dist_nord_km(idx2_1_ii)],'linewi',2,'color','g');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_if) dist_nord_km(idx2_1_if)],'linewi',2,'color','g');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_fi) dist_nord_km(idx2_1_fi)],'linewi',2,'color','g');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_ff) dist_nord_km(idx2_1_ff)],'linewi',2,'color','g');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_2_i) dist_nord_km(idx2_2_i)],'linewi',2,'color','b');
% plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_2_f) dist_nord_km(idx2_2_f)],'linewi',2,'color','b');
% 
% set(gca,'Ydir','reverse');
% dateFormat = 11;
% datetick('x',dateFormat)
% xlim([date_num_all(1,1) date_num_all(end)])
% x = colorbar;
% clim([-0.8 0.8])
% colormap(cmocean('balance','pivot',0))
% ylabel('Distance from Nordostrundingen [km]','FontSize',10)
% ylabel(x,"VT' [Sv]")
% title('Cross-shelf Volumetric Transport Anomalies for 2003-2019 (GLORYS)')
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/GLORYS quiver/vol tx monthly/heatmapA_kmean_sections.jpg']);
% 
% %%
% % figure
% % histogram(vx_g2) % (vx_g1) % (vx_g3)
% % ^ each of the three sections were tested for normality and visually
% % passed the test
% % [h,p,ksstat,cv] = kstest(vx_g1) < using the Kolmogorov-Smirnoff test
% % indicates that none of these distributions are normal...
% 
% % vx_g1_mean = (mean(vx_2D(idx3_1_i:idx3_1_f),1))';
% % vx_g2_mean = (mean(vx_2D(idx3_2_i:idx3_2_f),1))';
% % vx_g3_mean = (mean(vx_2D(idx3_3_i:idx3_3_f),1))';
% 
% vx_g1_1 = vol_tx_year(idx2_1_ii:idx2_1_if,:)';
% vx_g1_2 = vol_tx_year(idx2_1_fi:idx2_1_ff,:)';
% vx_g1 = cat(2,vx_g1_1,vx_g1_2);
% vx_g2 = vol_tx_year(idx2_2_i:idx2_2_f,:)';
% 
% vx_group = cat(2,vx_g1,vx_g2);
% group = NaN(92,1);
% group(idx2_1_ii:idx2_1_if) = 1;
% group(idx2_1_fi:idx2_1_ff) = 1;
% group(idx2_2_i:idx2_2_f) = 2;
% 
% [p,tbl,stats] = anova1(vx_group,group);
% 
% [p,tbl,stats] = kruskalwallis(vx_group,group);
% 
% % at least with two groups, these clusters pass both tests, in that they
% % are distinct populations

%% Looking at Auto-Correlation for the volumetric transport
[r,p] = corr(vx_yr_anom'); % this looks at the linear correlation between 
% the different areas on the shelf as we move through time (for all 17
% years)

figure
hold on
imagesc(dist_nord_km,dist_nord_km,r)

set(gca,'Ydir','reverse');
set(gca,'Xdir','reverse');
colormap(cmocean('balance','pivot',0))
colorbar
clim([-0.5 1])
ylabel('Distance from Nordostrundingen [km]','FontSize',12)
xlabel('Distance from Nordostrundingen [km]','FontSize',12)

title('VT Correlation along the 600m Isobath (Monthly, GLORYS)')

set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/vol tx monthly/corr_mo_res.jpg']);

%%
[r,p] = corr(vol_tx_600G'); % this looks at the linear correlation between 
% the different areas on the shelf as we move through time (for all 17
% years)

figure
hold on
imagesc(dist_nord_km,dist_nord_km,r)

set(gca,'Ydir','reverse');
set(gca,'Xdir','reverse');
colormap(cmocean('balance','pivot',0))
colorbar
clim([-1 1])
ylabel('Distance from Nordostrundingen [km]','FontSize',12)
xlabel('Distance from Nordostrundingen [km]','FontSize',12)

title('VT Correlation along the 600m Isobath (Yearly, GLORYS)')

set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/vol tx monthly/corr_yr_res.jpg']);

%% anomalies
[r,p] = corr(vol_tx_600G'); % this looks at the linear correlation between 
% the different areas on the shelf as we move through time (for all 17
% years)

figure
hold on
imagesc(dist_nord_km,dist_nord_km,r)

set(gca,'Ydir','reverse');
set(gca,'Xdir','reverse');
colormap(cmocean('balance','pivot',0))
colorbar
clim([-1 1])
ylabel('Distance from Nordostrundingen [km]','FontSize',12)
xlabel('Distance from Nordostrundingen [km]','FontSize',12)

title('VT Correlation along the 600m Isobath (Yearly, GLORYS)')

set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS quiver/vol tx monthly/corr_yr_res.jpg']);
%% Revisiting a Volumetric Transport Budget
% adding bottom and top of polygon

% % finding the closest coordinates at 74N
% w_74i = find(abs(lon2D-lon_74)==nanmin(abs(lon2D-lon_74)));
% e_74i = find(abs(lon2D-polySE_lon)==nanmin(abs(lon2D-polySE_lon)));
% w_74 = lon2D(w_74i);
% e_74 = lon2D(e_74i);
% 
% % finding closest coordinates at Nord lon
% w_nordi = find(abs(lon2D-lon_nord)==nanmin(abs(lon2D-lon_nord)));
% e_nordi = find(abs(lon2D-polyNE_lon)==nanmin(abs(lon2D-polyNE_lon)));
% w_nord = lon2D(w_nordi);
% e_nord = lon2D(e_nordi);
% 
% % closest latitude values
% n_nord = lat2D(idxlatN);
% s_74 = lat2D(idxlatS);

% testing
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
% hold on
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
% %m_pcolor(lonbath,latbath,bath')
% % m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
% m_plot(lon_600,lat_600,'-k','LineWidth',1)
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
%     'ytick',60:5:80,'xtick',-35:15:5);
% % m_plot(ll_600_poly(:,1),ll_600_poly(:,2),'-g','LineWidth',1)
% m_scatter([e_74,w_74,e_nord,w_nord],[s_74,s_74,n_nord,n_nord],'*r')

%% 
% applplying the polygon limits to the v velocities that are flowing out of
% the polygon and further down the shelf. Assessing the "offshore" tranport
% at this latitude to create more of a volume budget. This is a fairly
% preliminary volume budget...just trying it out

% v_74 = v_600_avg(w_74i:e_74i,idxlatS);
% v_offi = find(v_74<0);  % assuming that we're going straight across the 74
% % degree latitude line then all the -v will be considered to be "across" the
% % shelf. In reality, this flow is going outside of our polygon
% 
% lon_74_span = lon2D(w_74i:e_74i);
% dlon_km = NaN(size(lon_74_span));
% for nn = 1:length(lon_74_span)-1
%     dlon_km(nn) = haversine([s_74 lon_74_span(nn)],[s_74 lon_74_span(nn+1)]);
% end
% dlon_km(end) = dlon_km(end-1);
% dlon_m = 1000*dlon_km;
% 
% vx_74_out_ind = 600*v_74(v_offi).*dlon_m(v_offi);
% vx_74_out_total = sum(vx_74_out_ind)/10^6;
% 
% vx_out_poly = abs(vx_74_out_total-offshelf_vx_tot);
