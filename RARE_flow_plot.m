% RARE Flow Plot
% Created by: Sara Vianco
% Masters Thesis
% 20 May 2024

%% Preliminaries
close all
clear all

cd('/Volumes/TERABYTE/RARE/matching_theo/')
% cd('/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_quiver/')
% load('/Volumes/One Touch/Research Datasets/RARE/matching_theo/RARE_quiver_vars.mat');

addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/') % for cmocean
addpath('/Users/saravianco/Documents/Research/MATLAB/m_map/');
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/html/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/library/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/thermodynamics_from_t/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/pdf/");
addpath('/Users/saravianco/Documents/Research/MATLAB/GLORYS/GLORYS quiver/');
addpath('/Users/saravianco/Documents/MATLAB//haversine/') % this is for haversine

addpath('/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_VX_corrected/')

%% initial info/variable establishment
info = ncinfo("rare1.15.2_mn_ocean_reg_2003.nc");
info2 = ncinfo("rare1.15.2_mn_ocean_reg_2019.nc");

time = ncread('rare1.15.2_mn_ocean_reg_2003.nc','time');
lon = ncread('rare1.15.2_mn_ocean_reg_2003.nc','xt_ocean');
lat = ncread('rare1.15.2_mn_ocean_reg_2003.nc','yt_ocean');
z_star = ncread('rare1.15.2_mn_ocean_reg_2003.nc','sw_ocean');
z_star = z_star';

%% Establishing 3D Variables
% determining lat and lon limits for figure
latS_map = 73;
latN_map = 85;
lonW_map = -35;
lonE_map = 5;


%% testing out a new way to read in the files

% finding these lat/lon ranges in lat/lon variables
lat3D_i = find(lat>latS_map & lat<latN_map);
lon3D_i = find(lon<lonE_map & lon>lonW_map);
start3D = [lon3D_i(1) lat3D_i(1) 1 1];
count3D = [length(lon3D_i) length(lat3D_i) length(z_star) 12];

lon3D = lon(lon3D_i);
lat3D = lat(lat3D_i);

% Define the path to the directory where your netCDF files are located
directory_path = '/Volumes/TERABYTE/RARE/matching_theo/';

% Loop over the years from 2003 to 2017
for year = 2003:2019
    % Construct the filename based on the format
    filename = sprintf('rare1.15.2_mn_ocean_reg_%d.nc', year);

    % Create the full file path
    full_file_path = fullfile(directory_path, filename);

    % Check if the file exists
    if exist(full_file_path, 'file')
        % Read the netCDF file using ncread
        v(:,:,:,:,year-2002) = ncread(full_file_path,'v',start3D,count3D); 
        u(:,:,:,:,year-2002) = ncread(full_file_path,'u',start3D,count3D);
        time(:,year-2002) = ncread(full_file_path, 'time');
        ptemp(:,:,:,:,year-2002) = ncread(full_file_path,'temp',start3D,count3D);
        psal(:,:,:,:,year-2002) = ncread(full_file_path,'salt',start3D,count3D);

        % Perform operations on the data here
    else
        fprintf('File %s does not exist.\n', filename);
    end
end

[LON,LAT] = meshgrid(lon3D,lat3D);

%% reshaping u and v
size_u = size(u);
u_per = permute(u, [1, 2, 3, 4, 5]);
v_per = permute(v, [1, 2, 3, 4, 5]);

size_4D = [size_u(1:3), prod(size_u(4:5))];

u_4D = reshape(u_per, size_4D);
v_4D = reshape(v_per, size_4D);

%% reshaping practical salinity and potential temperature

size_psal = size(psal);
psal_per = permute(psal, [1, 2, 3, 4, 5]);
ptemp_per = permute(ptemp, [1, 2, 3, 4, 5]);

size_4D = [size_psal(1:3), prod(size_psal(4:5))];

psal_4D = reshape(psal_per, size_4D);
ptemp_4D = reshape(ptemp_per, size_4D);


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


%% adjusting the time

% making the date number into one that is readable by MATLAB, forming a
% date string afterwards

base = datenum(1979,12,31,0,0,0);
formatIn = 'mm/dd/yy';
   
date_num = datenum(time + base);
date_string = datestr(date_num,'mm/dd/yy');
dv = datevec(date_string,formatIn);

%% Interpolating the velocities to standardize the averge velocities in the upper 600 m
u_3D = nanmean(u_4D,4);
v_3D = nanmean(v_4D,4);
ptemp_3D = nanmean(ptemp_4D,4);

z_interp = [0:25:600];

for i = 1:length(lon3D)
    for j = 1:length(lat3D)
u_3D_interp(i,j,:) = interp1(z_star,squeeze(u_3D(i,j,:)),z_interp);
    end
end

for i = 1:length(lon3D)
    for j = 1:length(lat3D)
v_3D_interp(i,j,:) = interp1(z_star,squeeze(v_3D(i,j,:)),z_interp);
    end
end

for i = 1:length(lon3D)
    for j = 1:length(lat3D)
ptemp_3D_interp(i,j,:) = interp1(z_star,squeeze(ptemp_3D(i,j,:)),z_interp);
    end
end

%%
u_600_avg = nanmean(u_3D_interp,3);
v_600_avg = nanmean(v_3D_interp,3);
ptemp_50_avg = nanmean(ptemp_3D_interp,3);

speed_600_avg = sqrt(u_600_avg.^2+v_600_avg.^2);

%% Plotting the location of the vertical section plots

load('RARE_vx_vars.mat')

u_ref = NaN(size(LON));
v_ref = NaN(size(LON));
u_ref(41,50) = 0.25;
v_ref(41,50) = 0;

figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);

% m_pcolor(LON,LAT,speed_600_avg');
% cmocean('speed')
% x = colorbar;
% clim([0 0.2]);

% m_pcolor(LON,LAT,ptemp_50_avg');
% cmocean('thermal')
% x = colorbar;

% clim([-2 3]);

% m_contour(lonbath,latbath,bath',[100 350 600],'linew',1.5,'color',[0.5 0.5 0.5]);
m_pcolor(lonbath,latbath,bath');
x = colorbar;
set(x, 'YDir', 'reverse' );

cmocean('deep')
clim([0 5000])
x.Ticks = 0:1000:5000;

m_contour(lonbath,latbath,bath',[100 200 250 350 600 1000 2000 5000],'linew',1.5,'color',[0.5 0.5 0.5]);


% m_quiver(LON(2:1:end,35:7:395),LAT(2:1:end,35:7:395),u_600_avg(35:7:395,2:1:end)',...
%     v_600_avg(35:7:395,2:1:end)',2.5,'k','LineWidth',0.5,'Clipping','on');

% m_scatter(lon_nord_span,n_nord_lat,50,vx_nord_mean*1000,'Filled');
% cmocean('balance','pivot',0)
% clim([-100 100])
% 
% m_scatter(lon_74_span,s_74_lat,50,-vx_74_mean*1000,'Filled');
% cmocean('balance','pivot',0)
% clim([-100 100])
% 
% m_scatter(lon_600,lat_600,50,vx_avg*1000,'Filled');
% h = colorbar;
% clim([-100 100])
% cmocean('balance','pivot',0)

% 
% ylabel(x,'Potential Temperature [°C]')

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
% set(x,'fontsize',10,'fontname','helvetica','position',[0.7 0.55 0.01 0.3],'tickdir','out','ytick',0:0.05:0.2);
% m_quiver(LON,LAT,u_ref,v_ref,1,'k','LineWidth',0.7,'MaxHeadSize',2);

% m_line(-22:-13,70.95,'linewi',2,'color','r');
% m_line(-20:-8,73.95,'linewi',2,'color','r');
% m_line(-19:0,76.95,'linewi',2,'color','r');
% m_line(-18:0,78.55,'linewi',2,'color','r');
% m_line(-18:0,80.05,'linewi',2,'color','r');
% m_line(-13:0,81.45,'linewi',2,'color','r');

m_plot(-8:-2,78.5,'linewi',2,'color','r');

% dim = [0.38 0.35 .045 .06];
% str = '0.25 m/s';
% annotation('textbox',dim,'String',str,'FontSize',13,'FitBoxToText','off','Rotation',-14.5);

ylabel(x,['Depth [m]'],'FontSize',14)
%title('NEGS Region','FontSize',16)

set(gcf,'color','w');
imwrite(getframe(gcf).cdata,['/Users/saravianco/Documents/Research/Seperate Images/NEGS_bathy.jpg']);

%%
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[75 82],'lons',[-22 0]);
m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 200 250 350 600 1000 2000 5000],'linew',1.5,'color',[0.5 0.5 0.5]);
cmocean('deep')
h = colorbar;
set(h, 'YDir', 'reverse' );

clim([0 600])
h.Ticks = 0:200:600;

% m_plot([lon_79 lon_79],[lat_79S lat_79N],'linewi',2,'color','r');
% m_plot([-21 -19],[79.583 79.583],'linewi',2,'color','r');
m_plot(-8:-2,78.5,'linewi',2,'color','r');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
ylabel(h,'Depth [m]','FontSize',14)
% title('79NG Section','FontSize',16)
imwrite(getframe(gcf).cdata,['/Users/saravianco/Documents/Research/Seperate Images/NEGS_bathy_zoom.jpg']);

