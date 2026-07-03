% RARE T-S Diagram
% Created by: Sara Vianco
% Masters Thesis
% 16 Feb 2024

%% Preliminaries
close all
clear all

cd('/Volumes/TERABYTE/RARE/matching_theo/')

addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/') % for cmocean
addpath('/Users/saravianco/Documents/Research/MATLAB/m_map/');
addpath('/Users/saravianco/Documents/Research/MATLAB/');
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/html/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/library/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/thermodynamics_from_t/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/pdf/");
addpath('/Users/saravianco/Documents/Research/MATLAB/GLORYS/GLORYS quiver/');
addpath('/Users/saravianco/Documents/MATLAB/haversine/') % this is for haversine
addpath('/Users/saravianco/Documents/MATLAB/FACS/') % this is for the density scatterplot

%% initial info/variable establishment
info = ncinfo("rare1.15.2_mn_ocean_reg_2003.nc");

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
        psal(:,:,:,:,year-2002) = ncread(full_file_path,'salt',start3D,count3D);
        ptemp(:,:,:,:,year-2002) = ncread(full_file_path,'temp',start3D,count3D);
        time(:,year-2002) = ncread(full_file_path, 'time');

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

%% reshaping salinity and temperature
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

%% Determining the angle of the bathymetry
% C gives me the x and y coordinates of the 600m contour line
bath_600 = bath;
% manually reworking the bathymetry matrix to yield only the outer shelf
% 600m contour
bath_600(1:1057,1:422) = 0;
bath_600(1501:1801,1:62) = 1000;
bath_600(241:313,680:end) = 1000;
figure
[C,h] = contour(lonbath,latbath,bath_600',[600 600]);
x_600 = C(1,2:end);
y_600 = C(2,2:end);

% smoothing the x and y data (still in x and y coordinates
% to lat and lon...)
lat_600_sm = smooth(y_600,100)';
lon_600_sm = smooth(x_600,100)';

%% Identifying Polygon coordinates: Nordostrundingen latitude
% [x,y] = ginput(1);
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
load('nord_xy.mat');
[polyNW_lon,lat_nord] = m_xy2ll(x,y);

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
[polySW_74,lat_74] = m_xy2ll(x,y);

%% Curtailed 600m isobath

% cutting down the 600m isobath contour to only contain the section from
% 74N to Nord lat
lat_600_poly = lat_600_sm(1,nord_600_c:seven4_600_c);
lon_600_poly = lon_600_sm(1,nord_600_c:seven4_600_c);

%% Identifying Polygon coordinates: Nordostrundingen latitude
% [x,y] = ginput(1);
load('nord_xy.mat');
[polyNW_lon,lat_nord] = m_xy2ll(x,y);

% attempting to determine the intersection between the northern boundary
% and the 600m isobath
int_0 = lat_600_sm-lat_nord;

nord_600 = min(abs(int_0));
[nord_600_r,nord_600_c] = find(abs(int_0)==nord_600);

polyNE_lon = lon_600_sm(nord_600_r,nord_600_c);
polyNE_lat = lat_600_sm(nord_600_r,nord_600_c);

%% ID-ing cordinates: 74N latitude
% attempting to determine the intersection between the northern boundary
% and the 600m isobath
int_0 = lat_600_sm-74;

seven4_600 = min(abs(int_0));
[seven4_600_r,seven4_600_c] = find(abs(int_0)==seven4_600);

polySE_lon = lon_600_sm(seven4_600_r,seven4_600_c);
polySE_lat = lat_600_sm(seven4_600_r,seven4_600_c);

%% Finding the coordinates along the 600m isobath

% for the original coordinate system (lon3D and lat3D) I need to identify
% the coordinates of the 600m isobath, roughly

% first to constrain the original latitudes by the polygon limits
idxlatN = find(abs(lat3D-polyNE_lat)==nanmin(abs(lat3D-polyNE_lat)),1,'first');
idxlatS = find(abs(lat3D-polySE_lat)==nanmin(abs(lat3D-polySE_lat)),1,'first');
idxlonW = find(abs(lon3D-polySE_lon)==nanmin(abs(lon3D-polySE_lon)),1,'first');
idxlonE = find(abs(lon3D-polyNE_lon)==nanmin(abs(lon3D-polyNE_lon)),1,'first');

% The northern and southern limits of the area we want to look at have been
% set. 
%% a different method to match the 600m coordinates

% Reshape the grids into column vectors
LON_vector = LON(:);
LAT_vector = LAT(:);

% Combine model_lon_vector and model_lat_vector into a matrix for knnsearch
model_coordinates = [LON_vector, LAT_vector];

% Find the indices on the model grid closest to the bathymetric contour
indices_closest = knnsearch(model_coordinates, [lon_600_poly', lat_600_poly']);

% Extract corresponding latitude and longitude values
closest_lon_values = LON_vector(indices_closest);
closest_lat_values = LAT_vector(indices_closest);

% Combine longitude and latitude into a matrix
coordinates = [closest_lon_values, closest_lat_values];

% Find unique rows based on latitude
[unique_lat, idx, ~] = unique(coordinates(:, 2), 'stable');

% Extract corresponding longitude values
unique_lon = coordinates(idx, 1);

unique_coord = [unique_lon,unique_lat];

lon_600 = unique_lon;
lat_600 = unique_lat;

%% establishing the geographic distances
lon_km = NaN(length(lon_600)-1,1);

for jj = 1:(length(lon_600)-1)
lon_km(jj,1) = haversine([lat_600(jj),lon_600(jj)],[lat_600(jj+1),lon_600(jj+1)]); 
% haversign calculates the geometric distance in km between two coordinates
end

% lon_km(end,1) = lon_km(end-1,1);
lon_m = lon_km*1000;
lon_m(end+1,1) = lon_m(end,1);

% now to calculate the half-way distance for the volumetric transports
lon_diff = NaN(size(lon_m))';

for kk = 1:length(lon_m)-1
lon_diff(kk) = (lon_m(kk+1)/2)+(lon_m(kk)/2);
end
lon_diff(end) = lon_m(end)/2;

%% finding the vector angles
[x,y]  = m_ll2xy(lon_600_poly, lat_600_poly);

clear dx
clear dy

x_600 = x(idx);
y_600 = y(idx);

% calculating a simple slope 
for m = 2:(length(lon_600)-1)
dy(m) = (y_600(m+1)-y_600(m-1));
dx(m) = (x_600(m+1)-x_600(m-1));
end

ang_rough = atan2d(dy(2:end),dx(2:end));
vec_ang = NaN(1,77);
vec_ang(2:end-1) = -ang_rough;
vec_ang(1) = vec_ang(2);
vec_ang(end) = vec_ang(end-1);

%% Matching the unique lon and lat indicies of the 600m isobath

for i = 1:length(idx)
idx_lon(i) = find(abs(lon3D-unique_lon(i))==nanmin(abs(lon3D-unique_lon(i))));
idx_lat(i) = find(abs(lat3D-unique_lat(i))==nanmin(abs(lat3D-unique_lat(i))));
end

%% Reducing the temp and salinity data to only include the area within the polygon

% getting rid of regions north of Nord
psal_4D(:,idxlatN:end,:,:) = NaN;
% and south of 74
psal_4D(:,1:idxlatS,:,:) = NaN;

ptemp_4D(:,idxlatN:end,:,:) = NaN;
ptemp_4D(:,1:idxlatS,:,:) = NaN;

%% NaN-ing regions outside of the 600m isobath
for p = 1:length(idx_lat)
        psal_4D(idx_lon(p):end,idx_lat(p),:,:) = NaN;
        ptemp_4D(idx_lon(p):end,idx_lat(p),:,:) = NaN;
end

%% test figure
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
m_contourf(LON,LAT,psal_4D(:,:,10,200)');
cmocean('haline')
colorbar
clim([25 34])
% m_plot(lon_600,lat_600,'-g','LineWidth',1)
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);

%% calculating the density
p = NaN(length(lat3D),length(z_star));

for i = 1:length(lat3D)
p(i,:) = gsw_p_from_z(-z_star,lat3D(i));
end
p = p';


SA = NaN(size(psal_4D));
for i = 1:204
    for j = 1:length(z_star)
[SA(:,:,j,i), in_ocean] = gsw_SA_from_SP(psal_4D(:,:,j,i),p(j,:),lon3D,lat3D);
    end
end

CT = NaN(size(ptemp_4D));
for i = 1:204
    for j = 1:length(z_star)
CT(:,:,j,i) = gsw_CT_from_pt(SA(:,:,j,i),ptemp_4D(:,:,j,i));
    end
end

sigma0 = NaN(size(ptemp_4D));
for i = 1:204
    for j = 1:length(z_star)
sigma0(:,:,j,i) = gsw_sigma0(SA(:,:,j,i),CT(:,:,j,i));
    end
end

% rho_anom = rho_CT_exact-1000;

%% Reducing the psal and the ptemp variables

% only using values in the top 600 m of the water column and taking a time
% average

psal_3D = nanmean(psal_4D(:,:,1:21,:),4);
ptemp_3D = nanmean(ptemp_4D(:,:,1:21,:),4);

%% creating the density contours
% finding the min and max values of temperature and salinity
mint = min(CT,[],"all");
maxt = max(CT,[],"all");

mins = min(SA,[],"all");
maxs = max(SA,[],"all");

% creating an evanly spaced temp and salinity vector
tempL = linspace(mint-1,maxt+1,156);
saltL = linspace(mins-1,maxs+1,156);

[CT_ref,SA_ref] = meshgrid(tempL,saltL);

sigma0_ref = gsw_sigma0(SA_ref,CT_ref);
sigma0_5 = (gsw_rho(SA_ref,CT_ref,500))-1000;

%% 
dist_nord = cumsum(lon_m);
dist_nord_km = dist_nord/1000;

%% Dividing up the TS diagram by the different onshelf flow regimes %%

%% Interpolating the depths to avoid averaging

z_interp = [0:5:600];

% Interpolate the matrix onto the new grid
SA_interp = NaN(400,120,121,204);
CT_interp = NaN(400,120,121,204);

for jj = 1:204
SA_interp(:,:,:,jj) = interp3(lat3D',lon3D,z_star,SA(:,:,:,jj),lat3D',lon3D,z_interp);
CT_interp(:,:,:,jj) = interp3(lat3D',lon3D,z_star,CT(:,:,:,jj),lat3D',lon3D,z_interp);
end

%% Reducing the psal and the ptemp variables

% only using values in the top 600 m of the water column and taking a time
% average

SA_3D = nanmean(SA_interp,4);
CT_3D = nanmean(CT_interp,4);

%% Plotting the TS Diagram

figure('color','w','position',[1300 50 1500 900],'papersize',[9 15]);
hold on
[C,h] = contour(SA_ref,CT_ref,sigma0_ref,'Color',[0.5 0.5 0.5]);
clabel(C,h)

c = cmocean('thermal',121);

for ii = 1:121
scatter(SA_3D(:,:,ii),CT_3D(:,:,ii),10,c(ii,:),'o','filled')
end
cmocean('thermal')
x = colorbar('Direction','reverse');
clim([0 600])

% plotting the water masses based on Wilcox 2023
%scatter([35,35,34.8,34,31.4],[6.6,4,1.4,-1.5,-1.5],100,'k','filled');

xlim([30 36])
ylim([-2.3 7])

xlabel('Absolute Salinity [PSU]','FontSize',18)
ylabel('Conservative Temperature [^oC]','FontSize',18)
ylabel(x,'Depth [m]','FontSize',18)
title('2003-2019 RARE Averages','FontSize',20)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/TS Diagrams/ts_avg204.jpg']);

%% Upper Section: 81.65N -> 79.65N
% this is still in the upper 600 m and on the shelf, so we should be
% capturing the entirety of the water column

SA_3D_upper = mean(SA_3D(:,idx_lat(1:21),:),4);
CT_3D_upper = mean(CT_3D(:,idx_lat(1:21),:),4);

figure('color','w','position',[1300 50 1500 900],'papersize',[9 15]);
hold on
[C,h] = contour(SA_ref,CT_ref,sigma0_ref,'Color',[0.5 0.5 0.5]);
clabel(C,h)

c = cmocean('thermal',121);

for ii = 1:121
scatter(SA_3D_upper(:,:,ii),CT_3D_upper(:,:,ii),10,c(ii,:),'o','filled')
end
cmocean('thermal')
x = colorbar('Direction','reverse');
clim([0 600])

% %scatter([35,35,34.8,34,31.4],[6.6,4,1.4,-1.5,-1.5],100,'k','filled');

xlim([30 36])
ylim([-2.3 7])

xlabel('Absolute Salinity [PSU]','FontSize',18)
ylabel('Conservative Temperature [^oC]','FontSize',18)
ylabel(x,'Depth [m]','FontSize',18)

title('2003-2019 RARE Averages: Top Polygon Section ','FontSize',20)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/TS Diagrams/ts_avg_upper.jpg']);

%% Middle Section: 79.55N -> 77.05N
% this is still in the upper 600 m and on the shelf, so we should be
% capturing the entirety of the water column

SA_3D_middle = mean(SA_3D(:,idx_lat(22:47),:),4);
CT_3D_middle = mean(CT_3D(:,idx_lat(22:47),:),4);

figure('color','w','position',[1300 50 1500 900],'papersize',[9 15]);
hold on
[C,h] = contour(SA_ref,CT_ref,sigma0_ref,'Color',[0.5 0.5 0.5]);
clabel(C,h)

c = cmocean('thermal',121);

for ii = 1:121
scatter(SA_3D_middle(:,:,ii),CT_3D_middle(:,:,ii),10,c(ii,:),'o','filled')
end
cmocean('thermal')
x = colorbar('Direction','reverse');
clim([0 600])

% %scatter([35,35,34.8,34,31.4],[6.6,4,1.4,-1.5,-1.5],100,'k','filled');

xlim([30 36])
ylim([-2.3 4])

xlabel('Absolute Salinity [PSU]','FontSize',18)
ylabel('Conservative Temperature [^oC]','FontSize',18)
ylabel(x,'Depth [m]','FontSize',18)

title('2003-2019 RARE Averages: Middle Polygon Section ','FontSize',20)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/TS Diagrams/ts_avg_middle.jpg']);

%% Bottom Section: 76.95N -> 74.05N
% this is still in the upper 600 m and on the shelf, so we should be
% capturing the entirety of the water column

SA_3D_bottom = mean(SA_3D(:,idx_lat(48:end),:),4);
CT_3D_bottom = mean(CT_3D(:,idx_lat(48:end),:),4);

figure('color','w','position',[1300 50 1500 900],'papersize',[9 15]);
hold on
[C,h] = contour(SA_ref,CT_ref,sigma0_ref,'Color',[0.5 0.5 0.5]);
clabel(C,h)

c = cmocean('thermal',121);

for ii = 1:121
scatter(SA_3D_bottom(:,:,ii),CT_3D_bottom(:,:,ii),10,c(ii,:),'o','filled')
end
cmocean('thermal')
x = colorbar('Direction','reverse');
clim([0 600])

%scatter([35,35,34.8,34,31.4],[6.6,4,1.4,-1.5,-1.5],100,'k','filled');

xlim([30 36])
ylim([-2.3 4])

xlabel('Absolute Salinity [PSU]','FontSize',18)
ylabel('Conservative Temperature [^oC]','FontSize',18)
ylabel(x,'Depth [m]','FontSize',18)

title('2003-2019 RARE Averages: Bottom Polygon Section ','FontSize',20)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/TS Diagrams/ts_avg_bottom.jpg']);

%% Now to pivot to individual depths
% I'm going to try and separate the water column into individual depths
% determining the spacing between the different depths

% first to determine a rough estimate for the ekman depth
omega = (2*pi)/86400;
f = 2*omega*sind(78);
D_ek = sqrt(2*0.05/f);

%%
% separating the SA/CT variables
SA_3D_50 = mean(SA_3D(:,:,1:11),3,'omitnan');
SA_3D_200 = mean(SA_3D(:,:,12:41),3,'omitnan');
SA_3D_600 = mean(SA_3D(:,:,42:end),3,'omitnan');
SA_3D_full = mean(SA_3D,3,'omitnan');

CT_3D_50 = mean(CT_3D(:,:,1:11),3,'omitnan');
CT_3D_200 = mean(CT_3D(:,:,12:41),3,'omitnan');
CT_3D_600 = mean(CT_3D(:,:,42:end),3,'omitnan');
CT_3D_full = mean(CT_3D,3,'omitnan');

%% Plotting figures based on depths (reusing code for the different depths)
addpath('/Users/saravianco/Documents/MATLAB/densityplot/') % this is for the density scatterplot

% Define the number of bins for the 2D histogram
num_bins = 50;

% Create 2D histogram to estimate density
[counts, temperature_edges, salinity_edges] = histcounts2(CT_3D_50(:),...
    SA_3D_50(:), num_bins, 'Normalization', 'count');

% Compute density values from the counts
density = counts / numel(CT_3D_50) * 100; % Convert to percentage

% Flatten the data for scatter plot
temperature_flat = CT_3D_50(:);
salinity_flat = SA_3D_50(:);

% Convert temperature and salinity values to bin indices
temperature_bins = discretize(temperature_flat, temperature_edges);
salinity_bins = discretize(salinity_flat, salinity_edges);

% Remove NaN indices
valid_indices = ~isnan(temperature_bins) & ~isnan(salinity_bins);
temperature_bins = temperature_bins(valid_indices);
salinity_bins = salinity_bins(valid_indices);

%%
figure('color','w','position',[1300 50 1500 900],'papersize',[9 15]);
hold on
[C,h] = contour(SA_ref,CT_ref,sigma0_ref,'Color',[0.5 0.5 0.5]);
clabel(C,h,'FontSize',20,'Color',[0.5 0.5 0.5])

% Create scatter plot with density-based coloring
scatter(salinity_flat(valid_indices), temperature_flat(valid_indices),10,...
    density(sub2ind(size(density), temperature_bins, salinity_bins)), 'filled');
x = colorbar;
cmocean('thermal')
clim([0 0.3])

xlim([30 35.5])
ylim([-2.3 3])

% create the contours for the different water masses
% Polar Surface Water
yline(0,'-k','0','LineWidth',1)
% line([x_start, x_end], [y_value, y_value], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
yline(2,'-k','2','LineWidth',1)

% limiting the contours
% x_mask = (SA_ref >= x_start) & (X <= x_end);

[j,k] = contour(SA_ref,CT_ref,sigma0_ref,[27.7 27.7],'-k','LineWidth',1);
clabel(j,k,'FontSize',20)
[j,k] = contour(SA_ref,CT_ref,sigma0_ref,[27.97 27.97],'-k','LineWidth',1);
clabel(j,k,'FontSize',20)
[j,k] = contour(SA_ref,CT_ref,sigma0_ref,[29.97 29.97],'-k','LineWidth',1);
clabel(j,k,'FontSize',20)
[l,m] = contour(SA_ref,CT_ref,sigma0_5,[30.444 30.444],'-k','LineWidth',1);
clabel(l,m,'FontSize',20)

ax = gca; 
ax.FontSize = 20; 

xlabel('Absolute Salinity [PSU]','FontSize',25)
ylabel('Conservative Temperature [^oC]','FontSize',25)
ylabel(x,'Relative Percentage of Sample [%]','FontSize',25)
title('2003-2019 RARE Averages: Ekman Layer (0-50 m) ','FontSize',25)
% title('2003-2019 RARE Averages: Sub-Ekman Shelf Layer (50-200 m) ','FontSize',20)
% title('2003-2019 RARE Averages: Shelfbreak Layer (200-600 m) ','FontSize',20)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/TS Diagrams/ts_ekman.jpg']);
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/TS Diagrams/ts_shelf.jpg']);
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/TS Diagrams/ts_shelfbreak.jpg']);
