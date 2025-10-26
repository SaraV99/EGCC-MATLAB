% GLORYS Volume Transport re-worked with correct weighting for the
% changing depth levels
% Created by: Sara Vianco
% Masters Thesis
% 31 DEC 2023

%% preliminaries
close all
clear all

cd('/Volumes/TERABYTE/GLORYS/')
addpath(genpath('/Volumes/TERABYTE/GLORYS/'));
addpath('/Users/saravianco/Documents/Research/MATLAB/GLORYS/')
addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/') % for cmocean
addpath('/Users/saravianco/Documents/Research/MATLAB/m_map/');
addpath('/Users/saravianco/Documents/MATLAB//haversine/') % this is for haversine
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/html/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/library/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/thermodynamics_from_t/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/pdf/");
addpath('/Users/saravianco/Documents/Research/MATLAB/GLORYS/GLORYS quiver');
addpath('/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_quiver');

% this was not preliminarily added until later when I went to compare the
% Volume transport between RARE and GLORYS
% load('date_num.mat')


%% 
% gaining relevant info
info = ncinfo('mercatorglorys12v1_gl12_mean_200301.nc');

% single variables
lat = ncread('mercatorglorys12v1_gl12_mean_200301.nc','latitude');
lon = ncread('mercatorglorys12v1_gl12_mean_200301.nc','longitude');
z = ncread('mercatorglorys12v1_gl12_mean_200301.nc','depth');

%% Establishing 3D Variables
% determining lat and lon limits for figure
latS_map = 73;
latN_map = 85;
lonW_map = -35;
lonE_map = 5;

% finding these lat/lon ranges in lat/lon variables
lat3D_i = find(lat>latS_map & lat<latN_map);
lon3D_i = find(lon<lonE_map & lon>lonW_map);

start3D = [lon3D_i(1) lat3D_i(1) 1 1];
count3D = [length(lon3D_i) length(lat3D_i) length(z) 1];

lon3D = lon(lon3D_i);
lat3D = lat(lat3D_i);

[LON,LAT] = meshgrid(lon3D,lat3D);

%%
baseFolder = '/Volumes/TERABYTE/GLORYS/';

% Years range from 2003 to 2019
years = 2003:2019;

% Loop over each year
for year = years
    % Loop over each month
    for month = 1:12
        % Create the file name based on the format you provided
        fileName = sprintf('mercatorglorys12v1_gl12_mean_%04d%02d.nc', year, month);
        
        % Construct the full file path
        fullPath = fullfile(baseFolder, num2str(year), fileName);
        % Check if the file exists
        if exist(fullPath, 'file') == 2
            % Read data using ncread or your preferred method
            % For example, assuming you want to read a variable named 'temperature'
            pot_temp(:,:,:,month,year-2002) = ncread(fullPath,'thetao',start3D,count3D);
            prac_salinity(:,:,:,month,year-2002) = ncread(fullPath,'so',start3D,count3D);
            u(:,:,:,month,year-2002) = ncread(fullPath,'uo',start3D,count3D);
            v(:,:,:,month,year-2002) = ncread(fullPath,'vo',start3D,count3D);
            time_raw(month,year-2002) = ncread(fullPath,'time');
             else
            disp(['File not found: ', fullPath]);
        end
    end
end
            

%% reshaping u and v
size_u = size(u);
u_per = permute(u, [1, 2, 3, 4, 5]);
v_per = permute(v, [1, 2, 3, 4, 5]);

size_4D = [size_u(1:3), prod(size_u(4:5))];

u_4D = reshape(u_per, size_4D);
v_4D = reshape(v_per, size_4D);

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

%% Adjusting the time
time_pivot = datenum(1950,1,1,0,0,0);
time_matlab = (time_raw/24) + time_pivot;
time_matlab = floor(time_matlab);
date = datevec(time_matlab);

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

%% a testing figure for the polygon limits

figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
cmocean('deep')

b6 = m_plot(lon_600_poly,lat_600_poly,'-r','LineWidth',2);
m_plot([polyNW_lon polyNE_lon],[lat_nord lat_nord],'-r','LineWidth',2)
m_plot([polySW_74 polySE_lon],[polySE_lat polySE_lat],'-r','LineWidth',2)
m_plot(polyNE_lon,polyNE_lat,'*w')
m_plot(polySE_lon,polySE_lat,'*w')

% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
title('Shelf Polygon','FontSize',16)

%% Finding the coordinates of the quivers along the 600m isobath

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


% now to calculate the half-way distance for the Volume transports
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
vec_ang = NaN(1,length(idx));
vec_ang(2:end-1) = -ang_rough;
vec_ang(1) = vec_ang(2);
vec_ang(end) = vec_ang(end-1);

%% applying the unique lat and lon indicies to the u and v variables
% in order to do so, I first need to match the closest unique lon and lat
% values to the original lon and lat values

for i = 1:length(idx)
idx_lon(i) = find(abs(lon3D-unique_lon(i))==nanmin(abs(lon3D-unique_lon(i))));
idx_lat(i) = find(abs(lat3D-unique_lat(i))==nanmin(abs(lat3D-unique_lat(i))));
end

for j = 1:204
    for k = 1:50
u_600(:,k,j) = interp2(lat3D,lon3D,u_4D(:,:,k,j),lat3D(idx_lat),lon3D(idx_lon),'nearest');
    end
end

for j = 1:204
    for k = 1:50
v_600(:,k,j) = interp2(lat3D,lon3D,v_4D(:,:,k,j),lat3D(idx_lat),lon3D(idx_lon),'nearest');
    end
end

%% reshaping u_600 and v_600 so as to apply the rotation matrix

% Reshape the matrix into a single row vector
u_long = reshape(u_600, 1, []);
v_long = reshape(v_600, 1, []);

U = cat(1,u_long,v_long);

R = NaN(2,2,(length(u_long)));

vec_ang_mo = repmat(vec_ang,1,204*50);

for jj = 1:length(vec_ang_mo)
R(:,:,jj)=[cosd(vec_ang_mo(jj)), -sind(vec_ang_mo(jj)) ; sind(vec_ang_mo(jj)), cosd(vec_ang_mo(jj))];
U_r(:,jj) = R(:,:,jj)*U(:,jj);
end

clear u_600r v_600r
u_600r = U_r(1,:);
v_600r = U_r(2,:);

u_600r = reshape(u_600r,size(u_600));
v_600r = reshape(v_600r,size(v_600));

%% In order to make the GLORYS and RARE data/depths comparable...
% we need to interpolate a bit with the variables

z_int = NaN(51,1);
z_int(1:32) = z(1:32);
z_int(33) = 600;
z_int(34:end) = z(33:50);

% Interpolate the matrix onto the new grid
u_600ri = interp3(z,1:92,1:204,u_600r,z_int,1:92,1:204);
v_600ri = interp3(z,1:92,1:204,v_600r,z_int,1:92,1:204);

%% determining the onshore/off shore velocities
% curtailing v_600ri to only include the top 600m
v_600ri = v_600ri(:,1:33,:);

% I need to find the distances between the z cells at the midpoints
z_diff = NaN(size(z_int))';

for kk = 2:length(z_int)-1
z_diff(kk) = ((z_int(kk+1)+z_int(kk))/2)-((z_int(kk)+z_int(kk-1))/2);
end
z_diff(1) = ((z_int(2)+z_int(1))/2)-(z_int(1)/2);


% finding the indicies where the seasonal rotated v velocities are negative
% to diagnose onshelf transport
onshelf_idx = find(v_600ri<0);
v_onshelf = zeros(size(v_600ri));
v_onshelf(onshelf_idx) = v_600ri(onshelf_idx);

% this is calculating the individual Volume transport at each location
% along the shelf at each of the depths from surface to 600m

for i = 1:204
onshelf_vx(:,:,i) = z_diff(1:33).*squeeze(v_onshelf(:,:,i)).*lon_diff';
end

vx_sea = NaN(size(onshelf_vx));
vx_sea(onshelf_idx) = onshelf_vx((onshelf_idx));

for i = 1:204
    for j = 1:92
onshelf_600_sea(j,i) = (sum(onshelf_vx(j,:,i),2,'omitnan'))/10^6;
    end
end

% finding the indicies where the seasonal rotated v velocities are positive
% to diagnose offshelf transport
offshelf_idx = find(v_600ri>0);
v_offshelf = zeros(size(v_600ri));
v_offshelf(offshelf_idx) = v_600ri(offshelf_idx);

for i = 1:204
offshelf_vx(:,:,i) = z_diff(1:33).*squeeze(v_offshelf(:,:,i)).*lon_diff';
end

vx_sea(offshelf_idx) = offshelf_vx((offshelf_idx));

for i = 1:204
    for j = 1:92
offshelf_600_sea(j,i) = (sum(offshelf_vx(j,:,i),2,'omitnan'))/10^6;
    end
end

% calculating the Volume transport in the top 600m
for i = 1:204
    for j = 1:92
vx(j,i) = (sum(vx_sea(j,:,i),2,'omitnan'))/10^6;
    end
end

%% Imaging the onshore and offshore transport

dist_nord = cumsum(lon_m);
dist_nord_km = dist_nord/1000;

date_num_all = reshape(time_matlab,1,[]);

figure
imagesc(date_num_all,dist_nord_km,1000*vx)
set(gca,'Ydir','reverse');
dateFormat = 11;
datetick('x',dateFormat)
xlim([date_num_all(1,1) date_num_all(end)])
x = colorbar;
clim([-20 20])
colormap(cmocean('balance','pivot',0))
ylabel('Distance from Nordostrundingen [km]','FontSize',12)
ylabel(x,'Volume Transport [mSv]','FontSize',12)
xlabel('Year','FontSize',12)
title('Cross-shelf Volume Transport for 2003-2019 (GLORYS)','FontSize',14)
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_seasonal.jpg']);

%% Looking at VX anomalies

vx_avg = mean(vx,2);
vx_anom = vx - vx_avg;

figure
hold on
imagesc(date_num_all,dist_nord_km,1000*vx_anom);
set(gca,'Ydir','reverse');
dateFormat = 11;
datetick('x',dateFormat)
xlim([date_num_all(1,1) date_num_all(end)])
x = colorbar;
clim([-20 20])
colormap(cmocean('balance','pivot',0))
ylabel('Distance from Nordostrundingen [km]','FontSize',12)
ylabel(x,"Vol Tx' [mSv]",'FontSize',12)
title('Cross-shelf Volume Transport Anomalies for 2003-2019 (GLORYS)','FontSize',14)
set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_seasonal_anom.jpg']);

%%
years = [{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012',...
    '2013','2014','2015','2016','2017','2018','2019'}];
offshelf_sea_sum = sum(abs(offshelf_600_sea),1);
onshelf_sea_sum = sum(abs(onshelf_600_sea),1);
net_sea_vx = abs(offshelf_sea_sum)-abs(onshelf_sea_sum);
figure
hold on
plot([1:204],offshelf_sea_sum,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot([1:204],-onshelf_sea_sum,'-','Color',[0 0.4470 0.7410],'LineWidth',2)
plot([1:204],net_sea_vx,'-k','LineWidth',2)

legend('Offshelf Vol Tx','Onshelf Vol Tx','Net Vol Tx')
% ylim([])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)

title('Volume Transport Budget at the Outer Edge of the Polygon (GLORYS)')
ylabel('Vol Tx [Sv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_budget_shelf.jpg']);

%% Constructing a Volume Transport Budget
% adding bottom and top of polygon

% finding the closest coordinates at 74N
w_74i = find(abs(lon3D-polySW_74)==nanmin(abs(lon3D-polySW_74)));
e_74i = find(abs(lon3D-polySE_lon)==nanmin(abs(lon3D-polySE_lon)));
w_74_lon = lon3D(w_74i);
e_74_lon = lon3D(e_74i);

% finding closest coordinates at Nord lon
w_nordi = find(abs(lon3D-polyNW_lon)==nanmin(abs(lon3D-polyNW_lon)));
e_nordi = find(abs(lon3D-polyNE_lon)==nanmin(abs(lon3D-polyNE_lon)));
w_nord_lon = lon3D(w_nordi);
e_nord_lon = lon3D(e_nordi);

% closest latitude values
n_nord_lat = lat3D(idxlatN);
s_74_lat = lat3D(idxlatS);


%% Southern Edge of the Polygon
% applplying the polygon limits to the v velocities that are flowing out of
% the polygon and further down the shelf. Assessing the "offshore" tranport
% at this latitude to create more of a volume budget. This is a fairly
% preliminary volume budget...just trying it out

%reducing the depths to be the top 600m
v_74 = squeeze(v_4D(w_74i:e_74i,idxlatS,1:33,:));

lon_74_span = lon3D(w_74i:e_74i);
dxS_km = NaN(size(lon_74_span));
for nn = 2:length(lon_74_span)
    dxS_km(nn) = haversine([s_74_lat lon_74_span(nn-1)],[s_74_lat lon_74_span(nn)]);
end
dxS_km(1,1) = dxS_km(2,1);
dxS_m = 1000*dxS_km;
dxS_m(1) = dxS_m(1)/2;
dxS_m(end) = dxS_m(end)/2;

vx_600_S = NaN(size(v_74));
v_oni = find(v_74>0);  % assuming that we're going straight across the 74
% degree latitude line then all the +v will be considered to be "on" the
% shelf. In reality, this flow is going inside of our polygon
v_74_on = zeros(size(v_74));
v_74_on(v_oni) = v_74(v_oni);

for i = 1:204
onshelf_vx_74(:,:,i) = z_diff(:,1:33).*squeeze(v_74_on(:,:,i)).*dxS_m;
end

vx_600_S(v_oni) = onshelf_vx_74(v_oni);

v_offi = find(v_74<0);  
v_74_off = zeros(size(v_74));
v_74_off(v_offi) = v_74(v_offi);

for i = 1:204
offshelf_vx_74(:,:,i) = z_diff(:,1:33).*squeeze(v_74_off(:,:,i)).*dxS_m;
end
vx_600_S(v_offi) = offshelf_vx_74(v_offi);

% adding up all the depths to get the net vx at each gridcell across all
% months
for i = 1:204
    for j = 1:length(lon_74_span)
vx_74(j,i) = (sum(vx_600_S(j,:,i),2,'omitnan'))/10^6;
    end
end

%% Northern Edge of the Polygon

% applplying the polygon limits to the v velocities that are flowing out of
% the polygon and further down the shelf. Assessing the "offshore" tranport
% at this latitude to create more of a volume budget. This is a fairly
% preliminary volume budget...just trying it out

v_nord = squeeze(v_4D(w_nordi:e_nordi,idxlatN,1:33,:));

lon_nord_span = lon3D(w_nordi:e_nordi);

dxN_km = NaN(size(lon_nord_span));
for nn = 2:length(lon_nord_span)
    dxN_km(nn) = haversine([n_nord_lat lon_nord_span(nn-1)],[n_nord_lat lon_nord_span(nn)]);
end
dxN_km(1,1) = dxN_km(2,1);
dxN_m = 1000*dxN_km;
dxN_m(1) = dxN_m(1)/2;
dxN_m(end) = dxN_m(end)/2;

vx_600_N = NaN(size(v_nord));
v_oni = find(v_nord<0);  % assuming that we're going straight across the Nord
% latitude line then all the -v will be considered to be "on" the
% shelf. In reality, this flow is going inside of our polygon
v_nord_on = zeros(size(v_nord));
v_nord_on(v_oni) = v_nord(v_oni);

for i = 1:204
onshelf_vx_nord(:,:,i) = z_diff(:,1:33).*squeeze(v_nord_on(:,:,i)).*dxN_m;
end

vx_600_N(v_oni) = onshelf_vx_nord(v_oni);

v_offi = find(v_nord>0); % similarly, positive flow will be flowing outside
% of our polygon and thus will be considered offshelf

v_nord_off = zeros(size(v_nord));
v_nord_off(v_offi) = v_nord(v_offi);
for i = 1:204
offshelf_vx_nord(:,:,i) = z_diff(:,1:33).*squeeze(v_nord_off(:,:,i)).*dxN_m;
end
vx_600_N(v_offi) = offshelf_vx_nord(v_offi);

for i = 1:204
    for j = 1:length(lon_nord_span)
vx_nord(j,i) = (sum(vx_600_N(j,:,i),2,'omitnan'))/10^6;
    end
end

%% Summing up the Volume Transport to Construct the Volume Budget

% summing up 600m onshelf and offshelf VX at Nord: 2D
for i = 1:204
    for j = 1:length(lon_nord_span)
        onshelf_N_2D(j,i) = (sum(onshelf_vx_nord(j,:,i),2,'omitnan'))/10^6;
        offshelf_N_2D(j,i) = (sum(offshelf_vx_nord(j,:,i),2,'omitnan'))/10^6;
    end 
end

% summing up 600m onshelf and offshelf VX at 74: 2D
for i = 1:204
    for j = 1:length(lon_74_span)
        onshelf_S_2D(j,i) = (sum(onshelf_vx_74(j,:,i),2,'omitnan'))/10^6;
        offshelf_S_2D(j,i) = (sum(offshelf_vx_74(j,:,i),2,'omitnan'))/10^6;
    end 
end

% summing up 600m onshelf and offshelf VX at Nord: 1D
for i = 1:204
    onshelf_N_sea(i) = sum(onshelf_N_2D(:,i),1,'omitnan');
    offshelf_N_sea(i) = sum(offshelf_N_2D(:,i),1,'omitnan');
end

% summing up 600m onshelf and offshelf VX at 74: 1D
for i = 1:204
    onshelf_S_sea(i) = sum(onshelf_S_2D(:,i),1,'omitnan');
    offshelf_S_sea(i) = sum(offshelf_S_2D(:,i),1,'omitnan');
end

% now to add them together
net_onshelf_vx = abs(onshelf_S_sea)+abs(onshelf_N_sea)+abs(onshelf_sea_sum);
net_offshelf_vx = abs(offshelf_S_sea)+abs(offshelf_N_sea)+abs(offshelf_sea_sum);
net_vx = net_offshelf_vx-net_onshelf_vx;

test = mean(net_vx); % this is negative so it's net on

%%
figure
hold on
plot([1:204],net_offshelf_vx,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot([1:204],-net_onshelf_vx,'-','Color',[0 0.4470 0.7410],'LineWidth',2)
plot([1:204],net_vx,'-k','LineWidth',2)

legend('Offshelf Vol Tx','Onshelf Vol Tx','Net Vol Tx','location','southoutside')
% ylim([])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)

title('Volume Transport Budget of the 600m Polygon (GLORYS)')
ylabel('Vol Tx [Sv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_budget_600poly.jpg']);

%% creating the same figure but for the northern and southern edges of the polygon
%% Nord
net_N_vx = abs(offshelf_N_sea)-abs(onshelf_N_sea);

figure
hold on
plot([1:204],1000*abs(offshelf_N_sea),'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot([1:204],-1000*abs(onshelf_N_sea),'-','Color',[0 0.4470 0.7410],'LineWidth',2)
plot([1:204],1000*net_N_vx,'-k','LineWidth',2)

legend('Offshelf Vol Tx','Onshelf Vol Tx','Net Vol Tx','location','southeast')
% ylim([])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)

title('Volume Transport Budget: Northern Edge of Polygon (GLORYS)')
ylabel('Vol Tx [mSv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_budget_Nord.jpg']);

test_N = mean(net_N_vx);

%% 74N
net_S_vx = abs(onshelf_S_sea)-abs(offshelf_S_sea);

figure
hold on
plot([1:204],-abs(offshelf_S_sea),'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot([1:204],abs(onshelf_S_sea),'-','Color',[0 0.4470 0.7410],'LineWidth',2)
plot([1:204],net_S_vx,'-k','LineWidth',2)

legend('Offshelf Vol Tx','Onshelf Vol Tx','Net Vol Tx','location','southoutside')
% ylim([])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)

title('Volume Transport Budget: Southern Edge of Polygon (GLORYS)')
ylabel('Vol Tx [Sv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_budget_74.jpg']);

test_S = mean(net_S_vx);

%% Plotting the Net VX on each of the 3 Polygon Sides 
% (all 3 black lines on the same plot)

% calculating the mean of all of the Volume budgets on the polygon

mean_vx_600 = nanmean(net_sea_vx);
mean_vx_N = nanmean(net_N_vx);
mean_vx_S = nanmean(net_S_vx);
sum_vx = mean_vx_600+mean_vx_N-mean_vx_S;

figure
hold on
plot([1:204],net_sea_vx,'-','Color',[0 0 0],'LineWidth',2)
plot([1:204],net_N_vx,'-','Color',[0.35 0.35 0.35],'LineWidth',2)
plot([1:204],-net_S_vx,'-','Color',[0.65 0.65 0.65],'LineWidth',2)
yline(mean_vx_600,'Color',[0 0 0],'LineWidth',1);
yline(mean_vx_N,'Color',[0.35 0.35 0.35],'LineWidth',1);
yline(-mean_vx_S,'Color',[0.65 0.65 0.65],'LineWidth',1);


legend('Outer Edge (mean: -1.8968 Sv)','Northern Boundary (mean: 0.1621 Sv)',...
    'Southern Boundary (mean: 0.6304)','location','southoutside')

xticks(1:12:204);
xlim([1 204])
xticklabels(years)

ylim([-8 4])
title('Volume Transport Budget where (+) is Offshelf (GLORYS)')
subtitle('-1.2044 Sv')
ylabel('Vol Tx [Sv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_3net_polygon.jpg']);


%%

sum_vx_204 = net_sea_vx+net_N_vx-net_S_vx;

figure
hold on
plot([1:204],sum_vx_204,'-','Color',[0 0 0],'LineWidth',2)

xticks(1:12:204);
xlim([1 204])
xticklabels(years)

ylim([-5 4])
title('NEGS Volume Transport Budget from 81^oN to 74^oN (GLORYS)')
ylabel('Vol Tx [Sv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_3net_net_polygon.jpg']);


%% Plotting where we're taking these values from
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
% cmocean('deep')
m_plot(lon_600,lat_600,'-g','LineWidth',1)
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
m_scatter(lon_nord_span,n_nord_lat,10,'b','Filled');
m_scatter(lon_74_span,s_74_lat,10,'b','Filled');
m_scatter(lon_600,lat_600,10,'b','Filled');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_datapoints.jpg']);

%% "Heatmaps" at the northern and southern edges

% first I need to take the average Volume transport over all the months
vx_nord_mean = nanmean(vx_nord,2);
vx_74_mean = nanmean(vx_74,2);

figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[80 82],'lons',[-15 -5]);
%m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
cmocean('deep')
% m_plot(lon_600,lat_600,'-g','LineWidth',1)
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',80:0.5:82,'xtick',-15:1:-5);
m_scatter(lon_nord_span,n_nord_lat,40,vx_nord_mean*1000,'Filled');
clim([-20 20])
cmocean('balance','pivot',0)
h = colorbar;
ylabel(h,'Volume Transport [mSv]')

title('Average Cross-Shelf Volume Transport in the top 600m: Nordostrundingen (GLORYS)')

% positive indicates offshelf and negative is onshelf
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_avg_nord.jpg']);

%% southern edge
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[73 75],'lons',[-25 -15]);
%m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
cmocean('deep')
% m_plot(lon_600,lat_600,'-g','LineWidth',1)
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',73:0.5:75,'xtick',-25:1:-15);
m_scatter(lon_74_span,s_74_lat,40,-vx_74_mean*1000,'Filled');
clim([-50 50])
cmocean('balance','pivot',0)
h = colorbar;
ylabel(h,'Volume Transport [mSv]')

title('Average Cross-Shelf Volume Transport in the top 600m: 74^oN (GLORYS)')

% positive indicates offshelf and negative is onshelf
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_avg_74.jpg']);

%% eastern edge (600 m contour)

figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
cmocean('deep')
% m_plot(lon_600,lat_600,'-g','LineWidth',1)
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
m_scatter(lon_600,lat_600,40,vx_avg*1000,'Filled');
clim([-100 100])
cmocean('balance','pivot',0)
h = colorbar;
ylabel(h,'Volume Transport [mSv]')

title('Average Cross-Shelf Volume Transport in the top 600m: Eastern Edge (GLORYS)')

% positive indicates offshelf and negative is onshelf
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_avg_outer.jpg']);

%% Finding averags by month
% following the Spall 2023 paper, it may be helpful to note during which
% months we see more or less Volume transport onto/off of the shelf to
% note which processes are 

% finding the indicies of the individual months
jan_i = find(date(:,2) == 1);
feb_i = find(date(:,2) == 2);
mar_i = find(date(:,2) == 3);
apr_i = find(date(:,2) == 4);
may_i = find(date(:,2) == 5);
jun_i = find(date(:,2) == 6);
jul_i = find(date(:,2) == 7);
aug_i = find(date(:,2) == 8);
sep_i = find(date(:,2) == 9);
oct_i = find(date(:,2) == 10);
nov_i = find(date(:,2) == 11);
dec_i = find(date(:,2) == 12);

%% now to apply it to the VX data
% outer 600m isobath
jan_600_vx = nanmean(vx(:,jan_i),2);
feb_600_vx = nanmean(vx(:,feb_i),2);
mar_600_vx = nanmean(vx(:,mar_i),2);
apr_600_vx = nanmean(vx(:,apr_i),2);
may_600_vx = nanmean(vx(:,may_i),2);
jun_600_vx = nanmean(vx(:,jun_i),2);
jul_600_vx = nanmean(vx(:,jul_i),2);
aug_600_vx = nanmean(vx(:,aug_i),2);
sep_600_vx = nanmean(vx(:,sep_i),2);
oct_600_vx = nanmean(vx(:,oct_i),2);
nov_600_vx = nanmean(vx(:,nov_i),2);
dec_600_vx = nanmean(vx(:,dec_i),2);

vx_600_monthly = cat(2,jan_600_vx,feb_600_vx,mar_600_vx,apr_600_vx,may_600_vx,...
    jun_600_vx,jul_600_vx,aug_600_vx,sep_600_vx,oct_600_vx,nov_600_vx,dec_600_vx);


%% heatmap at the monthly resolution

months = {'Jan','Feb','Mar','Apr','May','Ju','Jul','Aug','Sep','Oct','Nov','Dec'};
figure
imagesc([1:12],dist_nord_km,1000*vx_600_monthly)
set(gca,'Ydir','reverse');
clim([-150 150])
colormap(cmocean('balance','pivot',0))

xticks(1:12);
xticklabels(months)
x = colorbar;
ylabel('Distance from Nordostrundingen [km]','FontSize',12)
ylabel(x,'Volume Transport [mSv]','FontSize',12)
title('Monthly Cross-shelf Volume Transport (GLORYS)','FontSize',14)
set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_heat_monthly.jpg']);

%% adding up the individual months and then using a scatterplot
% this may help to see where there's more onshelf and offshelf by month

jan_vx_sum = nansum(jan_600_vx);
feb_vx_sum = nansum(feb_600_vx);
mar_vx_sum = nansum(mar_600_vx);
apr_vx_sum = nansum(apr_600_vx);
may_vx_sum = nansum(may_600_vx);
jun_vx_sum = nansum(jun_600_vx);
jul_vx_sum = nansum(jul_600_vx);
aug_vx_sum = nansum(aug_600_vx);
sep_vx_sum = nansum(sep_600_vx);
oct_vx_sum = nansum(oct_600_vx);
nov_vx_sum = nansum(nov_600_vx);
dec_vx_sum = nansum(dec_600_vx);

vx_month = [jan_vx_sum,feb_vx_sum,mar_vx_sum,apr_vx_sum,may_vx_sum,jun_vx_sum,...
    jul_vx_sum,aug_vx_sum,sep_vx_sum,oct_vx_sum,nov_vx_sum,dec_vx_sum];

figure
scatter([1:12],vx_month,'ko','Filled')
xticks(1:12);
xticklabels(months)
ylim([-3 0])
ylabel('Volume Transport [Sv]')
title('Monthly Volume Transport (GLORYS)')
set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_monthly_scatter.jpg']);

%% NEW FROM RARE_VT %%

%% Plotting the scatterplot on a more seasonal resolution (also as a lineplot)

jan_vx = nansum(vx(:,jan_i),1);
feb_vx = nansum(vx(:,feb_i),1);
mar_vx = nansum(vx(:,mar_i),1);
apr_vx = nansum(vx(:,apr_i),1);
may_vx = nansum(vx(:,may_i),1);
jun_vx = nansum(vx(:,jun_i),1);
jul_vx = nansum(vx(:,jul_i),1);
aug_vx = nansum(vx(:,aug_i),1);
sep_vx = nansum(vx(:,sep_i),1);
oct_vx = nansum(vx(:,oct_i),1);
nov_vx = nansum(vx(:,nov_i),1);
dec_vx = nansum(vx(:,dec_i),1);

vx_year = cat(1,jan_vx,feb_vx,mar_vx,apr_vx,may_vx,jun_vx,...
    jul_vx,aug_vx,sep_vx,oct_vx,nov_vx,dec_vx);


c = cmocean('thermal',17);

figure
hold on
for K = 1:17
    plot([1:12],vx_year(:,K),'Color',c(K,:))
end
plot([1:12],vx_month,'k','LineWidth',2)

legend('2003','2004','2005','2006','2007','2008','2009','2010','2011',...
    '2012','2013','2014','2015','2016','2017','2018','2019','location','eastoutside')

xlim([1 12])
xticks(1:12);
xticklabels(months)
ylim([-10 3])
ylabel('VT [mSv]')
title('Cross-shelf Volume Transport for 2003-2019 (GLORYS)')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_yearly_lineplot.jpg']);

%% Cluster Analysis: Creating a Dendrogram 
% the purpose of this is to try and diagnose how many clusters there are
% within the Volume transport data

tree = linkage(vx,'average');

figure()
dendrogram(tree,25,'ColorThreshold','default')

%% K-means clustering

% 1 cluster
[idx1,C] = kmeans(vx,1);
[silh1,h] = silhouette(vx,idx1);

cluster1 = mean(silh1); % NaN

% 2 cluster
[idx2,C] = kmeans(vx,2);
[silh2,h] = silhouette(vx,idx2);

cluster2 = mean(silh2); 

% 3 cluster
[idx3,C] = kmeans(vx,3);
[silh3,h] = silhouette(vx,idx3);

cluster3 = mean(silh3); 

% 4 cluster
[idx4,C] = kmeans(vx,4);
[silh4,h] = silhouette(vx,idx4);

cluster4 = mean(silh4); 

% 5 cluster
[idx5,C] = kmeans(vx,5);
[silh5,h] = silhouette(vx,idx5);

cluster5 = mean(silh5); 

% 6 cluster
[idx6,C] = kmeans(vx,6);
[silh6,h] = silhouette(vx,idx6);

cluster6 = mean(silh6); 

% 7 cluster
[idx7,C] = kmeans(vx,7);
[silh7,h] = silhouette(vx,idx7);

cluster7 = mean(silh7); 

% 8 cluster
[idx8,C] = kmeans(vx,8);
[silh8,h] = silhouette(vx,idx8);

cluster8 = mean(silh8); 

% 9 cluster
[idx9,C] = kmeans(vx,9);
[silh9,h] = silhouette(vx,idx9);

cluster9 = mean(silh9); 

% 10 cluster
[idx10,C] = kmeans(vx,10);
[silh10,h] = silhouette(vx,idx10);

cluster10 = mean(silh10); 

%%
figure
scatter(idx3,dist_nord_km,'k','Filled');
set(gca,'Ydir','reverse');

title('K-Means Clustering Results (3 clusters, GLORYS)')
xlabel('Cluster Number','FontSize',10)
ylabel('Distance from Nordostrundingen [km]','FontSize',10)

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/3clusters.jpg']);

% it looks like there's two prominant clusters here

%% Plotting the Divergence of the Velocity on the Shelf

% Convert latitude and longitude to radians
latRad = deg2rad(LAT);
lonRad = deg2rad(LON);

% Interpolate the matrix onto the new grid
for i = 1:204
u_4Di(:,:,:,i) = interp3(1:143,1:479,z,u_4D(:,:,:,i),1:143,1:479,z_int);
v_4Di(:,:,:,i) = interp3(1:143,1:479,z,v_4D(:,:,:,i),1:143,1:479,z_int);
end

latRad = deg2rad(lat3D);
lonRad = deg2rad(lon3D);

% Calculating the divergence
for i = 1:204
    for j = 1:51
div_uv(:,:,j,i) = divergence(latRad,lonRad,u_4Di(:,:,j,i),v_4Di(:,:,j,i));
    end
end

%%
% Plotting Divergence
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);

m_contourf(LON,LAT,div_uv(:,:,20,1)',30,'linestyle','none')
cmocean('haline')
colorbar
clim([-30 30])
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');


%% For OSM and beyond, refining certain figures %%

%% Plotting the Velocity and speed plot with the heat scatter

u_3D = nanmean(u_4D,4);
u_600_avg = nanmean(u_3D(:,:,1:21),3);

v_3D = nanmean(v_4D,4);
v_600_avg = nanmean(v_3D(:,:,1:21),3);

speed_600_avg = sqrt(u_600_avg.^2+v_600_avg.^2);

%% Now to plot
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

m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');

m_scatter(lon_nord_span,n_nord_lat,50,vx_nord_mean*1000,'Filled');
cmocean('balance','pivot',0)
clim([-100 100])

m_scatter(lon_74_span,s_74_lat,50,-vx_74_mean*1000,'Filled');
cmocean('balance','pivot',0)
clim([-100 100])

m_scatter(lon_600,lat_600,50,vx_avg*1000,'Filled');
cmocean('balance','pivot',0)
h = colorbar;
clim([-100 100])

ylabel(h,'Volume Transport [mSv]')

m_quiver(LON(1:2:end,1:10:end),LAT(1:2:end,1:10:end),u_600_avg(1:10:end,1:2:end)',...
  v_600_avg(1:10:end,1:2:end)',2.5,'k','LineWidth',0.5,'Clipping','on','AutoScale','on');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
% set(x,'fontsize',10,'fontname','helvetica','position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',0:0.05:0.2);
m_quiver(LON,LAT,u_ref,v_ref,1,'k','LineWidth',0.7,'MaxHeadSize',2);

dim = [0.35 0.36 .045 .06];
str = '0.25 m/s';
annotation('textbox',dim,'String',str,'FontSize',13,'FitBoxToText','off','Rotation',-14.5);
m_plot([-8 -2],[78.9 78.9],'linewi',2,'color','r');

% ylabel(x,'Speed [m/s]')
title('Time Mean Flow and Cross-Shelf Volume Transport: GLORYS 2003-2019','FontSize',16)

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/flow_vx_heat.jpg']);

%% Refining the Volume Budget + Residuals Plot

%% calculating the correlations between different sides of the polygon and
% the residuals
[R_6N,P_6N] = corrcoef(net_sea_vx,net_N_vx);
[R_6S,P_6S] = corrcoef(net_sea_vx,-net_S_vx);
[R_NS,P_NS] = corrcoef(net_N_vx,-net_S_vx);
[R_6r,P_6r] = corrcoef(net_sea_vx,sum_vx_204);
[R_Nr,P_Nr] = corrcoef(net_N_vx,sum_vx_204);
[R_Sr,P_Sr] = corrcoef(-net_S_vx,sum_vx_204);


%%
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9])
hold on
plot([1:204],net_sea_vx,'-','Color',[0 0.5 0],'LineWidth',2)
plot([1:204],net_N_vx,'-','Color',[0 0 0.5],'LineWidth',2)
plot([1:204],-net_S_vx,'-','Color',[0.5 0 0],'LineWidth',2)
plot([1:204],sum_vx_204,'-','Color',[0.5 0.5 0.5],'LineWidth',2)


legend('Outer Edge (mean: -1.9303 Sv)','Northern Boundary (mean: 0.0582 Sv)',...
    'Southern Boundary (mean: 0.5926)','Residual (mean: -1.2795)','location','southoutside')
str = {'R_o_u_t_e_r_,_n_o_r_t_h_e_r_n = 0.3477','R_o_u_t_e_r_,_s_o_u_t_h_e_r_n = -0.8298',...
    'R_n_o_r_t_h_e_r_n_,_s_o_u_t_h_e_r_n = -0.4274','R_o_u_t_e_r_,_r_e_s_i_d_u_a_l = 0.9603',...
    'R_n_o_r_t_h_e_r_n_,_r_e_s_i_d_u_a_l = 0.3575','R_s_o_u_t_h_e_r_n_,_r_e_s_i_d_u_a_l = -0.6512'};
annotation('textbox', [0.8 0.3 0.2 0.05], 'String', str,'FitBoxToText','on',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontSize', 12, 'BackgroundColor', 'none');

% ylim([])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)
ylim([-8 4])
title('Volume Transport Budget and Residual (GLORYS)')
subtitle('(+) is Offshelf')
ylabel('Vol Tx [Sv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/vx_3sides_resids.jpg']);

%% Replotting the speed and velocity plot
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(LON,LAT,speed_600_avg');
cmocean('speed')
x = colorbar;
clim([0 0.2]);

m_contour(lonbath,latbath,bath',[100 350 600 1000],'linew',1.5,'color','w');

m_quiver(LON(1:2:end,1:10:end),LAT(1:2:end,1:10:end),u_600_avg(1:10:end,1:2:end)',...
  v_600_avg(1:10:end,1:2:end)',2.5,'k','LineWidth',0.5,'Clipping','on','AutoScale','on');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
% set(x,'fontsize',10,'fontname','helvetica','position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',0:0.05:0.2);
m_quiver(LON,LAT,u_ref,v_ref,1,'k','LineWidth',0.7,'MaxHeadSize',2);

dim = [0.35 0.36 .045 .06];
str = '0.25 m/s';
annotation('textbox',dim,'String',str,'FontSize',13,'FitBoxToText','off','Rotation',-14.5);
m_plot([-8 -2],[78.9 78.9],'linewi',2,'color','r');

ylabel(x,'Speed [m/s]','FontSize',13)
title('Time Mean Flow and Velocity: GLORYS 2003-2019','FontSize',16)

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/GLORYS_VX_corrected/flow_speed_corrected.jpg']);