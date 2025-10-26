% GLORYS Inital
% Created by: Sara Vianco
% Masters Thesis
% 29 June 2023

%% Preliminaries
close all
clear all

% Set working directory
cd('/Volumes/TERABYTE/GLORYS')
addpath(genpath('/Volumes/TERABYTE/GLORYS/'));
addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/')
addpath('/Users/saravianco/Documents/MATLAB/haversine/')

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

%% Adjusting the time
time_pivot = datenum(1950,1,1,0,0,0);
time_matlab = (time_raw/24) + time_pivot;
time_matlab = floor(time_matlab);
date = datevec(time_matlab);

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

%% Specifying FSAOO location
% condensing the lat and lon for the GLORYS coordinates as well
idxSy = find(abs(lat-latS_map)==nanmin(abs(lat-latS_map)),1,'first');
idxNy = find(abs(lat-latN_map)==nanmin(abs(lat-latN_map)),1,'first');
idxWx = find(abs(lon-lonW_map)==nanmin(abs(lon-lonW_map)),1,'first');
idxEx = find(abs(lon-lonE_map)==nanmin(abs(lon-lonE_map)),1,'first');

% setting the bounds of the bathymetry data to the area that we're
% interested in
lon = lon(idxWx:idxEx);
lat = lat(idxSy:idxNy);

fs_lat_i = find(abs(lat-78.9)==nanmin(abs(lat-78.9)),1,'first');
fs_lon_i = find(lon>=-8 & lon<=-2);

%% Interpolating Z to get a smooth value at 
z_int = NaN(51,1);
z_int(1:35) = z(1:35);
z_int(36) = 1000;
z_int(37:end) = z(36:50);

for i = 1:204
v_4D_interp(:,:,:,i) = interp3(1:143,1:479,z,v_4D(:,:,:,i),1:143,1:479,z_int);
end

%% Separating V into Monthly Means: FSAOO section

% January
jani = find(date(:,2) == 1);
jan_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,jani);
jan_v_mean = squeeze(mean(jan_v_total,4,'omitnan'));

% February
febi = find(date(:,2) == 2);
feb_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,febi);
feb_v_mean = squeeze(mean(feb_v_total,4,'omitnan'));

% March
mari = find(date(:,2) == 3);
mar_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,mari);
mar_v_mean = squeeze(mean(mar_v_total,4,'omitnan'));

% April
apri = find(date(:,2) == 4);
apr_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,apri);
apr_v_mean = squeeze(mean(apr_v_total,4,'omitnan'));

% May
mayi = find(date(:,2) == 5);
may_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,mayi);
may_v_mean = squeeze(mean(may_v_total,4,'omitnan'));

% June
juni = find(date(:,2) == 6);
jun_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,juni);
jun_v_mean = squeeze(mean(jun_v_total,4,'omitnan'));

% July
juli = find(date(:,2) == 7);
jul_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,juli);
jul_v_mean = squeeze(mean(jul_v_total,4,'omitnan'));

% August
augi = find(date(:,2) == 8);
aug_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,augi);
aug_v_mean = squeeze(mean(aug_v_total,4,'omitnan'));

% September
sepi = find(date(:,2) == 9);
sep_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,sepi);
sep_v_mean = squeeze(mean(sep_v_total,4,'omitnan'));

% October
octi = find(date(:,2) == 10);
oct_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,octi);
oct_v_mean = squeeze(mean(oct_v_total,4,'omitnan'));

% November
novi = find(date(:,2) == 11);
nov_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,novi);
nov_v_mean = squeeze(mean(nov_v_total,4,'omitnan'));

% December
deci = find(date(:,2) == 12);
dec_v_total = v_4D_interp(fs_lon_i,fs_lat_i,:,deci);
dec_v_mean = squeeze(mean(dec_v_total,4,'omitnan'));

%% Concatenate Velocity 
v_mean_GLORYS = cat(3,jan_v_mean,feb_v_mean,mar_v_mean,apr_v_mean,may_v_mean,jun_v_mean,...
    jul_v_mean,aug_v_mean,sep_v_mean,oct_v_mean,nov_v_mean,dec_v_mean);

%%
lat_fs = lat(fs_lat_i);
lon_fs = lon(fs_lon_i);

% [Z,LON_FS] = meshgrid(-z,lon(fs_lon_i));
% 
% % Density
% 
% % adding gibbs toolbox to path
% addpath('/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16');
% addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/html/");
% addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/library/");
% addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/thermodynamics_from_t/");
% addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/pdf/");
% 
% % calculating pressure variable from depth variable (z_cell) 
% p = gsw_p_from_z(-z,lat_fs);
% p = p';
% p_mat = repmat(p,length(lon_fs),1);
% 
% %% Practical Salinity Variable
% % January
% jani = find(date(:,2) == 1);
% jan_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,jani);
% jan_s_mean = squeeze(mean(jan_s_total,4,'omitnan'));
% 
% % February
% febi = find(date(:,2) == 1);
% feb_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,febi);
% feb_s_mean = squeeze(mean(feb_s_total,4,'omitnan'));
% 
% % March
% mari = find(date(:,2) == 3);
% mar_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,mari);
% mar_s_mean = squeeze(mean(mar_s_total,4,'omitnan'));
% 
% % April
% apri = find(date(:,2) == 4);
% apr_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,apri);
% apr_s_mean = squeeze(mean(apr_s_total,4,'omitnan'));
% 
% % May
% mayi = find(date(:,2) == 5);
% may_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,mayi);
% may_s_mean = squeeze(mean(may_s_total,4,'omitnan'));
% 
% % June
% juni = find(date(:,2) == 6);
% jun_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,juni);
% jun_s_mean = squeeze(mean(jun_s_total,4,'omitnan'));
% 
% % July
% juli = find(date(:,2) == 7);
% jul_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,juli);
% jul_s_mean = squeeze(mean(jul_s_total,4,'omitnan'));
% 
% % August
% augi = find(date(:,2) == 8);
% aug_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,augi);
% aug_s_mean = squeeze(mean(aug_s_total,4,'omitnan'));
% 
% % September
% sepi = find(date(:,2) == 9);
% sep_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,sepi);
% sep_s_mean = squeeze(mean(sep_s_total,4,'omitnan'));
% 
% % October
% octi = find(date(:,2) == 10);
% oct_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,octi);
% oct_s_mean = squeeze(mean(oct_s_total,4,'omitnan'));
% 
% % November
% novi = find(date(:,2) == 11);
% nov_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,novi);
% nov_s_mean = squeeze(mean(nov_s_total,4,'omitnan'));
% 
% % December
% deci = find(date(:,2) == 12);
% dec_s_total = prac_salinity(fs_lon_i,fs_lat_i,:,deci);
% dec_s_mean = squeeze(mean(dec_s_total,4,'omitnan'));
% 
% %% Potential Temperature Variable
% % January
% jani = find(date(:,2) == 1);
% jan_t_total = tempp(fs_lon_i,fs_lat_i,:,jani);
% jan_t_mean = squeeze(mean(jan_t_total,4,'omitnan'));
% 
% % February
% febi = find(date(:,2) == 1);
% feb_t_total = tempp(fs_lon_i,fs_lat_i,:,febi);
% feb_t_mean = squeeze(mean(feb_t_total,4,'omitnan'));
% 
% % March
% mari = find(date(:,2) == 3);
% mar_t_total = tempp(fs_lon_i,fs_lat_i,:,mari);
% mar_t_mean = squeeze(mean(mar_t_total,4,'omitnan'));
% 
% % April
% apri = find(date(:,2) == 4);
% apr_t_total = tempp(fs_lon_i,fs_lat_i,:,apri);
% apr_t_mean = squeeze(mean(apr_t_total,4,'omitnan'));
% 
% % May
% mayi = find(date(:,2) == 5);
% may_t_total = tempp(fs_lon_i,fs_lat_i,:,mayi);
% may_t_mean = squeeze(mean(may_t_total,4,'omitnan'));
% 
% % June
% juni = find(date(:,2) == 6);
% jun_t_total = tempp(fs_lon_i,fs_lat_i,:,juni);
% jun_t_mean = squeeze(mean(jun_t_total,4,'omitnan'));
% 
% % July
% juli = find(date(:,2) == 7);
% jul_t_total = tempp(fs_lon_i,fs_lat_i,:,juli);
% jul_t_mean = squeeze(mean(jul_t_total,4,'omitnan'));
% 
% % August
% augi = find(date(:,2) == 8);
% aug_t_total = tempp(fs_lon_i,fs_lat_i,:,augi);
% aug_t_mean = squeeze(mean(aug_t_total,4,'omitnan'));
% 
% % September
% sepi = find(date(:,2) == 9);
% sep_t_total = tempp(fs_lon_i,fs_lat_i,:,sepi);
% sep_t_mean = squeeze(mean(sep_t_total,4,'omitnan'));
% 
% % October
% octi = find(date(:,2) == 10);
% oct_t_total = tempp(fs_lon_i,fs_lat_i,:,octi);
% oct_t_mean = squeeze(mean(oct_t_total,4,'omitnan'));
% 
% % November
% novi = find(date(:,2) == 11);
% nov_t_total = tempp(fs_lon_i,fs_lat_i,:,novi);
% nov_t_mean = squeeze(mean(nov_t_total,4,'omitnan'));
% 
% % December
% deci = find(date(:,2) == 12);
% dec_t_total = tempp(fs_lon_i,fs_lat_i,:,deci);
% dec_t_mean = squeeze(mean(dec_t_total,4,'omitnan'));
% 
% %% Concatenating temp and sal
% s_mean = cat(3,jan_s_mean, feb_s_mean, mar_s_mean, apr_s_mean, may_s_mean, jun_s_mean,... 
%     jul_s_mean, aug_s_mean, sep_s_mean, oct_s_mean, nov_s_mean, dec_s_mean);
% 
% t_mean = cat(3,jan_t_mean, feb_t_mean, mar_t_mean, apr_t_mean, may_t_mean, jun_t_mean,... 
%     jul_t_mean, aug_t_mean, sep_t_mean, oct_t_mean, nov_t_mean, dec_t_mean);
% 
% %% Individual variables for density calculation
% 
% SA = NaN(size(s_mean));
% for i = 1:12
% [SA(:,:,i), in_ocean] = gsw_SA_from_SP(s_mean(:,:,i),p_mat,lon_fs,lat_fs);
% end
% 
% CT = NaN(size(t_mean));
% for i = 1:12
% CT(:,:,i) = gsw_CT_from_pt(SA(:,:,i),t_mean(:,:,i));
% end
% 
% %% Density Calculation
% 
% rho_CT_exact = NaN(size(t_mean));
% for i = 1:12
% rho_CT_exact(:,:,i) = gsw_rho_CT_exact(SA(:,:,i),CT(:,:,i),p_mat);
% end
% 
% rho_anom = rho_CT_exact-1000;
% 
% %% Velocity and Density Imaging
% % January
% figure
% hold on
% contourf(LON_FS,Z,jan_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,1),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(jan_v_mean))=0;
% imagesc(jan_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean January EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/01.jpg']);
% 
% % February
% figure
% hold on
% contourf(LON_FS,Z,feb_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,2),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(feb_v_mean));
% imAlpha(isnan(feb_v_mean))=0;
% imagesc(feb_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean February EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/02.jpg']);
% 
% % March
% figure
% hold on
% contourf(LON_FS,Z,mar_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,3),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(mar_v_mean));
% imAlpha(isnan(mar_v_mean))=0;
% imagesc(mar_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean March EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/03.jpg']);
% 
% % April
% figure
% hold on
% contourf(LON_FS,Z,apr_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,4),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(apr_v_mean));
% imAlpha(isnan(apr_v_mean))=0;
% imagesc(apr_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean April EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/04.jpg']);
% 
% % May
% figure
% hold on
% contourf(LON_FS,Z,may_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,5),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(may_v_mean));
% imAlpha(isnan(may_v_mean))=0;
% imagesc(may_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean May EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/05.jpg']);
% 
% % June
% figure
% hold on
% contourf(LON_FS,Z,jun_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,6),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(jun_v_mean));
% imAlpha(isnan(jun_v_mean))=0;
% imagesc(jun_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean June EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/06.jpg']);
% 
% % July
% figure
% hold on
% contourf(LON_FS,Z,jul_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,7),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(jul_v_mean));
% imAlpha(isnan(jul_v_mean))=0;
% imagesc(jul_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean July EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/07.jpg']);
% 
% % August
% figure
% hold on
% contourf(LON_FS,Z,aug_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,8),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(aug_v_mean));
% imAlpha(isnan(aug_v_mean))=0;
% imagesc(aug_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean August EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/08.jpg']);
% 
% % September
% figure
% hold on
% contourf(LON_FS,Z,sep_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,9),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(sep_v_mean));
% imAlpha(isnan(sep_v_mean))=0;
% imagesc(sep_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean September EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/09.jpg']);
% 
% % October
% figure
% hold on
% contourf(LON_FS,Z,oct_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,10),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(oct_v_mean));
% imAlpha(isnan(oct_v_mean))=0;
% imagesc(oct_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean October EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/10.jpg']);
% 
% % November
% figure
% hold on
% contourf(LON_FS,Z,nov_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,11),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(nov_v_mean));
% imAlpha(isnan(nov_v_mean))=0;
% imagesc(nov_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean November EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/11.jpg']);
% 
% % December
% figure
% hold on
% contourf(LON_FS,Z,dec_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(LON_FS,Z,rho_anom(:,:,12),[25:0.5:32.5],'-k');
% 
% imAlpha=ones(size(dec_v_mean));
% imAlpha(isnan(dec_v_mean))=0;
% imagesc(dec_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('Mean December EGS Velocity (GLORYS)','FontSize',16)
% xlabel('Longitude [ ^oE]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'GLORYS/monthly_v_fs/12.jpg']);

%% Volumetric Transport
lon_1 = -6;
lon_2 = -2;

lat_fs = 78.9;

lon1_idx = find(lon_fs==lon_1);
lon2_idx = find(lon_fs==lon_2);
gridx_diff = lon2_idx-lon1_idx;

lon_km = haversine([lat_fs,lon_1],[lat_fs,lon_2]);
lon_m = lon_km*1000;
dx = lon_m/gridx_diff; % horizontal distance of a gridcell

% I need to find the distances between the z cells at the midpoints
dz = NaN(size(z_int))';

for kk = 2:length(z_int)-1
dz(kk) = ((z_int(kk+1)+z_int(kk))/2)-((z_int(kk)+z_int(kk-1))/2);
end
dz(1) = ((z_int(2)+z_int(1))/2)-(z_int(1)/2);

dz = dz(1:36)';

%% Monthly VX

% January
jan_v_mat = squeeze(jan_v_total);
jan_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_jan_mat = NaN(1,17);

for i = 1:17
jan_egc_mat(:,:,i) = jan_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(jan_egc_mat(:,:,i)>0);
jan_egc_mat(row,col,i) = 0;
jan_egc_mat(isnan(jan_egc_mat)) = 0;
tx_jan_mat(1,i) = (sum(dx*jan_egc_mat(:,:,i)*dz))/10^6;
end

% February
feb_v_mat = squeeze(feb_v_total);
feb_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_feb_mat = NaN(1,17);

for i = 1:17
feb_egc_mat(:,:,i) = feb_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(feb_egc_mat(:,:,i)>0);
feb_egc_mat(row,col,i) = 0;
feb_egc_mat(isnan(feb_egc_mat)) = 0;
tx_feb_mat(1,i) = (sum(dx*feb_egc_mat(:,:,i)*dz))/10^6;
end

% March
mar_v_mat = squeeze(mar_v_total);
mar_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_mar_mat = NaN(1,17);

for i = 1:17
mar_egc_mat(:,:,i) = mar_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(mar_egc_mat(:,:,i)>0);
mar_egc_mat(row,col,i) = 0;
mar_egc_mat(isnan(mar_egc_mat)) = 0;
tx_mar_mat(1,i) = (sum(dx*mar_egc_mat(:,:,i)*dz))/10^6;
end

% April
apr_v_mat = squeeze(apr_v_total);
apr_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_apr_mat = NaN(1,17);

for i = 1:17
apr_egc_mat(:,:,i) = apr_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(apr_egc_mat(:,:,i)>0);
apr_egc_mat(row,col,i) = 0;
apr_egc_mat(isnan(apr_egc_mat)) = 0;
tx_apr_mat(1,i) = (sum(dx*apr_egc_mat(:,:,i)*dz))/10^6;
end

% May
may_v_mat = squeeze(may_v_total);
may_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_may_mat = NaN(1,17);

for i = 1:17
may_egc_mat(:,:,i) = may_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(may_egc_mat(:,:,i)>0);
may_egc_mat(row,col,i) = 0;
may_egc_mat(isnan(may_egc_mat)) = 0;
tx_may_mat(1,i) = (sum(dx*may_egc_mat(:,:,i)*dz))/10^6;
end

% June
jun_v_mat = squeeze(jun_v_total);
jun_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_jun_mat = NaN(1,17);

for i = 1:17
jun_egc_mat(:,:,i) = jun_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(jun_egc_mat(:,:,i)>0);
jun_egc_mat(row,col,i) = 0;
jun_egc_mat(isnan(jun_egc_mat)) = 0;
tx_jun_mat(1,i) = (sum(dx*jun_egc_mat(:,:,i)*dz))/10^6;
end

% July
jul_v_mat = squeeze(jul_v_total);
jul_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_jul_mat = NaN(1,17);

for i = 1:17
jul_egc_mat(:,:,i) = jul_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(jul_egc_mat(:,:,i)>0);
jul_egc_mat(row,col,i) = 0;
jul_egc_mat(isnan(jul_egc_mat)) = 0;
tx_jul_mat(1,i) = (sum(dx*jul_egc_mat(:,:,i)*dz))/10^6;
end

% August
aug_v_mat = squeeze(aug_v_total);
aug_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_aug_mat = NaN(1,17);

for i = 1:15
aug_egc_mat(:,:,i) = aug_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(aug_egc_mat(:,:,i)>0);
aug_egc_mat(row,col,i) = 0;
aug_egc_mat(isnan(aug_egc_mat)) = 0;
tx_aug_mat(1,i) = (sum(dx*aug_egc_mat(:,:,i)*dz))/10^6;
end

% September
sep_v_mat = squeeze(sep_v_total);
sep_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_sep_mat = NaN(1,17);

for i = 1:17
sep_egc_mat(:,:,i) = sep_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(sep_egc_mat(:,:,i)>0);
sep_egc_mat(row,col,i) = 0;
sep_egc_mat(isnan(sep_egc_mat)) = 0;
tx_sep_mat(1,i) = (sum(dx*sep_egc_mat(:,:,i)*dz))/10^6;
end

% October
oct_v_mat = squeeze(oct_v_total);
oct_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_oct_mat = NaN(1,17);

for i = 1:17
oct_egc_mat(:,:,i) = oct_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(oct_egc_mat(:,:,i)>0);
oct_egc_mat(row,col,i) = 0;
oct_egc_mat(isnan(oct_egc_mat)) = 0;
tx_oct_mat(1,i) = (sum(dx*oct_egc_mat(:,:,i)*dz))/10^6;
end

% November
nov_v_mat = squeeze(nov_v_total);
nov_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_nov_mat = NaN(1,17);

for i = 1:17
nov_egc_mat(:,:,i) = nov_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(nov_egc_mat(:,:,i)>0);
nov_egc_mat(row,col,i) = 0;
nov_egc_mat(isnan(nov_egc_mat)) = 0;
tx_nov_mat(1,i) = (sum(dx*nov_egc_mat(:,:,i)*dz))/10^6;
end

% December
dec_v_mat = squeeze(dec_v_total);
dec_egc_mat = NaN(length(lon_fs(lon1_idx:lon2_idx)),36,16);
tx_dec_mat = NaN(1,17);

for i = 1:17
dec_egc_mat(:,:,i) = dec_v_mat(lon1_idx:lon2_idx,1:36,i);
[row,col] = find(dec_egc_mat(:,:,i)>0);
dec_egc_mat(row,col,i) = 0;
dec_egc_mat(isnan(dec_egc_mat)) = 0;
tx_dec_mat(1,i) = (sum(dx*dec_egc_mat(:,:,i)*dz))/10^6;
end

%%
full_tx_GLORYS = cat(1,tx_jan_mat, tx_feb_mat, tx_mar_mat, tx_apr_mat, tx_may_mat, ...
    tx_jun_mat, tx_jul_mat, tx_aug_mat, tx_sep_mat, tx_oct_mat, tx_nov_mat, ...
    tx_dec_mat);


month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
c = cmocean('thermal',17);

%% Determining Errorbars/SD/95% Confidence bars (1.96*SD)
sd_vx_GLORYS = NaN(1,12);
mean_vx_GLORYS = NaN(1,12);

for jj = 1:12
sd_vx_GLORYS(1,jj) = std(full_tx_GLORYS(jj,:),'omitnan');
mean_vx_GLORYS(1,jj) = mean(full_tx_GLORYS(jj,:),'omitnan');
end

%% Figures

% Volumetric Transport for each year, at each month. The black line
% represents the average tx for that month
figure
hold on
for i = 1:17
plot(full_tx_GLORYS(:,i),'linewi',1,'Color',c(i,:))
end
plot(mean_vx_GLORYS,'-k','linewi',3)
xlim([1 12])
% ylim([-1 0.2])
xticks(1:12);
xticklabels(month)

legend('2003','2004','2005','2006','2007','2008','2009','2010','2011','2012' ...
    ,'2013','2014','2015','2016','2017','2018','2019','location','east outside')

ylabel('Transport [Sv]')
title('EGC Monthly Volumetric Transport 2003-2019 (GLORYS)')
saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/GLORYS/vol_tx/vx_all.jpg');

