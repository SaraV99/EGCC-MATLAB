% RARE full shelf--Redone and reduced to a yearly average
% Created by: Sara Vianco
% Masters Thesis
% 23 January 2024

%% Preliminaries
close all
clear all

cd('/Volumes/TERABYTE/RARE/matching_theo/')
addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/')

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

%% adjusting the time

% making the date number into one that is readable by MATLAB, forming a
% date string afterwards

base = datenum(1979,12,31,0,0,0);
formatIn = 'mm/dd/yy';
   
date_num = datenum(time + base);
date_string = datestr(date_num,'mm/dd/yy');
dv = datevec(date_string,formatIn);


%% now to reduce the varibles down to ths vertical sections

lat_7395_i = find(abs(lat3D-73.95)==nanmin(abs(lat3D-73.95)),1,'first');
lat_7695_i = find(abs(lat3D-76.95)==nanmin(abs(lat3D-76.95)),1,'first');
lat_7885_i = find(abs(lat3D-78.85)==nanmin(abs(lat3D-78.85)),1,'first');
lat_8005_i = find(abs(lat3D-80.05)==nanmin(abs(lat3D-80.05)),1,'first');
lat_8145_i = find(abs(lat3D-81.45)==nanmin(abs(lat3D-81.45)),1,'first');

lat_vert_sections = [lat3D(lat_7395_i),lat3D(lat_7695_i),lat3D(lat_7885_i),...
    lat3D(lat_8005_i),lat3D(lat_8145_i)];

%% Creating a time average of V, Psal, and Ptemp
v_3D = nanmean(v_4D,4);
psal_3D = nanmean(psal_4D,4);
ptemp_3D = nanmean(ptemp_4D,4);

%% Establishing Density Contours: Full Shelf %%

% adding gibbs toolbox to path
addpath('/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16');
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/html/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/library/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/thermodynamics_from_t/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/pdf/");

% calculating pressure variable from depth variable (z_cell)

for i = 1:length(lat_vert_sections)
p(:,i) = gsw_p_from_z(-z_star,lat_vert_sections(i));
end

% p = p';
p_mat_7395 = repmat(p(:,1),1,length(lon3D))';
p_mat_7695 = repmat(p(:,2),1,length(lon3D))';
p_mat_7885 = repmat(p(:,3),1,length(lon3D))';
p_mat_8005 = repmat(p(:,4),1,length(lon3D))';
p_mat_8145 = repmat(p(:,5),1,length(lon3D))';

%% Individual variables for density calculation

% @ 73.95N
[SA_7395, in_ocean] = gsw_SA_from_SP(squeeze(psal_3D(:,lat_7395_i,:)),p_mat_7395,lon3D,lat3D(lat_7395_i));
CT_7395 = gsw_CT_from_pt(SA_7395,squeeze(ptemp_3D(:,lat_7395_i,:)));

rho_7395 = gsw_rho_CT_exact(SA_7395,CT_7395,p_mat_7395);
alpha_7395 = rho_7395-1000;

% @ 76.95N
[SA_7695, in_ocean] = gsw_SA_from_SP(squeeze(psal_3D(:,lat_7695_i,:)),p_mat_7695,lon3D,lat3D(lat_7695_i));
CT_7695 = gsw_CT_from_pt(SA_7695,squeeze(ptemp_3D(:,lat_7695_i,:)));

rho_7695 = gsw_rho_CT_exact(SA_7695,CT_7695,p_mat_7695);
alpha_7695 = rho_7695-1000;

% @ 78.85N
[SA_7885, in_ocean] = gsw_SA_from_SP(squeeze(psal_3D(:,lat_7885_i,:)),p_mat_7885,lon3D,lat3D(lat_7885_i));
CT_7885 = gsw_CT_from_pt(SA_7885,squeeze(ptemp_3D(:,lat_7885_i,:)));

rho_7885 = gsw_rho_CT_exact(SA_7885,CT_7885,p_mat_7885);
alpha_7885 = rho_7885-1000;

% @ 80.05N
[SA_8005, in_ocean] = gsw_SA_from_SP(squeeze(psal_3D(:,lat_8005_i,:)),p_mat_8005,lon3D,lat3D(lat_8005_i));
CT_8005 = gsw_CT_from_pt(SA_8005,squeeze(ptemp_3D(:,lat_8005_i,:)));

rho_8005 = gsw_rho_CT_exact(SA_8005,CT_8005,p_mat_8005);
alpha_8005 = rho_8005-1000;

% @ 81.45N
[SA_8145, in_ocean] = gsw_SA_from_SP(squeeze(psal_3D(:,lat_8145_i,:)),p_mat_8145,lon3D,lat3D(lat_8145_i));
CT_8145 = gsw_CT_from_pt(SA_8145,squeeze(ptemp_3D(:,lat_8145_i,:)));

rho_8145 = gsw_rho_CT_exact(SA_8145,CT_8145,p_mat_8145);
alpha_8145 = rho_8145-1000;


%% Now to plot: these are overall averages

%% At 73.95N
shelf_lon_i = find(lon3D>-20 & lon3D<-8);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

v_7395 = squeeze(v_3D(shelf_lon_i,lat_7395_i,1:24));

figure
hold on
contourf(LON,Z,v_7395,8,'linestyle','none');
clim([-0.2 0.2])
colormap(cmocean('balance','pivot',0))
x = colorbar;
[C,d] = contour(LON,Z,alpha_7395(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_7395));
imAlpha(isnan(v_7395))=0;
imagesc(v_7395,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-20 -8])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'v [m/s]')
title('Mean NEGS V-Velocity and α at 73.95^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/7395.jpg']);

%% At 76.95
shelf_lon_i = find(lon3D>-19 & lon3D<0);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

v_7695 = squeeze(v_3D(shelf_lon_i,lat_7695_i,1:24));

figure
hold on
contourf(LON,Z,v_7695,8,'linestyle','none');
clim([-0.2 0.2])
colormap(cmocean('balance','pivot',0))
x = colorbar;
[C,d] = contour(LON,Z,alpha_7695(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_7695));
imAlpha(isnan(v_7695))=0;
imagesc(v_7695,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-19 -0])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'v [m/s]')
title('Mean NEGS V-Velocity and α at 76.95^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/7695.jpg']);

%% At 78.85N
shelf_lon_i = find(lon3D>-18 & lon3D<0);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

v_7885 = squeeze(v_3D(shelf_lon_i,lat_7885_i,1:24));

figure
hold on
contourf(LON,Z,v_7885,8,'linestyle','none');
clim([-0.2 0.2])
colormap(cmocean('balance','pivot',0))
x = colorbar;
[C,d] = contour(LON,Z,alpha_7885(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_7885));
imAlpha(isnan(v_7885))=0;
imagesc(v_7885,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-18 0])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'v [m/s]')
title('Mean NEGS V-Velocity and α at 78.85^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/7885.jpg']);

%% At 80.05N
shelf_lon_i = find(lon3D>-18 & lon3D<0);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

v_8005 = squeeze(v_3D(shelf_lon_i,lat_8005_i,1:24));

figure
hold on
contourf(LON,Z,v_8005,8,'linestyle','none');
clim([-0.2 0.2])
colormap(cmocean('balance','pivot',0))
x = colorbar;
[C,d] = contour(LON,Z,alpha_8005(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_8005));
imAlpha(isnan(v_8005))=0;
imagesc(v_8005,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-18 0])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'v [m/s]')
title('Mean NEGS V-Velocity and α at 80.05^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/8005.jpg']);

%% At 81.45N
shelf_lon_i = find(lon3D>-13 & lon3D<0);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

v_8145 = squeeze(v_3D(shelf_lon_i,lat_8145_i,1:24));

figure
hold on
contourf(LON,Z,v_8145,8,'linestyle','none');
clim([-0.2 0.2])
colormap(cmocean('balance','pivot',0))
x = colorbar;
[C,d] = contour(LON,Z,alpha_8145(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_8145));
imAlpha(isnan(v_8145))=0;
imagesc(v_8145,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-13 0])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'v [m/s]')
title('Mean NEGS V-Velocity and α at 81.45^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/8145.jpg']);

%% Now to plot the vertical sections of salinity and density

%% At 73.95N
shelf_lon_i = find(lon3D>-20 & lon3D<-8);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

figure
hold on
contourf(LON,Z,SA_7395(shelf_lon_i,1:24),12,'linestyle','none');

clim([31 36])
colormap(cmocean('haline'))
x = colorbar;
[C,d] = contour(LON,Z,alpha_7395(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_7395));
imAlpha(isnan(v_7395))=0;
imagesc(v_7395,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-20 -8])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'Absolute Salinity [PSU]')
title('Mean NEGS Salinity and α at 73.95^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/salinity:density/7395.jpg']);

%% At 76.95
shelf_lon_i = find(lon3D>-19 & lon3D<0);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

figure
hold on
contourf(LON,Z,SA_7695(shelf_lon_i,1:24),12,'linestyle','none');

clim([31 36])
colormap(cmocean('haline'))
x = colorbar;
[C,d] = contour(LON,Z,alpha_7695(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_7695));
imAlpha(isnan(v_7695))=0;
imagesc(v_7695,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-19 -0])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'Absolute Salinity [PSU]')
title('Mean NEGS Salinity and α at 76.95^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/salinity:density/7695.jpg']);

%% At 78.85
shelf_lon_i = find(lon3D>-18 & lon3D<0);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

figure
hold on
contourf(LON,Z,SA_7885(shelf_lon_i,1:24),12,'linestyle','none');

clim([31 36])
colormap(cmocean('haline'))
x = colorbar;
[C,d] = contour(LON,Z,alpha_7885(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_7885));
imAlpha(isnan(v_7885))=0;
imagesc(v_7885,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-18 -0])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'Absolute Salinity [PSU]')
title('Mean NEGS Salinity and α at 78.85^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/salinity:density/7885.jpg']);

%% At 80.05
shelf_lon_i = find(lon3D>-18 & lon3D<0);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

figure
hold on
contourf(LON,Z,SA_8005(shelf_lon_i,1:24),12,'linestyle','none');

clim([31 36])
colormap(cmocean('haline'))
x = colorbar;
[C,d] = contour(LON,Z,alpha_8005(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_8005));
imAlpha(isnan(v_8005))=0;
imagesc(v_8005,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-18 -0])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'Absolute Salinity [PSU]')
title('Mean NEGS Salinity and α at 80.05^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/salinity:density/8005.jpg']);

%% At 81.45
shelf_lon_i = find(lon3D>-13 & lon3D<0);
[Z,LON] = meshgrid(-z_star(1:24),lon3D(shelf_lon_i));

figure
hold on
contourf(LON,Z,SA_8145(shelf_lon_i,1:24),12,'linestyle','none');

clim([31 36])
colormap(cmocean('haline'))
x = colorbar;
[C,d] = contour(LON,Z,alpha_8145(shelf_lon_i,1:24),[25:0.5:32.5],'-k');
clabel(C,d);

imAlpha=ones(size(v_8145));
imAlpha(isnan(v_8145))=0;
imagesc(v_8145,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);
ylim([-1000 0])
xlim([-13 -0])

xlabel('Longitude [^oE]')
ylabel('Depth [m]')
ylabel(x,'Absolute Salinity [PSU]')
title('Mean NEGS Salinity and α at 81.45^oN');

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/salinity:density/8145.jpg']);