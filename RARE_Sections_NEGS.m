% RARE full shelf--Redone and reduced to a yearly average
% Created by: Sara Vianco
% Masters Thesis
% 23 January 2024

%% Preliminaries
close all
clear all

cd('/Volumes/TERABYTE/RARE/matching_theo/')
addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/')
addpath("/Users/saravianco/Documents/MATLAB/raacampbell-shadedErrorBar-aa6d919/"); % this is for the shaded error bar

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
latS_map = 70;
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

lat_7095_i = find(abs(lat3D-70.95)==nanmin(abs(lat3D-70.95)),1,'first');
lat_7395_i = find(abs(lat3D-73.95)==nanmin(abs(lat3D-73.95)),1,'first');
lat_7695_i = find(abs(lat3D-76.95)==nanmin(abs(lat3D-76.95)),1,'first');
lat_7885_i = find(abs(lat3D-78.85)==nanmin(abs(lat3D-78.85)),1,'first');
lat_8005_i = find(abs(lat3D-80.05)==nanmin(abs(lat3D-80.05)),1,'first');
lat_8145_i = find(abs(lat3D-81.45)==nanmin(abs(lat3D-81.45)),1,'first');

lat_vert_sections = [lat3D(lat_7095_i),lat3D(lat_7395_i),lat3D(lat_7695_i),lat3D(lat_7885_i),...
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
p_mat_7095 = repmat(p(:,1),1,length(lon3D))';
p_mat_7395 = repmat(p(:,1),1,length(lon3D))';
p_mat_7695 = repmat(p(:,2),1,length(lon3D))';
p_mat_7885 = repmat(p(:,3),1,length(lon3D))';
p_mat_8005 = repmat(p(:,4),1,length(lon3D))';
p_mat_8145 = repmat(p(:,5),1,length(lon3D))';

%% Individual variables for density calculation

% @ 70.95N
clear SA_7095
clear CT_7095

for i = 1:204
[SA_7095(:,:,:,i), in_ocean] = gsw_SA_from_SP(squeeze(psal_4D(:,lat_7095_i,:,i)),p_mat_7095,lon3D,lat3D(lat_7095_i));
CT_7095(:,:,:,i) = gsw_CT_from_pt(SA_7095(:,:,:,i),squeeze(ptemp_4D(:,lat_7095_i,:,i)));
end 

SA_7095 = squeeze(SA_7095);
CT_7095 = squeeze(CT_7095);

for i = 1:204
rho_7095(:,:,i) = gsw_rho_CT_exact(SA_7095(:,:,i),CT_7095(:,:,i),p_mat_7095);
alpha_7095 = rho_7095-1000;
end

% @ 73.95N
clear SA_7395
clear CT_7395

for i = 1:204
[SA_7395(:,:,:,i), in_ocean] = gsw_SA_from_SP(squeeze(psal_4D(:,lat_7395_i,:,i)),p_mat_7395,lon3D,lat3D(lat_7395_i));
CT_7395(:,:,:,i) = gsw_CT_from_pt(SA_7395(:,:,:,i),squeeze(ptemp_4D(:,lat_7395_i,:,i)));
end 

SA_7395 = squeeze(SA_7395);
CT_7395 = squeeze(CT_7395);

for i = 1:204
rho_7395(:,:,i) = gsw_rho_CT_exact(SA_7395(:,:,i),CT_7395(:,:,i),p_mat_7395);
alpha_7395 = rho_7395-1000;
end

% @ 76.95N
clear SA_7695
clear CT_7695

for i = 1:204
[SA_7695(:,:,:,i), in_ocean] = gsw_SA_from_SP(squeeze(psal_4D(:,lat_7695_i,:,i)),p_mat_7695,lon3D,lat3D(lat_7695_i));
CT_7695(:,:,:,i) = gsw_CT_from_pt(SA_7695(:,:,:,i),squeeze(ptemp_4D(:,lat_7695_i,:,i)));
end 

SA_7695 = squeeze(SA_7695);
CT_7695 = squeeze(CT_7695);

for i = 1:204
rho_7695(:,:,i) = gsw_rho_CT_exact(SA_7695(:,:,i),CT_7695(:,:,i),p_mat_7695);
alpha_7695 = rho_7695-1000;
end


% @ 78.85N
clear SA_7885
clear CT_7885

for i = 1:204
[SA_7885(:,:,:,i), in_ocean] = gsw_SA_from_SP(squeeze(psal_4D(:,lat_7885_i,:,i)),p_mat_7885,lon3D,lat3D(lat_7885_i));
CT_7885(:,:,:,i) = gsw_CT_from_pt(SA_7885(:,:,:,i),squeeze(ptemp_4D(:,lat_7885_i,:,i)));
end 

SA_7885 = squeeze(SA_7885);
CT_7885 = squeeze(CT_7885);

for i = 1:204
rho_7885(:,:,i) = gsw_rho_CT_exact(SA_7885(:,:,i),CT_7885(:,:,i),p_mat_7885);
alpha_7885 = rho_7885-1000;
end

% @ 80.05N
clear SA_8005
clear CT_8005

for i = 1:204
[SA_8005(:,:,:,i), in_ocean] = gsw_SA_from_SP(squeeze(psal_4D(:,lat_8005_i,:,i)),p_mat_8005,lon3D,lat3D(lat_8005_i));
CT_8005(:,:,:,i) = gsw_CT_from_pt(SA_8005(:,:,:,i),squeeze(ptemp_4D(:,lat_8005_i,:,i)));
end 

SA_8005 = squeeze(SA_8005);
CT_8005 = squeeze(CT_8005);

for i = 1:204
rho_8005(:,:,i) = gsw_rho_CT_exact(SA_8005(:,:,i),CT_8005(:,:,i),p_mat_8005);
alpha_8005 = rho_8005-1000;
end


% @ 81.45N
clear SA_8145
clear CT_8145

for i = 1:204
[SA_8145(:,:,:,i), in_ocean] = gsw_SA_from_SP(squeeze(psal_4D(:,lat_8145_i,:,i)),p_mat_8145,lon3D,lat3D(lat_8145_i));
CT_8145(:,:,:,i) = gsw_CT_from_pt(SA_8145(:,:,:,i),squeeze(ptemp_4D(:,lat_8145_i,:,i)));
end 

SA_8145 = squeeze(SA_8145);
CT_8145 = squeeze(CT_8145);

for i = 1:204
rho_8145(:,:,i) = gsw_rho_CT_exact(SA_8145(:,:,i),CT_8145(:,:,i),p_mat_8145);
alpha_8145 = rho_8145-1000;
end

%% Identifying the 34 isohaline %%
% this will delinate our Fresh Water that the EGC and EGCC are transporting

%% @ 70.95N
yn_34_7095 = NaN(size(SA_7095));
for i = 1:204
  for j = 1:38
      for k = 1:400
    if isnan(SA_7095(k,j,i)) || SA_7095(k,j,i)>34
        yn_34_7095(k,j,i) = 0;
    else
        yn_34_7095(k,j,i) = 1;
      end
    end
  end
end

%% @ 73.95N
yn_34_7395 = NaN(size(SA_7395));
for i = 1:204
  for j = 1:38
      for k = 1:400
    if isnan(SA_7395(k,j,i)) || SA_7395(k,j,i)>34
        yn_34_7395(k,j,i) = 0;
    else
        yn_34_7395(k,j,i) = 1;
      end
    end
  end
end

%% @ 76.95N
yn_34_7695 = NaN(size(SA_7695));
for i = 1:204
  for j = 1:38
      for k = 1:400
    if isnan(SA_7695(k,j,i)) || SA_7695(k,j,i)>34
        yn_34_7695(k,j,i) = 0;
    else
        yn_34_7695(k,j,i) = 1;
      end
    end
  end
end

%% @ 78.85N
yn_34_7885 = NaN(size(SA_7885));
for i = 1:204
  for j = 1:38
      for k = 1:400
    if isnan(SA_7885(k,j,i)) || SA_7885(k,j,i)>34
        yn_34_7885(k,j,i) = 0;
    else
        yn_34_7885(k,j,i) = 1;
      end
    end
  end
end
%% @ 80.05N
yn_34_8005 = NaN(size(SA_8005));
for i = 1:204
  for j = 1:38
      for k = 1:400
    if isnan(SA_8005(k,j,i)) || SA_8005(k,j,i)>34
        yn_34_8005(k,j,i) = 0;
    else
        yn_34_8005(k,j,i) = 1;
      end
    end
  end
end
%% @ 81.45N
yn_34_8145 = NaN(size(SA_8145));
for i = 1:204
  for j = 1:38
      for k = 1:400
    if isnan(SA_8145(k,j,i)) || SA_8145(k,j,i)>34
        yn_34_8145(k,j,i) = 0;
    else
        yn_34_8145(k,j,i) = 1;
      end
    end
  end
end

%% Now to calculate dz and dx for the volume transport
addpath('/Users/saravianco/Documents/MATLAB//haversine/') % this is for haversine

% @ 70.95
lon_km_7095 = NaN(length(lon3D)-1,1);
for jj = 1:(length(lon3D)-1)
lon_km_7095(jj,1) = haversine([lat3D(lat_7095_i),lon3D(jj)],[lat3D(lat_7095_i),lon3D(jj+1)]); 
% haversign calculates the geometric distance in km between two coordinates
end
lon_m_7095 = lon_km_7095*1000;
lon_m_7095(end+1,1) = lon_m_7095(end,1);

% now to calculate the half-way distance for the Volume transports
lon_diff_7095 = NaN(size(lon_m_7095))';

for kk = 1:length(lon_m_7095)-1
lon_diff_7095(kk) = (lon_m_7095(kk+1)/2)+(lon_m_7095(kk)/2);
end
lon_diff_7095(end) = lon_m_7095(end)/2;

% @ 73.95
lon_km_7395 = NaN(length(lon3D)-1,1);
for jj = 1:(length(lon3D)-1)
lon_km_7395(jj,1) = haversine([lat3D(lat_7395_i),lon3D(jj)],[lat3D(lat_7395_i),lon3D(jj+1)]); 
% haversign calculates the geometric distance in km between two coordinates
end
lon_m_7395 = lon_km_7395*1000;
lon_m_7395(end+1,1) = lon_m_7395(end,1);

% now to calculate the half-way distance for the Volume transports
lon_diff_7395 = NaN(size(lon_m_7395))';

for kk = 1:length(lon_m_7395)-1
lon_diff_7395(kk) = (lon_m_7395(kk+1)/2)+(lon_m_7395(kk)/2);
end
lon_diff_7395(end) = lon_m_7395(end)/2;

% @ 76.95
lon_km_7695 = NaN(length(lon3D)-1,1);
for jj = 1:(length(lon3D)-1)
lon_km_7695(jj,1) = haversine([lat3D(lat_7695_i),lon3D(jj)],[lat3D(lat_7695_i),lon3D(jj+1)]); 
% haversign calculates the geometric distance in km between two coordinates
end
lon_m_7695 = lon_km_7695*1000;
lon_m_7695(end+1,1) = lon_m_7695(end,1);

% now to calculate the half-way distance for the Volume transports
lon_diff_7695 = NaN(size(lon_m_7695))';

for kk = 1:length(lon_m_7695)-1
lon_diff_7695(kk) = (lon_m_7695(kk+1)/2)+(lon_m_7695(kk)/2);
end
lon_diff_7695(end) = lon_m_7695(end)/2;

% @ 78.85
lon_km_7885 = NaN(length(lon3D)-1,1);
for jj = 1:(length(lon3D)-1)
lon_km_7885(jj,1) = haversine([lat3D(lat_7885_i),lon3D(jj)],[lat3D(lat_7885_i),lon3D(jj+1)]); 
% haversign calculates the geometric distance in km between two coordinates
end
lon_m_7885 = lon_km_7885*1000;
lon_m_7885(end+1,1) = lon_m_7885(end,1);

% now to calculate the half-way distance for the Volume transports
lon_diff_7885 = NaN(size(lon_m_7885))';

for kk = 1:length(lon_m_7885)-1
lon_diff_7885(kk) = (lon_m_7885(kk+1)/2)+(lon_m_7885(kk)/2);
end
lon_diff_7885(end) = lon_m_7885(end)/2;

% @ 80.05
lon_km_8005 = NaN(length(lon3D)-1,1);
for jj = 1:(length(lon3D)-1)
lon_km_8005(jj,1) = haversine([lat3D(lat_8005_i),lon3D(jj)],[lat3D(lat_8005_i),lon3D(jj+1)]); 
% haversign calculates the geometric distance in km between two coordinates
end
lon_m_8005 = lon_km_8005*1000;
lon_m_8005(end+1,1) = lon_m_8005(end,1);

% now to calculate the half-way distance for the Volume transports
lon_diff_8005 = NaN(size(lon_m_8005))';

for kk = 1:length(lon_m_8005)-1
lon_diff_8005(kk) = (lon_m_8005(kk+1)/2)+(lon_m_8005(kk)/2);
end
lon_diff_8005(end) = lon_m_8005(end)/2;

% @ 81.45
lon_km_8145 = NaN(length(lon3D)-1,1);
for jj = 1:(length(lon3D)-1)
lon_km_8145(jj,1) = haversine([lat3D(lat_8145_i),lon3D(jj)],[lat3D(lat_8145_i),lon3D(jj+1)]); 
% haversign calculates the geometric distance in km between two coordinates
end
lon_m_8145 = lon_km_8145*1000;
lon_m_8145(end+1,1) = lon_m_8145(end,1);

% now to calculate the half-way distance for the Volume transports
lon_diff_8145 = NaN(size(lon_m_8145))';

for kk = 1:length(lon_m_8145)-1
lon_diff_8145(kk) = (lon_m_8145(kk+1)/2)+(lon_m_8145(kk)/2);
end
lon_diff_8145(end) = lon_m_8145(end)/2;

%% DZ
z_diff = NaN(size(z_star));

for kk = 2:length(z_star)-1
z_diff(kk) = ((z_star(kk+1)+z_star(1,kk))/2)-((z_star(kk)+z_star(1,kk-1))/2);
end
z_diff(1) = ((z_star(2)+z_star(1))/2)-(z_star(1)/2);
z_diff(end) = (z_star(end)-z_star(end-1));

%% Calculating the Volume Transport of the water fresher than 34
% at this point it isn't too standardized in terms of zonal extent, but
% since we're just trying to identify the EGC and the EGCC the (arbitrary)
% longitudes I chose for each section will suffice

addpath('/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_VxFW')

%% VX @ 70.95
load('600_7095.mat');

lon_7095_i = find(lon3D>-22 & lon3D<-13);
v_7095 = squeeze(v_4D(lon_7095_i,lat_7095_i,:,:));

clear vx_34_7095
for i = 1:204
vx_34_7095(:,:,i) = z_diff.*squeeze(v_7095(:,:,i)).*lon_diff_7095(lon_7095_i)';
end

% now to apply the mask of where there is salinities <= 34
vx_34_7095 = vx_34_7095.*yn_34_7095(lon_7095_i,:,:);

% lastly, to add up across the depths
vx_34_7095 = squeeze(sum(vx_34_7095,2,'omitnan'));

% convert to Sverdrups
vx_34_7095 = vx_34_7095/10^6;

% taking an average across the entire timeseries
vx_avg_7095 = mean(vx_34_7095,2,'omitnan');

figure
plot(lon3D(lon_7095_i),1000*vx_avg_7095,'-k','LineWidth',2)
xline([x70 x70],'-r','LineWidth',1)
ylim([-80 30])
xlabel('Longitude [^oE]')
ylabel('Volume Transport [mSv]')
title('Volume Transport of Fresh Water at 70.95 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/7095.jpg']);

%% VX @ 73.95
load('600_7395.mat');

lon_7395_i = find(lon3D>-20 & lon3D<-8);
v_7395 = squeeze(v_4D(lon_7395_i,lat_7395_i,:,:));

clear vx_34_7395
for i = 1:204
vx_34_7395(:,:,i) = z_diff.*squeeze(v_7395(:,:,i)).*lon_diff_7395(lon_7395_i)';
end

% now to apply the mask of where there is salinities <= 34
vx_34_7395 = vx_34_7395.*yn_34_7395(lon_7395_i,:,:);

% lastly, to add up across the depths
vx_34_7395 = squeeze(sum(vx_34_7395,2,'omitnan'));

% convert to Sverdrups
vx_34_7395 = vx_34_7395/10^6;

% taking an average across the entire timeseries
vx_avg_7395 = mean(vx_34_7395,2,'omitnan');

figure
plot(lon3D(lon_7395_i),1000*vx_avg_7395,'-k','LineWidth',2)
xline([x73 x73],'-r','LineWidth',1)
ylim([-40 30])
xlabel('Longitude [^oE]')
ylabel('Volume Transport [mSv]')
title('Volume Transport of Fresh Water at 73.95 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/7395.jpg']);

%% VX @ 76.95
load('600_7695.mat');

lon_7695_i = find(lon3D>-19 & lon3D<0);
v_7695 = squeeze(v_4D(lon_7695_i,lat_7695_i,:,:));

clear vx_34_7695
for i = 1:204
vx_34_7695(:,:,i) = z_diff.*squeeze(v_7695(:,:,i)).*lon_diff_7695(lon_7695_i)';
end

% now to apply the mask of where there is salinities <= 34
vx_34_7695 = vx_34_7695.*yn_34_7695(lon_7695_i,:,:);

% lastly, to add up across the depths
vx_34_7695 = squeeze(sum(vx_34_7695,2,'omitnan'));

% convert to Sverdrups
vx_34_7695 = vx_34_7695/10^6;

% taking an average across the entire timeseries
vx_avg_7695 = mean(vx_34_7695,2,'omitnan');

figure
plot(lon3D(lon_7695_i),1000*vx_avg_7695,'-k','LineWidth',2)
xline([x76 x76],'-r','LineWidth',1)
xlim([-19 0])
ylim([-40 30])

xlabel('Longitude [^oE]')
ylabel('Volume Transport [mSv]')
title('Volume Transport of Fresh Water at 76.95 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/7695.jpg']);

%% VX @ 78.85
load('600_7885.mat');

lon_7885_i = find(lon3D>-18 & lon3D<0);
v_7885 = squeeze(v_4D(lon_7885_i,lat_7885_i,:,:));

clear vx_34_7885
for i = 1:204
vx_34_7885(:,:,i) = z_diff.*squeeze(v_7885(:,:,i)).*lon_diff_7885(lon_7885_i)';
end

% now to apply the mask of where there is salinities <= 34
vx_34_7885 = vx_34_7885.*yn_34_7885(lon_7885_i,:,:);

% lastly, to add up across the depths
vx_34_7885 = squeeze(sum(vx_34_7885,2,'omitnan'));

% convert to Sverdrups
vx_34_7885 = vx_34_7885/10^6;

% taking an average across the entire timeseries
vx_avg_7885 = mean(vx_34_7885,2,'omitnan');

figure
plot(lon3D(lon_7885_i),1000*vx_avg_7885,'-k','LineWidth',2)
xline([x78 x78],'-r','LineWidth',1)
xlim([-18 0])
ylim([-40 30])

xlabel('Longitude [^oE]')
ylabel('Volume Transport [mSv]')
title('Volume Transport of Fresh Water at 78.85 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/7885.jpg']);

%% VX @ 80.05
load('600_8005.mat');

lon_8005_i = find(lon3D>-18 & lon3D<0);
v_8005 = squeeze(v_4D(lon_8005_i,lat_8005_i,:,:));

clear vx_34_8005
for i = 1:204
vx_34_8005(:,:,i) = z_diff.*squeeze(v_8005(:,:,i)).*lon_diff_8005(lon_8005_i)';
end

% now to apply the mask of where there is salinities <= 34
vx_34_8005 = vx_34_8005.*yn_34_8005(lon_8005_i,:,:);

% lastly, to add up across the depths
vx_34_8005 = squeeze(sum(vx_34_8005,2,'omitnan'));

% convert to Sverdrups
vx_34_8005 = vx_34_8005/10^6;

% taking an average across the entire timeseries
vx_avg_8005 = mean(vx_34_8005,2,'omitnan');

figure
plot(lon3D(lon_8005_i),1000*vx_avg_8005,'-k','LineWidth',2)
xline([x80 x80],'-r','LineWidth',1)
xlim([-18 0])
ylim([-40 30])

xlabel('Longitude [^oE]')
ylabel('Volume Transport [mSv]')
title('Volume Transport of Fresh Water at 80.05 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/8005.jpg']);

%% VX @ 81.45

load('600_8145.mat');

lon_8145_i = find(lon3D>-13 & lon3D<0);
v_8145 = squeeze(v_4D(lon_8145_i,lat_8145_i,:,:));

clear vx_34_8145
for i = 1:204
vx_34_8145(:,:,i) = z_diff.*squeeze(v_8145(:,:,i)).*lon_diff_8145(lon_8145_i)';
end

% now to apply the mask of where there is salinities <= 34
vx_34_8145 = vx_34_8145.*yn_34_8145(lon_8145_i,:,:);

% lastly, to add up across the depths
vx_34_8145 = squeeze(sum(vx_34_8145,2,'omitnan'));

% convert to Sverdrups
vx_34_8145 = vx_34_8145/10^6;

% taking an average across the entire timeseries
vx_avg_8145 = mean(vx_34_8145,2,'omitnan');

figure
plot(lon3D(lon_8145_i),1000*vx_avg_8145,'-k','LineWidth',2)
xlim([-13 0])
ylim([-40 30])

xline([x81 x81],'-r','LineWidth',1)
xlabel('Longitude [^oE]')
ylabel('Volume Transport [mSv]')
title('Volume Transport of Fresh Water at 81.45 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/8145.jpg']);


%% Plotting the Seasonal Changes in Volume Transport of Fresh Water

% finding the indicies of the individual months
jan_i = find(dv(:,2) == 1);
feb_i = find(dv(:,2) == 2);
mar_i = find(dv(:,2) == 3);
apr_i = find(dv(:,2) == 4);
may_i = find(dv(:,2) == 5);
jun_i = find(dv(:,2) == 6);
jul_i = find(dv(:,2) == 7);
aug_i = find(dv(:,2) == 8);
sep_i = find(dv(:,2) == 9);
oct_i = find(dv(:,2) == 10);
nov_i = find(dv(:,2) == 11);
dec_i = find(dv(:,2) == 12);

spring_i = [mar_i;apr_i;may_i];
summer_i = [jun_i;jul_i;aug_i];
fall_i = [sep_i;oct_i;nov_i];
winter_i = [dec_i;jan_i;feb_i];

%% Separating out at each latitude by season

% % @ 73.95N
% spring_vx_7395 = mean(vx_34_7395(:,spring_i),2,'omitnan');
% summer_vx_7395 = mean(vx_34_7395(:,summer_i),2,'omitnan');
% fall_vx_7395 = mean(vx_34_7395(:,fall_i),2,'omitnan');
% winter_vx_7395 = mean(vx_34_7395(:,winter_i),2,'omitnan');
% 
% figure
% plot(lon3D(lon_7395_i),1000*fall_vx_7395,'-k','LineWidth',2)
% xline([x73 x73],'-r','LineWidth',1)
% ylim([-45 30])
% xlabel('Longitude [^oE]')
% ylabel('Volume Transport [mSv]')
% title('Fall (SON) Volume Transport of Fresh Water at 73.95 ^oN (RARE)')
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_VxFW/fall_7395.jpg']);
% 
% %% @ 76.95N
% spring_vx_7695 = mean(vx_34_7695(:,spring_i),2,'omitnan');
% summer_vx_7695 = mean(vx_34_7695(:,summer_i),2,'omitnan');
% fall_vx_7695 = mean(vx_34_7695(:,fall_i),2,'omitnan');
% winter_vx_7695 = mean(vx_34_7695(:,winter_i),2,'omitnan');
% 
% figure
% plot(lon3D(lon_7695_i),1000*winter_vx_7695,'-k','LineWidth',2)
% xline([x76 x76],'-r','LineWidth',1)
% xlim([-19 0])
% ylim([-15 10])
% 
% xlabel('Longitude [^oE]')
% ylabel('Volume Transport [mSv]')
% % title('Spring (MAM) Volume Transport of Fresh Water at 76.95 ^oN (RARE)')
% % title('Summer (JJA) Volume Transport of Fresh Water at 76.95 ^oN (RARE)')
% % title('Fall (SON) Volume Transport of Fresh Water at 76.95 ^oN (RARE)')
% title('Winter (DJF) Volume Transport of Fresh Water at 76.95 ^oN (RARE)')
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_VxFW/winter_7695.jpg']);
% 
% %% @ 78.85N
% spring_vx_7885 = mean(vx_34_7885(:,spring_i),2,'omitnan');
% summer_vx_7885 = mean(vx_34_7885(:,summer_i),2,'omitnan');
% fall_vx_7885 = mean(vx_34_7885(:,fall_i),2,'omitnan');
% winter_vx_7885 = mean(vx_34_7885(:,winter_i),2,'omitnan');
% 
% figure
% plot(lon3D(lon_7885_i),1000*fall_vx_7885,'-k','LineWidth',2)
% xline([x78 x78],'-r','LineWidth',1)
% xlim([-18 0])
% ylim([-40 40])
% 
% xlabel('Longitude [^oE]')
% ylabel('Volume Transport [mSv]')
% % title('Spring (MAM) Volume Transport of Fresh Water at 78.85 ^oN (RARE)')
% % title('Summer (JJA) Volume Transport of Fresh Water at 78.85 ^oN (RARE)')
% title('Fall (SON) Volume Transport of Fresh Water at 78.85 ^oN (RARE)')
% % title('Winter (DJF) Volume Transport of Fresh Water at 78.85 ^oN (RARE)')
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_VxFW/fall_7885.jpg']);
% 
% %% @ 80.05N
% spring_vx_8005 = mean(vx_34_8005(:,spring_i),2,'omitnan');
% summer_vx_8005 = mean(vx_34_8005(:,summer_i),2,'omitnan');
% fall_vx_8005 = mean(vx_34_8005(:,fall_i),2,'omitnan');
% winter_vx_8005 = mean(vx_34_8005(:,winter_i),2,'omitnan');
% 
% figure
% plot(lon3D(lon_8005_i),1000*summer_vx_8005,'-k','LineWidth',2)
% xline([x80 x80],'-r','LineWidth',1)
% xlim([-18 0])
% ylim([-40 40])
% 
% xlabel('Longitude [^oE]')
% ylabel('Volume Transport [mSv]')
% % title('Spring (MAM) Volume Transport of Fresh Water at 80.05 ^oN (RARE)')
% title('Summer (JJA) Volume Transport of Fresh Water at 80.05 ^oN (RARE)')
% % title('Fall (SON) Volume Transport of Fresh Water at 80.05 ^oN (RARE)')
% % title('Winter (DJF) Volume Transport of Fresh Water at 80.05 ^oN (RARE)')
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_VxFW/summer_8005.jpg']);
% 
% %% @ 81.45N
% 
% spring_vx_8145 = mean(vx_34_8145(:,spring_i),2,'omitnan');
% summer_vx_8145 = mean(vx_34_8145(:,summer_i),2,'omitnan');
% fall_vx_8145 = mean(vx_34_8145(:,fall_i),2,'omitnan');
% winter_vx_8145 = mean(vx_34_8145(:,winter_i),2,'omitnan');
% 
% figure
% plot(lon3D(lon_8145_i),1000*fall_vx_8145,'-k','LineWidth',2)
% xline([x81 x81],'-r','LineWidth',1)
% xlim([-13 0])
% ylim([-15 10])
% 
% xlabel('Longitude [^oE]')
% ylabel('Volume Transport [mSv]')
% % title('Spring (MAM) Volume Transport of Fresh Water at 81.45 ^oN (RARE)')
% % title('Summer (JJA) Volume Transport of Fresh Water at 81.45 ^oN (RARE)')
% title('Fall (SON) Volume Transport of Fresh Water at 81.45 ^oN (RARE)')
% % title('Winter (DJF) Volume Transport of Fresh Water at 81.45 ^oN (RARE)')
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_VxFW/fall_8145.jpg']);

%% Now to plot: these are overall averages

%% @ 70.95N
% [Z,LON] = meshgrid(-z_star(1:24),lon3D(lon_7095_i));
% 
% v_7095_avg = squeeze(v_3D(lon_7095_i,lat_7095_i,1:24));
% alpha_7095_avg = mean(alpha_7095(lon_7095_i,1:24,:),3,'omitnan');
% 
% figure
% hold on
% contourf(LON,Z,v_7095_avg,8,'linestyle','none');
% clim([-0.2 0.2])
% colormap(cmocean('balance','pivot',0))
% x = colorbar;
% [C,d] = contour(LON,Z,alpha_7095_avg,[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% [C,d] = contour(LON,Z,mean(SA_7095(lon_7095_i,1:24),3,'omitnan'),[34 34],'-m','LineWidth',2);
% clabel(C,d,'Color','m');
% 
% % yline([-600 -600],'r','LineWidth',1)
% imAlpha=ones(size(v_7095_avg));
% imAlpha(isnan(v_7095_avg))=0;
% imagesc(v_7095_avg,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% ylim([-1000 0])
% xlim([-22 -14])
% 
% xlabel('Longitude [^oE]')
% ylabel('Depth [m]')
% ylabel(x,'V [m/s]')
% title('Mean NEGS V-Velocity, α, and 34 Isohaline at 70.95^oN (RARE)');
% 
% set(gcf,'color','w');
% 
% % [x70,y70] = ginput(1);
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE Vertical Sections/7095.jpg']);

% %% @ 73.95N
% [Z,LON] = meshgrid(-z_star(1:24),lon3D(lon_7395_i));
% 
% v_7395_avg = squeeze(v_3D(lon_7395_i,lat_7395_i,1:24));
% alpha_7395_avg = mean(alpha_7395(lon_7395_i,1:24,:),3,'omitnan');
% 
% figure
% hold on
% contourf(LON,Z,v_7395_avg,8,'linestyle','none');
% clim([-0.2 0.2])
% colormap(cmocean('balance','pivot',0))
% x = colorbar;
% [C,d] = contour(LON,Z,alpha_7395_avg,[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% [C,d] = contour(LON,Z,mean(SA_7395(lon_7395_i,1:24),3,'omitnan'),[34 34],'-m','LineWidth',2);
% clabel(C,d,'Color','m');
% 
% imAlpha=ones(size(v_7395_avg));
% imAlpha(isnan(v_7395_avg))=0;
% imagesc(v_7395_avg,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% ylim([-1000 0])
% xlim([-20 -8])
% 
% xlabel('Longitude [^oE]')
% ylabel('Depth [m]')
% ylabel(x,'V [m/s]')
% title('Mean NEGS V-Velocity, α, and 34 Isohaline at 73.95^oN (RARE)');
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE Vertical Sections/7395.jpg']);
% 
% %% @ 76.95
% [Z,LON] = meshgrid(-z_star(1:24),lon3D(lon_7695_i));
% 
% v_7695_avg = squeeze(v_3D(lon_7695_i,lat_7695_i,1:24));
% alpha_7695_avg = mean(alpha_7695(lon_7695_i,1:24,:),3,'omitnan');
% 
% figure
% hold on
% contourf(LON,Z,v_7695_avg,8,'linestyle','none');
% clim([-0.2 0.2])
% colormap(cmocean('balance','pivot',0))
% x = colorbar;
% [C,d] = contour(LON,Z,alpha_7695_avg,[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% [C,d] = contour(LON,Z,mean(SA_7695(lon_7695_i,1:24),3,'omitnan'),[34 34],'-m','LineWidth',2);
% clabel(C,d,'Color','m');
% 
% imAlpha=ones(size(v_7695_avg));
% imAlpha(isnan(v_7695_avg))=0;
% imagesc(v_7695_avg,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% ylim([-1000 0])
% xlim([-19 0])
% 
% xlabel('Longitude [^oE]')
% ylabel('Depth [m]')
% ylabel(x,'V [m/s]')
% title('Mean NEGS V-Velocity, α, and 34 Isohaline at 76.95^oN (RARE)');
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE Vertical Sections/7695.jpg']);
% 
% %% @ 78.85N
% [Z,LON] = meshgrid(-z_star(1:24),lon3D(lon_7885_i));
% 
% v_7885_avg = squeeze(v_3D(lon_7885_i,lat_7885_i,1:24));
% alpha_7885_avg = mean(alpha_7885(lon_7885_i,1:24,:),3,'omitnan');
% 
% figure
% hold on
% contourf(LON,Z,v_7885_avg,8,'linestyle','none');
% clim([-0.2 0.2])
% colormap(cmocean('balance','pivot',0))
% x = colorbar;
% [C,d] = contour(LON,Z,alpha_7885_avg,[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% [C,d] = contour(LON,Z,mean(SA_7885(lon_7885_i,1:24),3,'omitnan'),[34 34],'-m','LineWidth',2);
% clabel(C,d,'Color','m');
% 
% imAlpha=ones(size(v_7885_avg));
% imAlpha(isnan(v_7885_avg))=0;
% imagesc(v_7885_avg,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% ylim([-1000 0])
% xlim([-18 0])
% 
% xlabel('Longitude [^oE]')
% ylabel('Depth [m]')
% ylabel(x,'V [m/s]')
% title('Mean NEGS V-Velocity, α, and 34 Isohaline at 78.85^oN (RARE)');
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE Vertical Sections/7885.jpg']);
% 
% % @ 80.05N
% [Z,LON] = meshgrid(-z_star(1:24),lon3D(lon_8005_i));
% 
% v_8005_avg = squeeze(v_3D(lon_8005_i,lat_8005_i,1:24));
% alpha_8005_avg = mean(alpha_8005(lon_8005_i,1:24,:),3,'omitnan');
% 
% figure
% hold on
% contourf(LON,Z,v_8005_avg,8,'linestyle','none');
% clim([-0.2 0.2])
% colormap(cmocean('balance','pivot',0))
% x = colorbar;
% [C,d] = contour(LON,Z,alpha_8005_avg,[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% [C,d] = contour(LON,Z,mean(SA_8005(lon_8005_i,1:24),3,'omitnan'),[34 34],'-m','LineWidth',2);
% clabel(C,d,'Color','m');
% 
% imAlpha=ones(size(v_8005_avg));
% imAlpha(isnan(v_8005_avg))=0;
% imagesc(v_8005_avg,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% ylim([-1000 0])
% xlim([-18 0])
% 
% xlabel('Longitude [^oE]')
% ylabel('Depth [m]')
% ylabel(x,'V [m/s]')
% title('Mean NEGS V-Velocity, α, and 34 Isohaline at 80.05^oN (RARE)');
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE Vertical Sections/8005.jpg']);
% 
% %% @ 81.45N
% [Z,LON] = meshgrid(-z_star(1:24),lon3D(lon_8145_i));
% 
% v_8145_avg = squeeze(v_3D(lon_8145_i,lat_8145_i,1:24));
% alpha_8145_avg = mean(alpha_8145(lon_8145_i,1:24,:),3,'omitnan');
% 
% figure
% hold on
% contourf(LON,Z,v_8145_avg,8,'linestyle','none');
% clim([-0.2 0.2])
% colormap(cmocean('balance','pivot',0))
% x = colorbar;
% [C,d] = contour(LON,Z,alpha_8145_avg,[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% [C,d] = contour(LON,Z,mean(SA_8145(lon_8145_i,1:24),3,'omitnan'),[34 34],'-m','LineWidth',2);
% clabel(C,d,'Color','m');
% 
% imAlpha=ones(size(v_8145_avg));
% imAlpha(isnan(v_8145_avg))=0;
% imagesc(v_8145_avg,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% ylim([-1000 0])
% xlim([-13 0])
% 
% xlabel('Longitude [^oE]')
% ylabel('Depth [m]')
% ylabel(x,'V [m/s]')
% title('Mean NEGS V-Velocity, α, and 34 Isohaline at 81.45^oN (RARE)');
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE Vertical Sections/8145.jpg']);

%% @ 76.95N: Defining the EGC with Velocity

% focusing on just the 76.95N section right now. First to verticlly
% integrate over the top 1000m

for i = 1:204
    for j = 1:length(lon_7695_i)
v_integ_7695(j,i) = squeeze(sum((z_diff(1:24).*v_7695(j,1:24,i)),2,'omitnan'));
    end
end

% finding the minimum that we'll use as the core of the EGC
vint_avg_7695 = mean(v_integ_7695,2,'omitnan');
[core_egc,idx] = min(vint_avg_7695);

% finding the nearest longitude of the 600m isobath
idx600_7695 = find(abs(lon3D(lon_7695_i)-x76)==nanmin(abs(lon3D(lon_7695_i)-x76)),1,'first');

% now calculating the decay of the current using 1/e
core_egc_lim = core_egc*(1/exp(1));
idx_egc_limW = find(abs(vint_avg_7695(1:idx)-core_egc_lim)==nanmin(abs(vint_avg_7695(1:idx)-core_egc_lim)),1,'last');
idx_egc_limE = find(abs(vint_avg_7695(idx:end)-core_egc_lim)==nanmin(abs(vint_avg_7695(idx:end)-core_egc_lim)),1,'first');

lon_7695 = lon3D(lon_7695_i);

lon_egc_7695W = lon_7695(1:idx);
lon_egc_7695E = lon_7695(idx:end);

% plotting the time average for each longitude
figure
hold on
plot(lon_7695,vint_avg_7695,'-k','LineWidth',1)
scatter(lon3D(lon_7695_i(idx)),core_egc,200,'*b')
xline([x76 x76],'-r')
xline([lon_egc_7695E(idx_egc_limE) lon_egc_7695E(idx_egc_limE)],'-.b')
xline([lon_egc_7695W(idx_egc_limW) lon_egc_7695W(idx_egc_limW)],'-.b')

xlim([-19 0])
ylim([-120 25])

legend('','EGC Core','600 m isobath','','EGC limits','Location','southwest')

ylabel('Integrated Velocity [m^2 s^-^1]')
xlabel('Longitude [^oE]')
title('Integrated V-Velocity at 76.95 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Integrated Velocity/7695_avg.jpg']);

%% Now to find the EGCC
% should I look inshore of the 600m isobath? or onshore of the limit of the
% EGC... tried with inshore of the EGC decay and it wasn't low enough for
% the code to distinguish the onshelf max of southward flow. 
% Using the 600m isobath was also problamatic since there are moments where
% the max/e was met exactly, but visually it was offshore of the local
% velocity minimum (closest to 0)

% beginning with inshore of the EGC
% lonW_egc_76 = lon_decay_7695W(idx_egc_limW);

% first though I think I'm going to interpolate the top 150m so that it's
% on equal spacing

z_150 = [0:1:150];

v_7695_upper = interp3(z_star,lon_7695,[1:204],v_7695,z_150,lon_7695,[1:204]);
v_7695_upper_avg = mean(v_7695_upper,3);

% now to take the average in the top 150m and find the min inshore of the
% 600m isobath

v_7695_upper_avg = mean(v_7695_upper_avg,2,'omitnan');

[core_egcc,idx_cc] = min(v_7695_upper_avg(1:idx600_7695));

% identifying the min of velocity between the EGC and EGCC
[shelf_min,idx_min] = max(v_7695_upper_avg(idx_cc:idx));

% now calculating the decay of the current using 1/e
core_egcc_lim = core_egcc*(1/exp(1));
idx_egcc_limW = find(abs(v_7695_upper_avg(1:idx_cc)-core_egcc_lim)==nanmin(abs(v_7695_upper_avg(1:idx_cc)-core_egcc_lim)),1,'last');
idx_egcc_limE = find(abs(v_7695_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)==nanmin(abs(v_7695_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)),1,'first');

% applying the same limits on longitude so I can plot the x-lines
lon_egcc_7695W = lon_7695(1:idx_cc);
lon_egcc_7695E = lon_7695(idx_cc:end);

% and now identifying the max and limits of the EGCCC by looking for the
% velocity max inshore of the EGCC
[core_egccc,idx_ccc] = max(v_7695_upper_avg(1:idx_cc));

core_egccc_lim = core_egccc*(1/exp(1));
idx_egccc_limW = find(abs(v_7695_upper_avg(1:idx_ccc)-core_egccc_lim)==nanmin(abs(v_7695_upper_avg(1:idx_ccc)-core_egccc_lim)),1,'last');
idx_egccc_limE = find(abs(v_7695_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)==nanmin(abs(v_7695_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)),1,'first');

vx_avg_7695_egccc = sum(vx_avg_7695(idx_egccc_limW:(idx_ccc+idx_egccc_limE-1)),'all');

vx_avg_7695_egcc = sum(vx_avg_7695(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
vx_avg_7695_egc = sum(vx_avg_7695(idx_egc_limW:(idx+idx_egc_limE-1)),'all');

% applying the same limits on longitude 
lon_egccc_7695W = lon_7695(1:idx_ccc);
lon_egccc_7695E = lon_7695(idx_ccc:idx_cc);


%% plotting the average upp 150m flow with the EGC, EGCC, and EGCCC min

figure
hold on
plot(lon_7695,v_7695_upper_avg,'-k','LineWidth',1)
scatter(lon_7695(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_7695E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_7695(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x76 x76],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_7695E(idx_egc_limE) lon_egc_7695E(idx_egc_limE)],'-.b')
xline([lon_egc_7695W(idx_egc_limW) lon_egc_7695W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_7695E(idx_egcc_limE) lon_egcc_7695E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_7695W(idx_egcc_limW) lon_egcc_7695W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_7695E(idx_egccc_limE) lon_egccc_7695E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_7695W(idx_egccc_limW) lon_egccc_7695W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

xlim([-19 0])
ylim([-0.20 0.1])

legend('','EGCC Core','Velocity min','EGCCC Core','600 m isobath','','EGC limits','','',...
    '','EGCC limits','','','','EGCCC limits','Location','northeast')

ylabel('Velocity [m s^-^1]')
xlabel('Longitude [^oE]')
title('Average Upper 150m V-Velocity at 76.95 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Upper Velocity/7695_avg.jpg']);

%% Combining the Integrated Velocities and Upper 150m Velocities 76.95°N

yyaxis left
hold on
plot(lon_7695,vint_avg_7695,'--k','LineWidth',1)
scatter(lon3D(lon_7695_i(idx)),core_egc,200,'*b')
xline([x76 x76],'-r')
xline([lon_egc_7695E(idx_egc_limE) lon_egc_7695E(idx_egc_limE)],'-.b')
xline([lon_egc_7695W(idx_egc_limW) lon_egc_7695W(idx_egc_limW)],'-.b')

ylabel('Integrated Velocity [m^2 s^-^1]','FontSize',16)
ylim([-120 25])

% Plot data for the right y-axis
yyaxis right
plot(lon_7695,v_7695_upper_avg,'-k','LineWidth',1)
scatter(lon_7695(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_7695E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_7695(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x76 x76],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_7695E(idx_egc_limE) lon_egc_7695E(idx_egc_limE)],'-.b')
xline([lon_egc_7695W(idx_egc_limW) lon_egc_7695W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_7695E(idx_egcc_limE) lon_egcc_7695E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_7695W(idx_egcc_limW) lon_egcc_7695W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_7695E(idx_egccc_limE) lon_egccc_7695E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_7695W(idx_egccc_limW) lon_egccc_7695W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

ylabel('Velocity [m s^-^1]','FontSize',16)
ylim([-0.20 0.1])

legend('Integrated V','','','','','','','','Average 150 m V','Location','southwest','FontSize',12);

% Add common labels and title
xlim([-19 0])
xlabel('Longitude [^oE]')
title('76.95°N','FontSize',16);

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Combined/7695.jpg']);


%% Applying limits to the Volume Transport of Fresh Water: 76.95N

% these are the volumes (on average, seasonally, etc) of water 34 or
% fresher that the egcc is transporting at 76.95N (in Sv).

% reassigning the values of these variables to include the entire
% timeseries

% first to redefine the 150m average
v_7695_upper_avg = squeeze(mean(v_7695_upper,2,'omitnan'));

% then, in order to accurately define the limits of the various currents,
% I'm going to smooth the data for the EGCC and EGCCC limits

v_7695_smooth = NaN(size(v_7695_upper_avg));
for i = 1:204
v_7695_smooth(7:end,i) = smooth(v_7695_upper_avg(7:end,i),10);
end


% redefining some of the EGC variables using integrated velocity
for i = 1:204
% EGC
[core_egc(i),idx(i)] = min(v_integ_7695(:,i));
core_egc_lim(i) = core_egc(i)*(1/exp(1));

idx_egc_limW(i) = find(abs(v_integ_7695(1:idx(i),i)-core_egc_lim(i))==nanmin(abs(v_integ_7695(1:idx(i),i)-core_egc_lim(i))),1,'last');
idx_egc_limE(i) = find(abs(v_integ_7695(idx(i):end,i)-core_egc_lim(i))==nanmin(abs(v_integ_7695(idx(i):end,i)-core_egc_lim(i))),1,'first');
end

for i = 1:204
% EGCC
[core_egcc(i),idx_cc(i)] = min(v_7695_smooth(1:idx600_7695,i));
[shelf_min(i),idx_min(i)] = max(v_7695_upper_avg(idx_cc(i):idx(i),i));
[core_egccc(i),idx_ccc(i)] = max(v_7695_smooth(1:idx_cc(i),i));

% we need to find where the current crosses over to zero
v_7695_abs = abs(v_7695_upper_avg);
[shelf_0(i),idx_0(i)] = min(v_7695_abs(idx_ccc(i):idx_cc(i),i));

% now calculating the decay of the EGCC using 1/e
core_egcc_lim(i) = core_egcc(i)*(1/exp(1));
idx_egcc_limW(i) = find(abs(v_7695_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))==nanmin(abs(v_7695_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))),1,'last');
idx_egcc_limE(i) = find(abs(v_7695_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))==nanmin(abs(v_7695_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))),1,'first');
end

for i = 1:204
% EGCCC
core_egccc_lim(i) = core_egccc(i)*(1/exp(1));
idx_egccc_limW(i) = find(abs(v_7695_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))==nanmin(abs(v_7695_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))),1,'last');
idx_egccc_limE(i) = find(abs(v_7695_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))==nanmin(abs(v_7695_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))),1,'first');


% redefining the decay limits
idx_egccc_limE(i) = idx_ccc(i)+idx_egccc_limE(i)-1;
idx_egcc_limW(i) = idx_ccc(i)+idx_0(i)+idx_egcc_limW(i)-1;
idx_egcc_limE(i) = idx_cc(i)+idx_egcc_limE(i)-1;
idx_egc_limE(i) = idx(i)+idx_egc_limE(i)-1;
idx_min(i) = idx_cc(i)+idx_min(i);

vx_7695_egccc(i) = sum(vx_34_7695(idx_egccc_limW(i):idx_egccc_limE(i),i),1);
vx_7695_egcc(i) = sum(vx_34_7695(idx_egcc_limW(i):idx_egcc_limE(i),i),1);
vx_7695_egc(i) = sum(vx_34_7695(idx_egc_limW(i):idx_egc_limE(i),i),1);

end

vx_avg_7695_egc = mean(vx_7695_egc,'omitnan');
vx_avg_7695_egcc = mean(vx_7695_egcc,'omitnan');
vx_avg_7695_egccc = mean(vx_7695_egccc,'omitnan');

% winter_vx_7695_egcc = sum(winter_vx_7695(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% spring_vx_7695_egcc = sum(spring_vx_7695(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% summer_vx_7695_egcc = sum(summer_vx_7695(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% fall_vx_7695_egcc = sum(fall_vx_7695(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');

% imaging this transports

years = [{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012',...
    '2013','2014','2015','2016','2017','2018','2019'}];

figure
hold on
plot([1:204],vx_7695_egccc,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:204],vx_7695_egcc,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:204],vx_7695_egc,'Color',[0.7 0 0.4],'LineWidth',2)

yline([vx_avg_7695_egccc vx_avg_7695_egccc],':','Color',[0.5 0.9 0.8],'LineWidth',2)
yline([vx_avg_7695_egcc vx_avg_7695_egcc],':','Color',[0.5 0.5 1],'LineWidth',2)
yline([vx_avg_7695_egc vx_avg_7695_egc],':','Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12:204);
xlim([1 204])
xticklabels(years)

ylim([-2 1])

ylabel('Volume Transport [Sv]')
title('Volume Transport of Fresh Water by the East Greenland Current System at 76.95 ^oN (RARE)')
legend('','','',sprintf('Average_E_G_C_C_C = %.4f Sv', vx_avg_7695_egccc),'',...
    sprintf('Average_E_G_C_C = %.4f Sv', vx_avg_7695_egcc),'',...
    sprintf('Average_E_G_C = %.4f Sv', vx_avg_7695_egc),'location','southeast')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/EGCC/7695_avg.jpg']);

%% Now to look at the other sections %%

%% @ 70.95N: Defining the EGC with Velocity

% focusing on just the 70.95N section right now. First to verticlly
% integrate over the top 1000m

for i = 1:204
    for j = 1:length(lon_7095_i)
v_integ_7095(j,i) = squeeze(sum((z_diff(1:24).*v_7095(j,1:24,i)),2,'omitnan'));
    end
end

% finding the minimum that we'll use as the core of the EGC
vint_avg_7095 = mean(v_integ_7095,2,'omitnan');
[core_egc,idx] = min(vint_avg_7095);

% finding the nearest longitude of the 600m isobath
idx600_7095 = find(abs(lon3D(lon_7095_i)-x70)==nanmin(abs(lon3D(lon_7095_i)-x70)),1,'first');

% now calculating the decay of the current using 1/e
core_egc_lim = core_egc*(1/exp(1));
idx_egc_limW = find(abs(vint_avg_7095(1:idx)-core_egc_lim)==nanmin(abs(vint_avg_7095(1:idx)-core_egc_lim)),1,'last');
idx_egc_limE = find(abs(vint_avg_7095(idx:end)-core_egc_lim)==nanmin(abs(vint_avg_7095(idx:end)-core_egc_lim)),1,'first');

lon_7095 = lon3D(lon_7095_i);

lon_egc_7095W = lon_7095(1:idx);
lon_egc_7095E = lon_7095(idx:end);

% plotting the time average for each longitude
figure
hold on
plot(lon_7095,vint_avg_7095,'-k','LineWidth',1)
scatter(lon3D(lon_7095_i(idx)),core_egc,200,'*b')
xline([x70 x70],'-r')
xline([lon_egc_7095E(idx_egc_limE) lon_egc_7095E(idx_egc_limE)],'-.b')
xline([lon_egc_7095W(idx_egc_limW) lon_egc_7095W(idx_egc_limW)],'-.b')

xlim([-22 -13])
ylim([-120 25])

legend('','EGC Core','600 m isobath','','EGC limits','Location','southwest')

ylabel('Integrated Velocity [m^2 s^-^1]')
xlabel('Longitude [^oE]')
title('Integrated V-Velocity at 70.95 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Integrated Velocity/7095_avg.jpg']);

%% Now to find the EGCC
% should I look inshore of the 600m isobath? or onshore of the limit of the
% EGC... tried with inshore of the EGC decay and it wasn't low enough for
% the code to distinguish the onshelf max of southward flow. 
% Using the 600m isobath was also problamatic since there are moments where
% the max/e was met exactly, but visually it was offshore of the local
% velocity minimum (closest to 0)

% beginning with inshore of the EGC
% lonW_egc_70 = lon_decay_7095W(idx_egc_limW);

% first though I think I'm going to interpolate the top 150m so that it's
% on equal spacing

z_150 = [0:1:150];

v_7095_upper = interp3(z_star,lon_7095,[1:204],v_7095,z_150,lon_7095,[1:204]);
v_7095_upper_avg = mean(v_7095_upper,3);

% now to take the average in the top 150m and find the min inshore of the
% 600m isobath

v_7095_upper_avg = mean(v_7095_upper_avg,2,'omitnan');

[core_egcc,idx_cc] = min(v_7095_upper_avg(1:idx600_7095));

% identifying the min of velocity between the EGC and EGCC
[shelf_min,idx_min] = max(v_7095_upper_avg(idx_cc:idx));

% now calculating the decay of the current using 1/e
core_egcc_lim = core_egcc*(1/exp(1));
idx_egcc_limW = find(abs(v_7095_upper_avg(1:idx_cc)-core_egcc_lim)==nanmin(abs(v_7095_upper_avg(1:idx_cc)-core_egcc_lim)),1,'last');
idx_egcc_limE = find(abs(v_7095_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)==nanmin(abs(v_7095_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)),1,'first');

% applying the same limits on longitude so I can plot the x-lines
lon_egcc_7095W = lon_7095(1:idx_cc);
lon_egcc_7095E = lon_7095(idx_cc:end);

% % and now identifying the max and limits of the EGCCC by looking for the
% % velocity max inshore of the EGCC
% [core_egccc,idx_ccc] = max(v_7095_upper_avg(1:idx_cc));
% 
% core_egccc_lim = core_egccc*(1/exp(1));
% idx_egccc_limW = find(abs(v_7095_upper_avg(1:idx_ccc)-core_egccc_lim)==nanmin(abs(v_7095_upper_avg(1:idx_ccc)-core_egccc_lim)),1,'last');
% idx_egccc_limE = find(abs(v_7095_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)==nanmin(abs(v_7095_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)),1,'first');
% 
% vx_avg_7095_egccc = sum(vx_avg_7095(idx_egccc_limW:(idx_ccc+idx_egccc_limE-1)),'all');

vx_avg_7095_egcc = sum(vx_avg_7095(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
vx_avg_7095_egc = sum(vx_avg_7095(idx_egc_limW:(idx+idx_egc_limE-1)),'all');

% applying the same limits on longitude 
lon_egccc_7095W = lon_7095(1:idx_ccc);
lon_egccc_7095E = lon_7095(idx_ccc:idx_cc);


%% plotting the average upp 150m flow with the EGC, EGCC, and EGCCC min

figure
hold on
plot(lon_7095,v_7095_upper_avg,'-k','LineWidth',1)
scatter(lon_7095(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_7095E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
% scatter(lon_7095(idx_ccc),core_egccc,80,'^',...
%     'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x70 x70],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_7095E(idx_egc_limE) lon_egc_7095E(idx_egc_limE)],'-.b')
xline([lon_egc_7095W(idx_egc_limW) lon_egc_7095W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_7095E(idx_egcc_limE) lon_egcc_7095E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_7095W(idx_egcc_limW) lon_egcc_7095W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% % EGCCC limits
% xline([lon_egccc_7095E(idx_egccc_limE) lon_egccc_7095E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
% xline([lon_egccc_7095W(idx_egccc_limW) lon_egccc_7095W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

xlim([-22 -13])
ylim([-0.20 0.1])

legend('','EGCC Core','Velocity min','600 m isobath','','EGC limits','','',...
    '','EGCC limits','','','','Location','northeast')

ylabel('Velocity [m s^-^1]')
xlabel('Longitude [^oE]')
title('Average Upper 150m V-Velocity at 70.95 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Upper Velocity/7095_avg.jpg']);

%% Combining the Integrated Velocities and Upper 150m Velocities 70.95°N

yyaxis left
hold on
plot(lon_7095,vint_avg_7095,'--k','LineWidth',1)
scatter(lon3D(lon_7095_i(idx)),core_egc,200,'*b')
xline([x70 x70],'-r')
xline([lon_egc_7095E(idx_egc_limE) lon_egc_7095E(idx_egc_limE)],'-.b')
xline([lon_egc_7095W(idx_egc_limW) lon_egc_7095W(idx_egc_limW)],'-.b')

ylabel('Integrated Velocity [m^2 s^-^1]','FontSize',16)
ylim([-120 25])

% Plot data for the right y-axis
yyaxis right
plot(lon_7095,v_7095_upper_avg,'-k','LineWidth',1)
scatter(lon_7095(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_7095E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
% scatter(lon_7095(idx_ccc),core_egccc,80,'^',...
%     'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x70 x70],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_7095E(idx_egc_limE) lon_egc_7095E(idx_egc_limE)],'-.b')
xline([lon_egc_7095W(idx_egc_limW) lon_egc_7095W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_7095E(idx_egcc_limE) lon_egcc_7095E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_7095W(idx_egcc_limW) lon_egcc_7095W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
%xline([lon_egccc_7095E(idx_egccc_limE) lon_egccc_7095E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
%xline([lon_egccc_7095W(idx_egccc_limW) lon_egccc_7095W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

ylabel('Velocity [m s^-^1]','FontSize',16)
ylim([-0.20 0.1])

legend('Integrated V','','','','','','','','Average 150 m V','Location','southwest','FontSize',12);

% Add common labels and title
xlim([-22 -13])
xlabel('Longitude [^oE]')
title('70.95°N','FontSize',16);

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Combined/7095.jpg']);


%% Applying limits to the Volume Transport of Fresh Water: 70.95N

% these are the volumes (on average, seasonally, etc) of water 34 or
% fresher that the egcc is transporting at 70.95N (in Sv).

% reassigning the values of these variables to include the entire
% timeseries

% first to redefine the 150m average
v_7095_upper_avg = squeeze(mean(v_7095_upper,2,'omitnan'));

% then, in order to accurately define the limits of the various currents,
% I'm going to smooth the data for the EGCC and EGCCC limits

v_7095_smooth = NaN(size(v_7095_upper_avg));
for i = 1:204
v_7095_smooth(7:end,i) = smooth(v_7095_upper_avg(7:end,i),10);
end


% redefining some of the EGC variables using integrated velocity
for i = 1:204
% EGC
[core_egc(i),idx(i)] = min(v_integ_7095(:,i));
core_egc_lim(i) = core_egc(i)*(1/exp(1));

idx_egc_limW(i) = find(abs(v_integ_7095(1:idx(i),i)-core_egc_lim(i))==nanmin(abs(v_integ_7095(1:idx(i),i)-core_egc_lim(i))),1,'last');
idx_egc_limE(i) = find(abs(v_integ_7095(idx(i):end,i)-core_egc_lim(i))==nanmin(abs(v_integ_7095(idx(i):end,i)-core_egc_lim(i))),1,'first');
end

for i = 1:204
% EGCC
[core_egcc(i),idx_cc(i)] = min(v_7095_smooth(1:idx600_7095,i));
[shelf_min(i),idx_min(i)] = max(v_7095_upper_avg(idx_cc(i):idx(i),i));
% [core_egccc(i),idx_ccc(i)] = max(v_7095_smooth(1:idx_cc(i),i));

% we need to find where the current crosses over to zero
v_7095_abs = abs(v_7095_upper_avg);
[shelf_0(i),idx_0(i)] = min(v_7095_abs(1:idx_cc(i),i));

% now calculating the decay of the EGCC using 1/e
core_egcc_lim(i) = core_egcc(i)*(1/exp(1));
idx_egcc_limW(i) = find(abs(v_7095_smooth(idx_0(i):idx_cc(i),i)-core_egcc_lim(i))==nanmin(abs(v_7095_smooth(idx_0(i):idx_cc(i),i)-core_egcc_lim(i))),1,'last');
idx_egcc_limE(i) = find(abs(v_7095_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))==nanmin(abs(v_7095_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))),1,'first');
end

for i = 1:204
% % EGCCC
% core_egccc_lim(i) = core_egccc(i)*(1/exp(1));
% idx_egccc_limW(i) = find(abs(v_7095_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))==nanmin(abs(v_7095_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))),1,'last');
% idx_egccc_limE(i) = find(abs(v_7095_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))==nanmin(abs(v_7095_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))),1,'first');
% 

% redefining the decay limits
% idx_egccc_limE(i) = idx_ccc(i)+idx_egccc_limE(i)-1;
idx_egcc_limW(i) = idx_0(i)+idx_egcc_limW(i)-1;
idx_egcc_limE(i) = idx_cc(i)+idx_egcc_limE(i)-1;
idx_egc_limE(i) = idx(i)+idx_egc_limE(i)-1;
idx_min(i) = idx_cc(i)+idx_min(i);

% vx_7095_egccc(i) = sum(vx_34_7095(idx_egccc_limW(i):idx_egccc_limE(i),i),1);
vx_7095_egcc(i) = sum(vx_34_7095(idx_egcc_limW(i):idx_egcc_limE(i),i),1);
vx_7095_egc(i) = sum(vx_34_7095(idx_egc_limW(i):idx_egc_limE(i),i),1);

end

vx_avg_7095_egc = mean(vx_7095_egc,'omitnan');
vx_avg_7095_egcc = mean(vx_7095_egcc,'omitnan');
% vx_avg_7095_egccc = mean(vx_7095_egccc,'omitnan');

% winter_vx_7095_egcc = sum(winter_vx_7095(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% spring_vx_7095_egcc = sum(spring_vx_7095(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% summer_vx_7095_egcc = sum(summer_vx_7095(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% fall_vx_7095_egcc = sum(fall_vx_7095(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');

% imaging this transports

years = [{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012',...
    '2013','2014','2015','2016','2017','2018','2019'}];

figure
hold on
% plot([1:204],vx_7095_egccc,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:204],vx_7095_egcc,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:204],vx_7095_egc,'Color',[0.7 0 0.4],'LineWidth',2)

% yline([vx_avg_7095_egccc vx_avg_7095_egccc],':','Color',[0.5 0.9 0.8],'LineWidth',2)
yline([vx_avg_7095_egcc vx_avg_7095_egcc],':','Color',[0.5 0.5 1],'LineWidth',2)
yline([vx_avg_7095_egc vx_avg_7095_egc],':','Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12:204);
xlim([1 204])
xticklabels(years)

ylim([-2 1])

ylabel('Volume Transport [Sv]')
title('Volume Transport of Fresh Water by the East Greenland Current System at 70.95 ^oN (RARE)')
legend('','','',sprintf('Average_E_G_C_C = %.4f Sv', vx_avg_7095_egcc),'',...
    sprintf('Average_E_G_C = %.4f Sv', vx_avg_7095_egc),'location','southeast')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/EGCC/7095_avg.jpg']);

%% @ 73.95N: Defining the EGC with Velocity

% focusing on just the 73.95N section right now. First to verticlly
% integrate over the top 1000m

for i = 1:204
    for j = 1:length(lon_7395_i)
v_integ_7395(j,i) = squeeze(sum((z_diff(1:24).*v_7395(j,1:24,i)),2,'omitnan'));
    end
end

% finding the minimum that we'll use as the core of the EGC
vint_avg_7395 = mean(v_integ_7395,2,'omitnan');
[core_egc,idx] = min(vint_avg_7395);

% finding the nearest longitude of the 600m isobath
idx600_7395 = find(abs(lon3D(lon_7395_i)-x73)==nanmin(abs(lon3D(lon_7395_i)-x73)),1,'first');

% now calculating the decay of the current using 1/e
core_egc_lim = core_egc*(1/exp(1));
idx_egc_limW = find(abs(vint_avg_7395(1:idx)-core_egc_lim)==nanmin(abs(vint_avg_7395(1:idx)-core_egc_lim)),1,'last');
idx_egc_limE = find(abs(vint_avg_7395(idx:end)-core_egc_lim)==nanmin(abs(vint_avg_7395(idx:end)-core_egc_lim)),1,'first');


lon_7395 = lon3D(lon_7395_i);

lon_egc_7395W = lon_7395(1:idx);
lon_egc_7395E = lon_7395(idx:end);

% plotting the time average for each longitude
figure
hold on
plot(lon_7395,vint_avg_7395,'-k','LineWidth',1)
scatter(lon3D(lon_7395_i(idx)),core_egc,200,'*b')
xline([x73 x73],'-r')
xline([lon_egc_7395E(idx_egc_limE) lon_egc_7395E(idx_egc_limE)],'-.b')
xline([lon_egc_7395W(idx_egc_limW) lon_egc_7395W(idx_egc_limW)],'-.b')

xlim([-20 -8])
ylim([-120 25])

legend('','EGC Core','600 m isobath','','EGC limits','Location','southwest')

ylabel('Integrated Velocity [m^2 s^-^1]')
xlabel('Longitude [^oE]')
title('Integrated V-Velocity at 73.95 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Integrated Velocity/7395_avg.jpg']);

%% Now to find the EGCC

z_150 = [0:1:150];

v_7395_upper = interp3(z_star,lon_7395,[1:204],v_7395,z_150,lon_7395,[1:204]);
v_7395_upper_avg = mean(v_7395_upper,3);

% now to take the average in the top 150m and find the min inshore of the EGC

v_7395_upper_avg = mean(v_7395_upper_avg,2,'omitnan');

[core_egcc,idx_cc] = min(v_7395_upper_avg(1:idx600_7395));
[shelf_min,idx_min] = max(v_7395_upper_avg(idx_cc:idx)); % this is where a definition 
% comes into play. Finding this minimum based on what we're calling the
% shelf current. Is it confined to inshore of the 600m isobath? Or just
% inshore of the EGC decay limit?

% now calculating the decay of the current using 1/e
core_egcc_lim = core_egcc*(1/exp(1));
idx_egcc_limW = find(abs(v_7395_upper_avg(1:idx_cc)-core_egcc_lim)==nanmin(abs(v_7395_upper_avg(1:idx_cc)-core_egcc_lim)),1,'last');
idx_egcc_limE = find(abs(v_7395_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)==nanmin(abs(v_7395_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)),1,'first');

% applying the same limits on longitude so I can plot the x-lines
lon_egcc_7395W = lon_7395(1:idx_cc);
lon_egcc_7395E = lon_7395(idx_cc:end);
% and now identifying the max and limits of the EGCCC by looking for the
% velocity max inshore of the EGCC
[core_egccc,idx_ccc] = max(v_7395_upper_avg(1:idx_cc));

core_egccc_lim = core_egccc*(1/exp(1));
idx_egccc_limW = find(abs(v_7395_upper_avg(1:idx_ccc)-core_egccc_lim)==nanmin(abs(v_7395_upper_avg(1:idx_ccc)-core_egccc_lim)),1,'last');
idx_egccc_limE = find(abs(v_7395_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)==nanmin(abs(v_7395_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)),1,'first');

vx_avg_7395_egccc = sum(vx_avg_7395(idx_egccc_limW:(idx_ccc+idx_egccc_limE-1)),'all');

vx_avg_7395_egcc = sum(vx_avg_7395(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
vx_avg_7395_egc = sum(vx_avg_7395(idx_egc_limW:(idx+idx_egc_limE-1)),'all');

% applying the same limits on longitude 
lon_egccc_7395W = lon_7395(1:idx_ccc);
lon_egccc_7395E = lon_7395(idx_ccc:idx_cc);


%% plotting the average upp 150m flow with the EGC, EGCC, and EGCCC min

figure
hold on
plot(lon_7395,v_7395_upper_avg,'-k','LineWidth',1)
scatter(lon_7395(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_7395E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_7395(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x73 x73],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_7395E(idx_egc_limE) lon_egc_7395E(idx_egc_limE)],'-.b')
xline([lon_egc_7395W(idx_egc_limW) lon_egc_7395W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_7395E(idx_egcc_limE) lon_egcc_7395E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_7395W(idx_egcc_limW) lon_egcc_7395W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_7395E(idx_egccc_limE) lon_egccc_7395E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_7395W(idx_egccc_limW) lon_egccc_7395W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

xlim([-20 -8])
ylim([-0.20 0.1])

legend('','EGCC Core','Velocity min','EGCCC Core','600 m isobath','','EGC limits','','',...
    '','EGCC limits','','','','EGCCC limits','Location','northeast')

ylabel('Velocity [m s^-^1]')
xlabel('Longitude [^oE]')
title('Average Upper 150m V-Velocity at 73.95 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Upper Velocity/7395_avg.jpg']);

%% Combining the Integrated Velocities and Upper 150m Velocities 73.95°N

yyaxis left
hold on
plot(lon_7395,vint_avg_7395,'--k','LineWidth',1)
scatter(lon3D(lon_7395_i(idx)),core_egc,200,'*b')
xline([x73 x73],'-r')
xline([lon_egc_7395E(idx_egc_limE) lon_egc_7395E(idx_egc_limE)],'-.b')
xline([lon_egc_7395W(idx_egc_limW) lon_egc_7395W(idx_egc_limW)],'-.b')

ylabel('Integrated Velocity [m^2 s^-^1]','FontSize',16)
ylim([-120 25])

% Plot data for the right y-axis
yyaxis right
plot(lon_7395,v_7395_upper_avg,'-k','LineWidth',1)
scatter(lon_7395(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_7395E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_7395(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x73 x73],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_7395E(idx_egc_limE) lon_egc_7395E(idx_egc_limE)],'-.b')
xline([lon_egc_7395W(idx_egc_limW) lon_egc_7395W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_7395E(idx_egcc_limE) lon_egcc_7395E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_7395W(idx_egcc_limW) lon_egcc_7395W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_7395E(idx_egccc_limE) lon_egccc_7395E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_7395W(idx_egccc_limW) lon_egccc_7395W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

ylabel('Velocity [m s^-^1]','FontSize',16)
ylim([-0.20 0.1])

legend('Integrated V','','','','','','','','Average 150 m V','Location','southwest','FontSize',12);

% Add common labels and title
xlim([-20 -8])
xlabel('Longitude [^oE]')
title('73.95°N','FontSize',16);

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Combined/7395.jpg']);


%% Applying limits to the Volume Transport of Fresh Water: 73.95N

% these are the volumes (on average, seasonally, etc) of water 34 or
% fresher that the egcc is transporting at 73.95N (in Sv).

% reassigning the values of these variables to include the entire
% timeseries

% first to redefine the 150m average
v_7395_upper_avg = squeeze(mean(v_7395_upper,2,'omitnan'));

% then, in order to accurately define the limits of the various currents,
% I'm going to smooth the data for the EGCC and EGCCC limits

v_7395_smooth = NaN(size(v_7395_upper_avg));
for i = 1:204
v_7395_smooth(7:end,i) = smooth(v_7395_upper_avg(7:end,i),10);
end


% redefining some of the EGC variables using integrated velocity
for i = 1:204
% EGC
[core_egc(i),idx(i)] = min(v_integ_7395(:,i));
core_egc_lim(i) = core_egc(i)*(1/exp(1));

idx_egc_limW(i) = find(abs(v_integ_7395(1:idx(i),i)-core_egc_lim(i))==nanmin(abs(v_integ_7395(1:idx(i),i)-core_egc_lim(i))),1,'last');
idx_egc_limE(i) = find(abs(v_integ_7395(idx(i):end,i)-core_egc_lim(i))==nanmin(abs(v_integ_7395(idx(i):end,i)-core_egc_lim(i))),1,'first');
end

for i = 1:204
% EGCC
[core_egcc(i),idx_cc(i)] = min(v_7395_smooth(1:idx600_7395,i));
[shelf_min(i),idx_min(i)] = max(v_7395_upper_avg(idx_cc(i):idx(i),i));
[core_egccc(i),idx_ccc(i)] = max(v_7395_smooth(1:idx_cc(i),i));

% we need to find where the current crosses over to zero
v_7395_abs = abs(v_7395_upper_avg);
[shelf_0(i),idx_0(i)] = min(v_7395_abs(idx_ccc(i):idx_cc(i),i));

% now calculating the decay of the EGCC using 1/e
core_egcc_lim(i) = core_egcc(i)*(1/exp(1));
idx_egcc_limW(i) = find(abs(v_7395_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))==nanmin(abs(v_7395_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))),1,'last');
idx_egcc_limE(i) = find(abs(v_7395_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))==nanmin(abs(v_7395_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))),1,'first');
end

for i = 1:204
% EGCCC
core_egccc_lim(i) = core_egccc(i)*(1/exp(1));
idx_egccc_limW(i) = find(abs(v_7395_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))==nanmin(abs(v_7395_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))),1,'last');
idx_egccc_limE(i) = find(abs(v_7395_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))==nanmin(abs(v_7395_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))),1,'first');


% redefining the decay limits
idx_egccc_limE(i) = idx_ccc(i)+idx_egccc_limE(i)-1;
idx_egcc_limW(i) = idx_ccc(i)+idx_0(i)+idx_egcc_limW(i)-1;
idx_egcc_limE(i) = idx_cc(i)+idx_egcc_limE(i)-1;
idx_egc_limE(i) = idx(i)+idx_egc_limE(i)-1;
idx_min(i) = idx_cc(i)+idx_min(i);

vx_7395_egccc(i) = sum(vx_34_7395(idx_egccc_limW(i):idx_egccc_limE(i),i),1);
vx_7395_egcc(i) = sum(vx_34_7395(idx_egcc_limW(i):idx_egcc_limE(i),i),1);
vx_7395_egc(i) = sum(vx_34_7395(idx_egc_limW(i):idx_egc_limE(i),i),1);

end

vx_avg_7395_egc = mean(vx_7395_egc,'omitnan');
vx_avg_7395_egcc = mean(vx_7395_egcc,'omitnan');
vx_avg_7395_egccc = mean(vx_7395_egccc,'omitnan');

% winter_vx_7395_egcc = sum(winter_vx_7395(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% spring_vx_7395_egcc = sum(spring_vx_7395(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% summer_vx_7395_egcc = sum(summer_vx_7395(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% fall_vx_7395_egcc = sum(fall_vx_7395(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');

% imaging this transports

years = [{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012',...
    '2013','2014','2015','2016','2017','2018','2019'}];

figure
hold on
plot([1:204],vx_7395_egccc,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:204],vx_7395_egcc,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:204],vx_7395_egc,'Color',[0.7 0 0.4],'LineWidth',2)

yline([vx_avg_7395_egccc vx_avg_7395_egccc],':','Color',[0.5 0.9 0.8],'LineWidth',2)
yline([vx_avg_7395_egcc vx_avg_7395_egcc],':','Color',[0.5 0.5 1],'LineWidth',2)
yline([vx_avg_7395_egc vx_avg_7395_egc],':','Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12:204);
xlim([1 204])
xticklabels(years)

ylim([-2 1])

ylabel('Volume Transport [Sv]')
title('Volume Transport of Fresh Water by the East Greenland Current System at 73.95 ^oN (RARE)')
legend('','','',sprintf('Average_E_G_C_C_C = %.4f Sv', vx_avg_7395_egccc),'',...
    sprintf('Average_E_G_C_C = %.4f Sv', vx_avg_7395_egcc),'',...
    sprintf('Average_E_G_C = %.4f Sv', vx_avg_7395_egc),'location','southeast')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/EGCC/7395_avg.jpg']);

%% @ 78.85N: Defining the EGC with Velocity 

% focusing on just the 78.85N section right now. First to verticlly
% integrate over the top 1000m

for i = 1:204
    for j = 1:length(lon_7885_i)
v_integ_7885(j,i) = squeeze(sum((z_diff(1:24).*v_7885(j,1:24,i)),2,'omitnan'));
    end
end

% finding the minimum that we'll use as the core of the EGC
vint_avg_7885 = mean(v_integ_7885,2,'omitnan');
[core_egc,idx] = min(vint_avg_7885);

% finding the nearest longitude of the 600m isobath
idx600_7885 = find(abs(lon3D(lon_7885_i)-x78)==nanmin(abs(lon3D(lon_7885_i)-x78)),1,'first');

% now calculating the decay of the current using 1/e
core_egc_lim = core_egc*(1/exp(1));
idx_egc_limW = find(abs(vint_avg_7885(1:idx)-core_egc_lim)==nanmin(abs(vint_avg_7885(1:idx)-core_egc_lim)),1,'last');
idx_egc_limE = find(abs(vint_avg_7885(idx:end)-core_egc_lim)==nanmin(abs(vint_avg_7885(idx:end)-core_egc_lim)),1,'first');

lon_7885 = lon3D(lon_7885_i);

lon_egc_7885W = lon_7885(1:idx);
lon_egc_7885E = lon_7885(idx:end);

% plotting the time average for each longitude
figure
hold on
plot(lon_7885,vint_avg_7885,'-k','LineWidth',1)
scatter(lon3D(lon_7885_i(idx)),core_egc,200,'*b')
xline([x78 x78],'-r')
xline([lon_egc_7885E(idx_egc_limE) lon_egc_7885E(idx_egc_limE)],'-.b')
xline([lon_egc_7885W(idx_egc_limW) lon_egc_7885W(idx_egc_limW)],'-.b')

xlim([-18 0])
ylim([-120 25])

legend('','EGC Core','600 m isobath','','EGC limits','Location','southwest')

ylabel('Integrated Velocity [m^2 s^-^1]')
xlabel('Longitude [^oE]')
title('Integrated V-Velocity at 78.85 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Integrated Velocity/7885_avg.jpg']);

%% Now to find the EGCC

z_150 = [0:1:150];

v_7885_upper = interp3(z_star,lon_7885,[1:204],v_7885,z_150,lon_7885,[1:204]);
v_7885_upper_avg = mean(v_7885_upper,3);

% now to take the average in the top 150m and find the min inshore of the EGC

v_7885_upper_avg = mean(v_7885_upper_avg,2,'omitnan');

[core_egcc,idx_cc] = min(v_7885_upper_avg(1:idx600_7885));
[shelf_min,idx_min] = max(v_7885_upper_avg(idx_cc:idx));

% now calculating the decay of the current using 1/e
core_egcc_lim = core_egcc*(1/exp(1));
idx_egcc_limW = find(abs(v_7885_upper_avg(1:idx_cc)-core_egcc_lim)==nanmin(abs(v_7885_upper_avg(1:idx_cc)-core_egcc_lim)),1,'last');
idx_egcc_limE = find(abs(v_7885_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)==nanmin(abs(v_7885_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)),1,'first');

% applying the same limits on longitude so I can plot the x-lines
lon_egcc_7885W = lon_7885(1:idx_cc);
lon_egcc_7885E = lon_7885(idx_cc:end);
% and now identifying the max and limits of the EGCCC by looking for the
% velocity max inshore of the EGCC
[core_egccc,idx_ccc] = max(v_7885_upper_avg(1:idx_cc));

core_egccc_lim = core_egccc*(1/exp(1));
idx_egccc_limW = find(abs(v_7885_upper_avg(1:idx_ccc)-core_egccc_lim)==nanmin(abs(v_7885_upper_avg(1:idx_ccc)-core_egccc_lim)),1,'last');
idx_egccc_limE = find(abs(v_7885_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)==nanmin(abs(v_7885_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)),1,'first');

vx_avg_7885_egccc = sum(vx_avg_7885(idx_egccc_limW:(idx_ccc+idx_egccc_limE-1)),'all');

vx_avg_7885_egcc = sum(vx_avg_7885(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
vx_avg_7885_egc = sum(vx_avg_7885(idx_egc_limW:(idx+idx_egc_limE-1)),'all');

% applying the same limits on longitude 
lon_egccc_7885W = lon_7885(1:idx_ccc);
lon_egccc_7885E = lon_7885(idx_ccc:idx_cc);


%% plotting the average upp 150m flow with the EGC, EGCC, and EGCCC min

figure
hold on
plot(lon_7885,v_7885_upper_avg,'-k','LineWidth',1)
scatter(lon_7885(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_7885E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_7885(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x78 x78],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_7885E(idx_egc_limE) lon_egc_7885E(idx_egc_limE)],'-.b')
xline([lon_egc_7885W(idx_egc_limW) lon_egc_7885W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_7885E(idx_egcc_limE) lon_egcc_7885E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_7885W(idx_egcc_limW) lon_egcc_7885W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_7885E(idx_egccc_limE) lon_egccc_7885E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_7885W(idx_egccc_limW) lon_egccc_7885W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

xlim([-18 0])
ylim([-0.20 0.1])

legend('','EGCC Core','Velocity min','EGCCC Core','600 m isobath','','EGC limits','','',...
    '','EGCC limits','','','','EGCCC limits','Location','northeast')

ylabel('Velocity [m s^-^1]')
xlabel('Longitude [^oE]')
title('Average Upper 150m V-Velocity at 78.85 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Upper Velocity/7885_avg.jpg']);


%% Combining the Integrated Velocities and Upper 150m Velocities 78.85°N

yyaxis left
hold on
plot(lon_7885,vint_avg_7885,'--k','LineWidth',1)
scatter(lon3D(lon_7885_i(idx)),core_egc,200,'*b')
xline([x78 x78],'-r')
xline([lon_egc_7885E(idx_egc_limE) lon_egc_7885E(idx_egc_limE)],'-.b')
xline([lon_egc_7885W(idx_egc_limW) lon_egc_7885W(idx_egc_limW)],'-.b')

ylabel('Integrated Velocity [m^2 s^-^1]','FontSize',16)
ylim([-120 25])

% Plot data for the right y-axis
yyaxis right
plot(lon_7885,v_7885_upper_avg,'-k','LineWidth',1)
scatter(lon_7885(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_7885E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_7885(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x78 x78],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_7885E(idx_egc_limE) lon_egc_7885E(idx_egc_limE)],'-.b')
xline([lon_egc_7885W(idx_egc_limW) lon_egc_7885W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_7885E(idx_egcc_limE) lon_egcc_7885E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_7885W(idx_egcc_limW) lon_egcc_7885W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_7885E(idx_egccc_limE) lon_egccc_7885E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_7885W(idx_egccc_limW) lon_egccc_7885W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

ylabel('Velocity [m s^-^1]','FontSize',16)
ylim([-0.20 0.1])

legend('Integrated V','','','','','','','','Average 150 m V','Location','southwest','FontSize',12);

% Add common labels and title
xlim([-18 0])
xlabel('Longitude [^oE]')
title('78.85°N','FontSize',16);

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Combined/7885.jpg']);

%% Applying limits to the Volume Transport of Fresh Water: 78.85N

% these are the volumes (on average, seasonally, etc) of water 34 or
% fresher that the egcc is transporting at 78.85N (in Sv).

% reassigning the values of these variables to include the entire
% timeseries

% first to redefine the 150m average
v_7885_upper_avg = squeeze(mean(v_7885_upper,2,'omitnan'));

% then, in order to accurately define the limits of the various currents,
% I'm going to smooth the data for the EGCC and EGCCC limits

v_7885_smooth = NaN(size(v_7885_upper_avg));
for i = 1:204
v_7885_smooth(7:end,i) = smooth(v_7885_upper_avg(7:end,i),10);
end


% redefining some of the EGC variables using integrated velocity
for i = 1:204
% EGC
[core_egc(i),idx(i)] = min(v_integ_7885(:,i));
core_egc_lim(i) = core_egc(i)*(1/exp(1));

idx_egc_limW(i) = find(abs(v_integ_7885(1:idx(i),i)-core_egc_lim(i))==nanmin(abs(v_integ_7885(1:idx(i),i)-core_egc_lim(i))),1,'last');
idx_egc_limE(i) = find(abs(v_integ_7885(idx(i):end,i)-core_egc_lim(i))==nanmin(abs(v_integ_7885(idx(i):end,i)-core_egc_lim(i))),1,'first');
end

for i = 1:204
% EGCC
[core_egcc(i),idx_cc(i)] = min(v_7885_smooth(1:idx600_7885,i));
[shelf_min(i),idx_min(i)] = max(v_7885_upper_avg(idx_cc(i):idx(i),i));
[core_egccc(i),idx_ccc(i)] = max(v_7885_smooth(1:idx_cc(i),i));

% we need to find where the current crosses over to zero
v_7885_abs = abs(v_7885_upper_avg);
[shelf_0(i),idx_0(i)] = min(v_7885_abs(idx_ccc(i):idx_cc(i),i));

% now calculating the decay of the EGCC using 1/e
core_egcc_lim(i) = core_egcc(i)*(1/exp(1));
idx_egcc_limW(i) = find(abs(v_7885_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))==nanmin(abs(v_7885_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))),1,'last');
idx_egcc_limE(i) = find(abs(v_7885_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))==nanmin(abs(v_7885_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))),1,'first');
end

for i = 1:204
% EGCCC
core_egccc_lim(i) = core_egccc(i)*(1/exp(1));
idx_egccc_limW(i) = find(abs(v_7885_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))==nanmin(abs(v_7885_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))),1,'last');
idx_egccc_limE(i) = find(abs(v_7885_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))==nanmin(abs(v_7885_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))),1,'first');


% redefining the decay limits
idx_egccc_limE(i) = idx_ccc(i)+idx_egccc_limE(i)-1;
idx_egcc_limW(i) = idx_ccc(i)+idx_0(i)+idx_egcc_limW(i)-1;
idx_egcc_limE(i) = idx_cc(i)+idx_egcc_limE(i)-1;
idx_egc_limE(i) = idx(i)+idx_egc_limE(i)-1;
idx_min(i) = idx_cc(i)+idx_min(i);

vx_7885_egccc(i) = sum(vx_34_7885(idx_egccc_limW(i):idx_egccc_limE(i),i),1);
vx_7885_egcc(i) = sum(vx_34_7885(idx_egcc_limW(i):idx_egcc_limE(i),i),1);
vx_7885_egc(i) = sum(vx_34_7885(idx_egc_limW(i):idx_egc_limE(i),i),1);

end

vx_avg_7885_egc = mean(vx_7885_egc,'omitnan');
vx_avg_7885_egcc = mean(vx_7885_egcc,'omitnan');
vx_avg_7885_egccc = mean(vx_7885_egccc,'omitnan');

% winter_vx_7885_egcc = sum(winter_vx_7885(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% spring_vx_7885_egcc = sum(spring_vx_7885(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% summer_vx_7885_egcc = sum(summer_vx_7885(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% fall_vx_7885_egcc = sum(fall_vx_7885(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');

% imaging this transports

years = [{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012',...
    '2013','2014','2015','2016','2017','2018','2019'}];

figure
hold on
plot([1:204],vx_7885_egccc,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:204],vx_7885_egcc,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:204],vx_7885_egc,'Color',[0.7 0 0.4],'LineWidth',2)

yline([vx_avg_7885_egccc vx_avg_7885_egccc],':','Color',[0.5 0.9 0.8],'LineWidth',2)
yline([vx_avg_7885_egcc vx_avg_7885_egcc],':','Color',[0.5 0.5 1],'LineWidth',2)
yline([vx_avg_7885_egc vx_avg_7885_egc],':','Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12:204);
xlim([1 204])
xticklabels(years)

ylim([-2 1])

ylabel('Volume Transport [Sv]')
title('Volume Transport of Fresh Water by the East Greenland Current System at 78.85 ^oN (RARE)')
legend('','','',sprintf('Average_E_G_C_C_C = %.4f Sv', vx_avg_7885_egccc),'',...
    sprintf('Average_E_G_C_C = %.4f Sv', vx_avg_7885_egcc),'',...
    sprintf('Average_E_G_C = %.4f Sv', vx_avg_7885_egc),'location','southeast')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/EGCC/7885_avg.jpg']);

%% @ 80.05N: Defining the EGC with Velocity

% focusing on just the 80.05N section right now. First to verticlly
% integrate over the top 1000m

for i = 1:204
    for j = 1:length(lon_8005_i)
v_integ_8005(j,i) = squeeze(sum((z_diff(1:24).*v_8005(j,1:24,i)),2,'omitnan'));
    end
end

% finding the minimum that we'll use as the core of the EGC
vint_avg_8005 = mean(v_integ_8005,2,'omitnan');
[core_egc,idx] = min(vint_avg_8005);

% finding the nearest longitude of the 600m isobath
idx600_8005 = find(abs(lon3D(lon_8005_i)-x80)==nanmin(abs(lon3D(lon_8005_i)-x80)),1,'first');

% now calculating the decay of the current using 1/e
core_egc_lim = core_egc*(1/exp(1));
idx_egc_limW = find(abs(vint_avg_8005(1:idx)-core_egc_lim)==nanmin(abs(vint_avg_8005(1:idx)-core_egc_lim)),1,'last');
idx_egc_limE = find(abs(vint_avg_8005(idx:end)-core_egc_lim)==nanmin(abs(vint_avg_8005(idx:end)-core_egc_lim)),1,'first');

lon_8005 = lon3D(lon_8005_i);

lon_egc_8005W = lon_8005(1:idx);
lon_egc_8005E = lon_8005(idx:end);

% plotting the time average for each longitude
figure
hold on
plot(lon_8005,vint_avg_8005,'-k','LineWidth',1)
scatter(lon3D(lon_8005_i(idx)),core_egc,200,'*b')
xline([x80 x80],'-r')
xline([lon_egc_8005E(idx_egc_limE) lon_egc_8005E(idx_egc_limE)],'-.b')
xline([lon_egc_8005W(idx_egc_limW) lon_egc_8005W(idx_egc_limW)],'-.b')

xlim([-18 0])
ylim([-120 25])

legend('','AOO Core','600 m isobath','','AOO limits','Location','southwest')

ylabel('Integrated Velocity [m^2 s^-^1]')
xlabel('Longitude [^oE]')
title('Integrated V-Velocity at 80.05 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Integrated Velocity/8005_avg.jpg']);

%% Now to find the EGCC

z_150 = [0:1:150];

v_8005_upper = interp3(z_star,lon_8005,[1:204],v_8005,z_150,lon_8005,[1:204]);
v_8005_upper_avg = mean(v_8005_upper,3);

% now to take the average in the top 150m and find the min inshore of the EGC

v_8005_upper_avg = mean(v_8005_upper_avg,2,'omitnan');

[core_egcc,idx_cc] = min(v_8005_upper_avg(1:idx600_8005));
[shelf_min,idx_min] = max(v_8005_upper_avg(idx_cc:idx));

% now calculating the decay of the current using 1/e
core_egcc_lim = core_egcc*(1/exp(1));
idx_egcc_limW = find(abs(v_8005_upper_avg(1:idx_cc)-core_egcc_lim)==nanmin(abs(v_8005_upper_avg(1:idx_cc)-core_egcc_lim)),1,'last');
idx_egcc_limE = find(abs(v_8005_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)==nanmin(abs(v_8005_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)),1,'first');

% applying the same limits on longitude so I can plot the x-lines
lon_egcc_8005W = lon_8005(1:idx_cc);
lon_egcc_8005E = lon_8005(idx_cc:end);

% and now identifying the max and limits of the EGCCC by looking for the
% velocity max inshore of the EGCC
[core_egccc,idx_ccc] = max(v_8005_upper_avg(1:idx_cc));

core_egccc_lim = core_egccc*(1/exp(1));
idx_egccc_limW = find(abs(v_8005_upper_avg(1:idx_ccc)-core_egccc_lim)==nanmin(abs(v_8005_upper_avg(1:idx_ccc)-core_egccc_lim)),1,'last');
idx_egccc_limE = find(abs(v_8005_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)==nanmin(abs(v_8005_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)),1,'first');

vx_avg_8005_egccc = sum(vx_avg_8005(idx_egccc_limW:(idx_ccc+idx_egccc_limE-1)),'all');

vx_avg_8005_egcc = sum(vx_avg_8005(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
vx_avg_8005_egc = sum(vx_avg_8005(idx_egc_limW:(idx+idx_egc_limE-1)),'all');

% applying the same limits on longitude 
lon_egccc_8005W = lon_8005(1:idx_ccc);
lon_egccc_8005E = lon_8005(idx_ccc:idx_cc);


%% plotting the average upp 150m flow with the EGC, EGCC, and EGCCC min

figure
hold on
plot(lon_8005,v_8005_upper_avg,'-k','LineWidth',1)
scatter(lon_8005(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_8005E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_8005(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x80 x80],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_8005E(idx_egc_limE) lon_egc_8005E(idx_egc_limE)],'-.b')
xline([lon_egc_8005W(idx_egc_limW) lon_egc_8005W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_8005E(idx_egcc_limE) lon_egcc_8005E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_8005W(idx_egcc_limW) lon_egcc_8005W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_8005E(idx_egccc_limE) lon_egccc_8005E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_8005W(idx_egccc_limW) lon_egccc_8005W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

xlim([-18 0])
ylim([-0.20 0.1])

legend('','EGCC Core','Velocity min','EGCCC Core','600 m isobath','','EGC limits','','',...
    '','EGCC limits','','','','EGCCC limits','Location','northeast')

ylabel('Velocity [m s^-^1]')
xlabel('Longitude [^oE]')
title('Average Upper 150m V-Velocity at 80.05 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Upper Velocity/8005_avg.jpg']);

%% Combining the Integrated Velocities and Upper 150m Velocities 80.05°N

yyaxis left
hold on
plot(lon_8005,vint_avg_8005,'--k','LineWidth',1)
scatter(lon3D(lon_8005_i(idx)),core_egc,200,'*b')
xline([x80 x80],'-r')
xline([lon_egc_8005E(idx_egc_limE) lon_egc_8005E(idx_egc_limE)],'-.b')
xline([lon_egc_8005W(idx_egc_limW) lon_egc_8005W(idx_egc_limW)],'-.b')

ylabel('Integrated Velocity [m^2 s^-^1]','FontSize',16)
ylim([-120 25])

% Plot data for the right y-axis
yyaxis right
plot(lon_8005,v_8005_upper_avg,'-k','LineWidth',1)
scatter(lon_8005(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_8005E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_8005(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x80 x80],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_8005E(idx_egc_limE) lon_egc_8005E(idx_egc_limE)],'-.b')
xline([lon_egc_8005W(idx_egc_limW) lon_egc_8005W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_8005E(idx_egcc_limE) lon_egcc_8005E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_8005W(idx_egcc_limW) lon_egcc_8005W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_8005E(idx_egccc_limE) lon_egccc_8005E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_8005W(idx_egccc_limW) lon_egccc_8005W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

ylabel('Velocity [m s^-^1]','FontSize',16)
ylim([-0.20 0.1])

legend('Integrated V','','','','','','','','Average 150 m V','Location','southwest','FontSize',12);

% Add common labels and title
xlim([-18 0])
xlabel('Longitude [^oE]')
title('80.05°N','FontSize',16);

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Combined/8005.jpg']);

%% Applying limits to the Volume Transport of Fresh Water: 80.05N

% these are the volumes (on average, seasonally, etc) of water 34 or
% fresher that the egcc is transporting at 80.05N (in Sv).

% reassigning the values of these variables to include the entire
% timeseries

% first to redefine the 150m average
v_8005_upper_avg = squeeze(mean(v_8005_upper,2,'omitnan'));

% then, in order to accurately define the limits of the various currents,
% I'm going to smooth the data for the EGCC and EGCCC limits

v_8005_smooth = NaN(size(v_8005_upper_avg));
for i = 1:204
v_8005_smooth(7:end,i) = smooth(v_8005_upper_avg(7:end,i),10);
end


% redefining some of the EGC variables using integrated velocity
for i = 1:204
% EGC
[core_egc(i),idx(i)] = min(v_integ_8005(:,i));
core_egc_lim(i) = core_egc(i)*(1/exp(1));

idx_egc_limW(i) = find(abs(v_integ_8005(1:idx(i),i)-core_egc_lim(i))==nanmin(abs(v_integ_8005(1:idx(i),i)-core_egc_lim(i))),1,'last');
idx_egc_limE(i) = find(abs(v_integ_8005(idx(i):end,i)-core_egc_lim(i))==nanmin(abs(v_integ_8005(idx(i):end,i)-core_egc_lim(i))),1,'first');
end


for i = 1:204
% EGCC
[core_egcc(i),idx_cc(i)] = min(v_8005_smooth(1:idx_egc_limW(i),i));
[shelf_min(i),idx_min(i)] = max(v_8005_upper_avg(idx_cc(i):idx(i),i));
[core_egccc(i),idx_ccc(i)] = max(v_8005_smooth(1:idx_cc(i),i));

% we need to find where the current crosses over to zero
v_8005_abs = abs(v_8005_upper_avg);
[shelf_0(i),idx_0(i)] = min(v_8005_abs(idx_ccc(i):idx_cc(i),i));

% now calculating the decay of the EGCC using 1/e
core_egcc_lim(i) = core_egcc(i)*(1/exp(1));
idx_egcc_limW(i) = find(abs(v_8005_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))==nanmin(abs(v_8005_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))),1,'last');
idx_egcc_limE(i) = find(abs(v_8005_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))==nanmin(abs(v_8005_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))),1,'first');
end

for i = 1:204
% EGCCC
core_egccc_lim(i) = core_egccc(i)*(1/exp(1));
idx_egccc_limW(i) = find(abs(v_8005_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))==nanmin(abs(v_8005_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))),1,'last');
idx_egccc_limE(i) = find(abs(v_8005_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))==nanmin(abs(v_8005_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))),1,'first');


% redefining the decay limits
idx_egccc_limE(i) = idx_ccc(i)+idx_egccc_limE(i)-1;
idx_egcc_limW(i) = idx_ccc(i)+idx_0(i)+idx_egcc_limW(i)-1;
idx_egcc_limE(i) = idx_cc(i)+idx_egcc_limE(i)-1;
idx_egc_limE(i) = idx(i)+idx_egc_limE(i)-1;
idx_min(i) = idx_cc(i)+idx_min(i);

vx_8005_egccc(i) = sum(vx_34_8005(idx_egccc_limW(i):idx_egccc_limE(i),i),1);
vx_8005_egcc(i) = sum(vx_34_8005(idx_egcc_limW(i):idx_egcc_limE(i),i),1);
vx_8005_egc(i) = sum(vx_34_8005(idx_egc_limW(i):idx_egc_limE(i),i),1);

end

vx_avg_8005_egc = mean(vx_8005_egc,'omitnan');
vx_avg_8005_egcc = mean(vx_8005_egcc,'omitnan');
vx_avg_8005_egccc = mean(vx_8005_egccc,'omitnan');

% winter_vx_8005_egcc = sum(winter_vx_8005(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% spring_vx_8005_egcc = sum(spring_vx_8005(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% summer_vx_8005_egcc = sum(summer_vx_8005(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% fall_vx_8005_egcc = sum(fall_vx_8005(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');

% imaging this transports

years = [{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012',...
    '2013','2014','2015','2016','2017','2018','2019'}];

figure
hold on
plot([1:204],vx_8005_egccc,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:204],vx_8005_egcc,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:204],vx_8005_egc,'Color',[0.7 0 0.4],'LineWidth',2)

yline([vx_avg_8005_egccc vx_avg_8005_egccc],':','Color',[0.5 0.9 0.8],'LineWidth',2)
yline([vx_avg_8005_egcc vx_avg_8005_egcc],':','Color',[0.5 0.5 1],'LineWidth',2)
yline([vx_avg_8005_egc vx_avg_8005_egc],':','Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12:204);
xlim([1 204])
xticklabels(years)

ylim([-2 1])

ylabel('Volume Transport [Sv]')
title('Volume Transport of Fresh Water by the East Greenland Current System at 80.05 ^oN (RARE)')
legend('','','',sprintf('Average_E_G_C_C_C = %.4f Sv', vx_avg_8005_egccc),'',...
    sprintf('Average_E_G_C_C = %.4f Sv', vx_avg_8005_egcc),'',...
    sprintf('Average_E_G_C = %.4f Sv', vx_avg_8005_egc),'location','southeast')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/EGCC/8005_avg.jpg']);

%% @ 81.45N: Defining the EGC with Velocity

% focusing on just the 81.45N section right now. First to verticlly
% integrate over the top 1000m

for i = 1:204
    for j = 1:length(lon_8145_i)
v_integ_8145(j,i) = squeeze(sum((z_diff(1:24).*v_8145(j,1:24,i)),2,'omitnan'));
    end
end

% finding the minimum that we'll use as the core of the EGC
vint_avg_8145 = mean(v_integ_8145,2,'omitnan');
[core_egc,idx] = min(vint_avg_8145);

% finding the nearest longitude of the 600m isobath
idx600_8145 = find(abs(lon3D(lon_8145_i)-x81)==nanmin(abs(lon3D(lon_8145_i)-x81)),1,'first');

% now calculating the decay of the current using 1/e
core_egc_lim = core_egc*(1/exp(1));
idx_egc_limW = find(abs(vint_avg_8145(1:idx)-core_egc_lim)==nanmin(abs(vint_avg_8145(1:idx)-core_egc_lim)),1,'last');
idx_egc_limE = find(abs(vint_avg_8145(idx:end)-core_egc_lim)==nanmin(abs(vint_avg_8145(idx:end)-core_egc_lim)),1,'first');

lon_8145 = lon3D(lon_8145_i);

lon_egc_8145W = lon_8145(1:idx);
lon_egc_8145E = lon_8145(idx:end);

% plotting the time average for each longitude
figure
hold on
plot(lon_8145,vint_avg_8145,'-k','LineWidth',1)
scatter(lon3D(lon_8145_i(idx)),core_egc,200,'*b')
xline([x81 x81],'-r')
xline([lon_egc_8145E(idx_egc_limE) lon_egc_8145E(idx_egc_limE)],'-.b')
xline([lon_egc_8145W(idx_egc_limW) lon_egc_8145W(idx_egc_limW)],'-.b')

xlim([-13 0])
ylim([-80 40])

legend('','AOO Core','600 m isobath','','AOO limits','Location','southwest')

ylabel('Integrated Velocity [m^2 s^-^1]')
xlabel('Longitude [^oE]')
title('Integrated V-Velocity at 81.45 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Integrated Velocity/8145_avg.jpg']);

%% Now to find the EGCC

z_150 = [0:1:150];

v_8145_upper = interp3(z_star,lon_8145,[1:204],v_8145,z_150,lon_8145,[1:204]);
v_8145_upper_avg = mean(v_8145_upper,3);

% now to take the average in the top 150m and find the min inshore of the EGC

v_8145_upper_avg = mean(v_8145_upper_avg,2,'omitnan');

[core_egcc,idx_cc] = min(v_8145_upper_avg(1:idx600_8145));
[shelf_min,idx_min] = max(v_8145_upper_avg(idx_cc:idx));

% now calculating the decay of the current using 1/e
core_egcc_lim = core_egcc*(1/exp(1));
idx_egcc_limW = find(abs(v_8145_upper_avg(1:idx_cc)-core_egcc_lim)==nanmin(abs(v_8145_upper_avg(1:idx_cc)-core_egcc_lim)),1,'last');
idx_egcc_limE = find(abs(v_8145_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)==nanmin(abs(v_8145_upper_avg(idx_cc:(idx_cc+idx_min))-core_egcc_lim)),1,'first');

% applying the same limits on longitude so I can plot the x-lines
lon_egcc_8145W = lon_8145(1:idx_cc);
lon_egcc_8145E = lon_8145(idx_cc:end);

% and now identifying the max and limits of the EGCCC by looking for the
% velocity max inshore of the EGCC
[core_egccc,idx_ccc] = max(v_8145_upper_avg(1:idx_cc));

core_egccc_lim = core_egccc*(1/exp(1));
idx_egccc_limW = find(abs(v_8145_upper_avg(1:idx_ccc)-core_egccc_lim)==nanmin(abs(v_8145_upper_avg(1:idx_ccc)-core_egccc_lim)),1,'last');
idx_egccc_limE = find(abs(v_8145_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)==nanmin(abs(v_8145_upper_avg(idx_ccc:idx_cc)-core_egccc_lim)),1,'first');

vx_avg_8145_egccc = sum(vx_avg_8145(idx_egccc_limW:(idx_ccc+idx_egccc_limE-1)),'all');

vx_avg_8145_egcc = sum(vx_avg_8145(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
vx_avg_8145_egc = sum(vx_avg_8145(idx_egc_limW:(idx+idx_egc_limE-1)),'all');

% applying the same limits on longitude 
lon_egccc_8145W = lon_8145(1:idx_ccc);
lon_egccc_8145E = lon_8145(idx_ccc:idx_cc);


%% plotting the average upp 150m flow with the EGC, EGCC, and EGCCC min

figure
hold on
plot(lon_8145,v_8145_upper_avg,'-k','LineWidth',1)
scatter(lon_8145(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_8145E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_8145(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x81 x81],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_8145E(idx_egc_limE) lon_egc_8145E(idx_egc_limE)],'-.b')
xline([lon_egc_8145W(idx_egc_limW) lon_egc_8145W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_8145E(idx_egcc_limE) lon_egcc_8145E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_8145W(idx_egcc_limW) lon_egcc_8145W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_8145E(idx_egccc_limE) lon_egccc_8145E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_8145W(idx_egccc_limW) lon_egccc_8145W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

xlim([-13 0])
ylim([-0.20 0.1])

legend('','EGCC Core','Velocity min','EGCCC Core','600 m isobath','','EGC limits','','',...
    '','EGCC limits','','','','EGCCC limits','Location','northeast')

ylabel('Velocity [m s^-^1]')
xlabel('Longitude [^oE]')
title('Average Upper 150m V-Velocity at 81.45 ^oN (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Upper Velocity/8145_avg.jpg']);

% now I'll be going back to the individual sections and applying the
% definitions of the EGCC (and EGC maybe) to calculate the VT of Fresh Water

%% Combining the Integrated Velocities and Upper 150m Velocities 81.45°N

yyaxis left
hold on
plot(lon_8145,vint_avg_8145,'--k','LineWidth',1)
scatter(lon3D(lon_8145_i(idx)),core_egc,200,'*b')
xline([x81 x81],'-r')
xline([lon_egc_8145E(idx_egc_limE) lon_egc_8145E(idx_egc_limE)],'-.b')
xline([lon_egc_8145W(idx_egc_limW) lon_egc_8145W(idx_egc_limW)],'-.b')

ylabel('Integrated Velocity [m^2 s^-^1]','FontSize',16)
ylim([-80 40])

% Plot data for the right y-axis
yyaxis right
plot(lon_8145,v_8145_upper_avg,'-k','LineWidth',1)
scatter(lon_8145(idx_cc),core_egcc,100,'pentagram',...
    'MarkerEdgeColor',[0.1 0.7 0.2],'LineWidth',1.5) % the EGCC core
scatter(lon_egcc_8145E(idx_min),shelf_min,80,'o',...
    'MarkerEdgeColor',[0.5 .1 .1],'LineWidth',1.5) % the shelf min
scatter(lon_8145(idx_ccc),core_egccc,80,'^',...
    'MarkerEdgeColor',[0.8 0 1],'LineWidth',1.5) % EGCCC core

xline([x81 x81],'-r','LineWidth',1.5) % 600m isobath

% EGC limits
xline([lon_egc_8145E(idx_egc_limE) lon_egc_8145E(idx_egc_limE)],'-.b')
xline([lon_egc_8145W(idx_egc_limW) lon_egc_8145W(idx_egc_limW)],'-.b')

% EGCC limits
xline([lon_egcc_8145E(idx_egcc_limE) lon_egcc_8145E(idx_egcc_limE)],':','Color',[0.1 0.7 0.2],'LineWidth',2)
xline([lon_egcc_8145W(idx_egcc_limW) lon_egcc_8145W(idx_egcc_limW)],':','Color',[0.1 0.7 0.2],'LineWidth',2)

% EGCCC limits
xline([lon_egccc_8145E(idx_egccc_limE) lon_egccc_8145E(idx_egccc_limE)],'--','Color',[0.8 0 1],'LineWidth',1)
xline([lon_egccc_8145W(idx_egccc_limW) lon_egccc_8145W(idx_egccc_limW)],'--','Color',[0.8 0 1],'LineWidth',1)

ylabel('Velocity [m s^-^1]','FontSize',16)
ylim([-0.20 0.1])

legend('Integrated V','','','','','','','','Average 150 m V','Location','southwest','FontSize',12);

% Add common labels and title
xlim([-13 0])
xlabel('Longitude [^oE]')
title('81.45°N','FontSize',16);

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE Vertical Sections/Combined/8145.jpg']);


%% Applying limits to the Volume Transport of Fresh Water: 81.45N

% these are the volumes (on average, seasonally, etc) of water 34 or
% fresher that the egcc is transporting at 81.45N (in Sv).

% reassigning the values of these variables to include the entire
% timeseries

% first to redefine the 150m average
v_8145_upper_avg = squeeze(mean(v_8145_upper,2,'omitnan'));

% then, in order to accurately define the limits of the various currents,
% I'm going to smooth the data for the EGCC and EGCCC limits

v_8145_smooth = NaN(size(v_8145_upper_avg));
for i = 1:204
v_8145_smooth(7:end,i) = smooth(v_8145_upper_avg(7:end,i),10);
end


% redefining some of the EGC variables using integrated velocity
for i = 1:204
% EGC
[core_egc(i),idx(i)] = min(v_integ_8145(:,i));
core_egc_lim(i) = core_egc(i)*(1/exp(1));

idx_egc_limW(i) = find(abs(v_integ_8145(1:idx(i),i)-core_egc_lim(i))==nanmin(abs(v_integ_8145(1:idx(i),i)-core_egc_lim(i))),1,'last');
idx_egc_limE(i) = find(abs(v_integ_8145(idx(i):end,i)-core_egc_lim(i))==nanmin(abs(v_integ_8145(idx(i):end,i)-core_egc_lim(i))),1,'first');
end

for i = 1:204
% EGCC
[core_egcc(i),idx_cc(i)] = min(v_8145_smooth(1:idx600_8145,i));
[shelf_min(i),idx_min(i)] = max(v_8145_upper_avg(idx_cc(i):idx(i),i));
[core_egccc(i),idx_ccc(i)] = max(v_8145_smooth(1:idx_cc(i),i));

% we need to find where the current crosses over to zero
v_8145_abs = abs(v_8145_upper_avg);
[shelf_0(i),idx_0(i)] = min(v_8145_abs(idx_ccc(i):idx_cc(i),i));

% now calculating the decay of the EGCC using 1/e
core_egcc_lim(i) = core_egcc(i)*(1/exp(1));
idx_egcc_limW(i) = find(abs(v_8145_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))==nanmin(abs(v_8145_smooth((idx_ccc(i)+idx_0(i)-1):idx_cc(i),i)-core_egcc_lim(i))),1,'last');
idx_egcc_limE(i) = find(abs(v_8145_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))==nanmin(abs(v_8145_smooth(idx_cc(i):(idx_cc(i)+idx_min(i)-1),i)-core_egcc_lim(i))),1,'first');
end

for i = 1:204
% EGCCC
core_egccc_lim(i) = core_egccc(i)*(1/exp(1));
idx_egccc_limW(i) = find(abs(v_8145_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))==nanmin(abs(v_8145_upper_avg(1:idx_ccc(i),i)-core_egccc_lim(i))),1,'last');
idx_egccc_limE(i) = find(abs(v_8145_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))==nanmin(abs(v_8145_upper_avg(idx_ccc(i):(idx_ccc(i)+idx_0(i)),i)-core_egccc_lim(i))),1,'first');


% redefining the decay limits
idx_egccc_limE(i) = idx_ccc(i)+idx_egccc_limE(i)-1;
idx_egcc_limW(i) = idx_ccc(i)+idx_0(i)+idx_egcc_limW(i)-1;
idx_egcc_limE(i) = idx_cc(i)+idx_egcc_limE(i)-1;
idx_egc_limE(i) = idx(i)+idx_egc_limE(i)-1;
idx_min(i) = idx_cc(i)+idx_min(i);

vx_8145_egccc(i) = sum(vx_34_8145(idx_egccc_limW(i):idx_egccc_limE(i),i),1);
vx_8145_egcc(i) = sum(vx_34_8145(idx_egcc_limW(i):idx_egcc_limE(i),i),1);
vx_8145_egc(i) = sum(vx_34_8145(idx_egc_limW(i):idx_egc_limE(i),i),1);

end

vx_avg_8145_egc = mean(vx_8145_egc,'omitnan');
vx_avg_8145_egcc = mean(vx_8145_egcc,'omitnan');
vx_avg_8145_egccc = mean(vx_8145_egccc,'omitnan');

% winter_vx_8145_egcc = sum(winter_vx_8145(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% spring_vx_8145_egcc = sum(spring_vx_8145(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% summer_vx_8145_egcc = sum(summer_vx_8145(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');
% fall_vx_8145_egcc = sum(fall_vx_8145(idx_egcc_limW:(idx_cc+idx_egcc_limE-1)),'all');

% imaging this transports

years = [{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012',...
    '2013','2014','2015','2016','2017','2018','2019'}];

figure
hold on
plot([1:204],vx_8145_egccc,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:204],vx_8145_egcc,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:204],vx_8145_egc,'Color',[0.7 0 0.4],'LineWidth',2)

yline([vx_avg_8145_egccc vx_avg_8145_egccc],':','Color',[0.5 0.9 0.8],'LineWidth',2)
yline([vx_avg_8145_egcc vx_avg_8145_egcc],':','Color',[0.5 0.5 1],'LineWidth',2)
yline([vx_avg_8145_egc vx_avg_8145_egc],':','Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12:204);
xlim([1 204])
xticklabels(years)

ylim([-2 1])

ylabel('Volume Transport [Sv]')
title('Volume Transport of Fresh Water by the East Greenland Current System at 81.45 ^oN (RARE)')
legend('','','',sprintf('Average_E_G_C_C_C = %.4f Sv', vx_avg_8145_egccc),'',...
    sprintf('Average_E_G_C_C = %.4f Sv', vx_avg_8145_egcc),'',...
    sprintf('Average_E_G_C = %.4f Sv', vx_avg_8145_egc),'location','southeast')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/EGCC/8145_avg.jpg']);

%% Correlating the EGCC Vx and EGC Vx
% [R_7095,P_7095] = corrcoef(vx_7095_egcc_sea,vx_7095_egc_sea);
% [R_7395,P_7395] = corrcoef(vx_7395_egcc_sea,vx_7395_egc_sea);
% [R_7695,P_7695] = corrcoef(vx_7695_egcc_sea,vx_7695_egc_sea);
% [R_7885,P_7885] = corrcoef(vx_7885_egcc_sea,vx_7885_egc_sea);
% [R_8005,P_8005] = corrcoef(vx_8005_egcc_sea,vx_8005_egc_sea);
% [R_8145,P_8145] = corrcoef(vx_8145_egcc_sea,vx_8145_egc_sea);

%% Total Vx of Fresh Water across the 5 sections
% however, I think I'm going to limit the total fwt to the EGC and EGCC
% system

vx_FW_8145 = sum(vx_34_8145,1);
vx_FW_8005 = sum(vx_34_8005,1);
vx_FW_7885 = sum(vx_34_7885,1);
vx_FW_7695 = sum(vx_34_7695,1);
vx_FW_7395 = sum(vx_34_7395,1);
vx_FW_7095 = sum(vx_34_7095,1);

%% now to determine the percent of VX of FW that each current contributes at
% each latitude (FWP is Fresh Water percentage)

% @ 81.45
egcc_fwp_8145 = (vx_8145_egcc./vx_FW_8145);
egc_fwp_8145 = (vx_8145_egc./vx_FW_8145);
egccc_fwp_8145 = (vx_8145_egccc./vx_FW_8145);

shelf_fwp_8145 = (vx_8145_egcc+vx_8145_egccc)./vx_FW_8145;
% resid_fwp_8145 = (vx_FW_8145-((vx_8145_egc)+(vx_8145_egcc)+(vx_8145_egccc)))./vx_FW_8145;

% the average for all
egccc_fwp_8145avg = nanmean(egccc_fwp_8145);
egcc_fwp_8145avg = nanmean(egcc_fwp_8145);
egc_fwp_8145avg = nanmean(egc_fwp_8145);

shelf_fwp_8145avg = nanmean(shelf_fwp_8145);
% resid_fwp_8145avg = nanmean(resid_fwp_8145);
vx_FW_8145avg = nanmean(vx_FW_8145);

% @ 80.05
egcc_fwp_8005 = (vx_8005_egcc./vx_FW_8005);
egc_fwp_8005 = (vx_8005_egc./vx_FW_8005);
egccc_fwp_8005 = (vx_8005_egccc./vx_FW_8005);

shelf_fwp_8005 = (vx_8005_egcc+vx_8005_egccc)./vx_FW_8005;

% resid_fwp_8005 = (vx_FW_8005-((vx_8005_egc)+(vx_8005_egcc)+(vx_8005_egccc)))./vx_FW_8005;

% the average for all
egccc_fwp_8005avg = nanmean(egccc_fwp_8005);
egcc_fwp_8005avg = nanmean(egcc_fwp_8005);
egc_fwp_8005avg = nanmean(egc_fwp_8005);

shelf_fwp_8005avg = nanmean(shelf_fwp_8005);
% resid_fwp_8005avg = nanmean(resid_fwp_8005);
vx_FW_8005avg = nanmean(vx_FW_8005);

% @ 78.85
egcc_fwp_7885 = (vx_7885_egcc./vx_FW_7885);
egc_fwp_7885 = (vx_7885_egc./vx_FW_7885);
egccc_fwp_7885 = (vx_7885_egccc./vx_FW_7885);

shelf_fwp_7885 = (vx_7885_egcc+vx_7885_egccc)./vx_FW_7885;

% resid_fwp_7885 = (vx_FW_7885-((vx_7885_egc)+(vx_7885_egcc)+(vx_7885_egccc)))./vx_FW_7885;

% the average for all
egccc_fwp_7885avg = nanmean(egccc_fwp_7885);
egcc_fwp_7885avg = nanmean(egcc_fwp_7885);
egc_fwp_7885avg = nanmean(egc_fwp_7885);

shelf_fwp_7885avg = nanmean(shelf_fwp_7885);
% resid_fwp_7885avg = nanmean(resid_fwp_7885);
vx_FW_7885avg = nanmean(vx_FW_7885);

% @ 76.95
egcc_fwp_7695 = (vx_7695_egcc./vx_FW_7695);
egc_fwp_7695 = (vx_7695_egc./vx_FW_7695);
egccc_fwp_7695 = (vx_7695_egccc./vx_FW_7695);

shelf_fwp_7695 = (vx_7695_egcc+vx_7695_egccc)./vx_FW_7695;

% resid_fwp_7695 = (vx_FW_7695-((vx_7695_egc)+(vx_7695_egcc)+(vx_7695_egccc)))./vx_FW_7695;

% the average for all
egccc_fwp_7695avg = nanmean(egccc_fwp_7695);
egcc_fwp_7695avg = nanmean(egcc_fwp_7695);
egc_fwp_7695avg = nanmean(egc_fwp_7695);

shelf_fwp_7695avg = nanmean(shelf_fwp_7695);
% resid_fwp_7695avg = nanmean(resid_fwp_7695);
vx_FW_7695avg = nanmean(vx_FW_7695);

% @ 73.95
egcc_fwp_7395 = (vx_7395_egcc./vx_FW_7395);
egc_fwp_7395 = (vx_7395_egc./vx_FW_7395);
egccc_fwp_7395 = (vx_7395_egccc./vx_FW_7395);

shelf_fwp_7395 = (vx_7395_egcc+vx_7395_egccc)./vx_FW_7395;

% resid_fwp_7395 = (vx_FW_7395-((vx_7395_egc)+(vx_7395_egcc)+(vx_7395_egccc)))./vx_FW_7395;

% the average for all
egccc_fwp_7395avg = nanmean(egccc_fwp_7395);
egcc_fwp_7395avg = nanmean(egcc_fwp_7395);
egc_fwp_7395avg = nanmean(egc_fwp_7395);

shelf_fwp_7395avg = nanmean(shelf_fwp_7395);
% resid_fwp_7395avg = nanmean(resid_fwp_7395);
vx_FW_7395avg = nanmean(vx_FW_7395);

% @ 70.95
egcc_fwp_7095 = (vx_7095_egcc./vx_FW_7095);
egc_fwp_7095 = (vx_7095_egc./vx_FW_7095);
% egccc_fwp_7095 = (vx_7095_egccc./vx_FW_7095);

shelf_fwp_7095 = (vx_7095_egcc)./vx_FW_7095;

% resid_fwp_7095 = (vx_FW_7095-((vx_7095_egc)+(vx_7095_egcc)+(vx_7095_egccc)))./vx_FW_7095;

% the average for all
% egccc_fwp_7095avg = nanmean(egccc_fwp_7095);
egcc_fwp_7095avg = nanmean(egcc_fwp_7095);
egc_fwp_7095avg = nanmean(egc_fwp_7095);

shelf_fwp_7095avg = nanmean(shelf_fwp_7095);
% resid_fwp_7095avg = nanmean(resid_fwp_7095);
vx_FW_7095avg = nanmean(vx_FW_7095);

% now to plot all of these values
lats = [81.45 80.05 78.85 76.95 73.95 70.95];
egcc_fwp_avg = [egcc_fwp_8145avg, egcc_fwp_8005avg, egcc_fwp_7885avg, egcc_fwp_7695avg, egcc_fwp_7395avg, egcc_fwp_7095avg];
egc_fwp_avg = [egc_fwp_8145avg, egc_fwp_8005avg, egc_fwp_7885avg, egc_fwp_7695avg, egc_fwp_7395avg, egc_fwp_7095avg];
% egccc_fwp_avg = [egccc_fwp_8145avg, egccc_fwp_8005avg, egccc_fwp_7885avg, egccc_fwp_7695avg, egccc_fwp_7395avg, egccc_fwp_7095avg];
shelf_fwp_avg = [shelf_fwp_8145avg, shelf_fwp_8005avg, shelf_fwp_7885avg, shelf_fwp_7695avg, shelf_fwp_7395avg, shelf_fwp_7095avg];
% resid_fwp_avg = [resid_fwp_8145avg, resid_fwp_8005avg, resid_fwp_7885avg, resid_fwp_7695avg, resid_fwp_7395avg, resid_fwp_7095avg];


%%
figure
hold on
plot(egc_fwp_avg,[1:6],'-*','Color',[0.7 0 0.4],'LineWidth',1,'MarkerSize',7,'MarkerFaceColor',[0.7 0 0.4],'MarkerEdgeColor',[0.7 0 0.4]);
plot(shelf_fwp_avg,[1:6],'-pentagram','Color',[0.5 0.5 1],'LineWidth',1,'MarkerSize',7,'MarkerFaceColor',[0.5 0.5 1],'MarkerEdgeColor',[0.5 0.5 1]);
% plot(resid_fwp_avg,[1:6],'-square','Color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',7,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);

yticks(1:6);
ylim([0.5 6.5])
yticklabels(lats)
set(gca,'YDir','reverse')
xlim([0 1])

legend('EGC','Shelf')
xlabel('Percentage')
ylabel('Latitude [^oN]')
title('Average Percentage of Fresh Water Volume Transport (RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/VxFw_percent_avg.jpg']);

%% Seasonality of FW VX

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

%% @ 81.45

% establishing overall averages
vx_8145_egccc_sea = [nanmean(vx_8145_egccc(jan_i)), nanmean(vx_8145_egccc(feb_i)),...
    nanmean(vx_8145_egccc(mar_i)), nanmean(vx_8145_egccc(apr_i)), nanmean(vx_8145_egccc(may_i)),...
    nanmean(vx_8145_egccc(jun_i)), nanmean(vx_8145_egccc(jul_i)), nanmean(vx_8145_egccc(aug_i)),...
    nanmean(vx_8145_egccc(sep_i)), nanmean(vx_8145_egccc(oct_i)), nanmean(vx_8145_egccc(nov_i)),...
    nanmean(vx_8145_egccc(dec_i))];

vx_8145_egcc_sea = [nanmean(vx_8145_egcc(jan_i)), nanmean(vx_8145_egcc(feb_i)),...
    nanmean(vx_8145_egcc(mar_i)), nanmean(vx_8145_egcc(apr_i)), nanmean(vx_8145_egcc(may_i)),...
    nanmean(vx_8145_egcc(jun_i)), nanmean(vx_8145_egcc(jul_i)), nanmean(vx_8145_egcc(aug_i)),...
    nanmean(vx_8145_egcc(sep_i)), nanmean(vx_8145_egcc(oct_i)), nanmean(vx_8145_egcc(nov_i)),...
    nanmean(vx_8145_egcc(dec_i))];

vx_8145_egc_sea = [nanmean(vx_8145_egc(jan_i)), nanmean(vx_8145_egc(feb_i)),...
    nanmean(vx_8145_egc(mar_i)), nanmean(vx_8145_egc(apr_i)), nanmean(vx_8145_egc(may_i)),...
    nanmean(vx_8145_egc(jun_i)), nanmean(vx_8145_egc(jul_i)), nanmean(vx_8145_egc(aug_i)),...
    nanmean(vx_8145_egc(sep_i)), nanmean(vx_8145_egc(oct_i)), nanmean(vx_8145_egc(nov_i)),...
    nanmean(vx_8145_egc(dec_i))];

vx_8145_sea = [nanmean(vx_FW_8145(jan_i)), nanmean(vx_FW_8145(feb_i)),...
    nanmean(vx_FW_8145(mar_i)), nanmean(vx_FW_8145(apr_i)), nanmean(vx_FW_8145(may_i)),...
    nanmean(vx_FW_8145(jun_i)), nanmean(vx_FW_8145(jul_i)), nanmean(vx_FW_8145(aug_i)),...
    nanmean(vx_FW_8145(sep_i)), nanmean(vx_FW_8145(oct_i)), nanmean(vx_FW_8145(nov_i)),...
    nanmean(vx_FW_8145(dec_i))];

% now calculating the variation in the seasonal cycle

vx_8145_egccc_mat = [(vx_8145_egccc(jan_i)); (vx_8145_egccc(feb_i));...
    (vx_8145_egccc(mar_i)); (vx_8145_egccc(apr_i)); (vx_8145_egccc(may_i));...
    (vx_8145_egccc(jun_i)); (vx_8145_egccc(jul_i)); (vx_8145_egccc(aug_i));...
    (vx_8145_egccc(sep_i)); (vx_8145_egccc(oct_i)); (vx_8145_egccc(nov_i));...
    (vx_8145_egccc(dec_i))];

vx_8145_egcc_mat = [(vx_8145_egcc(jan_i)); (vx_8145_egcc(feb_i));...
    (vx_8145_egcc(mar_i)); (vx_8145_egcc(apr_i)); (vx_8145_egcc(may_i));...
    (vx_8145_egcc(jun_i)); (vx_8145_egcc(jul_i)); (vx_8145_egcc(aug_i));...
    (vx_8145_egcc(sep_i)); (vx_8145_egcc(oct_i)); (vx_8145_egcc(nov_i));...
    (vx_8145_egcc(dec_i))];

vx_8145_egc_mat = [(vx_8145_egc(jan_i)); (vx_8145_egc(feb_i));...
    (vx_8145_egc(mar_i)); (vx_8145_egc(apr_i)); (vx_8145_egc(may_i));...
    (vx_8145_egc(jun_i)); (vx_8145_egc(jul_i)); (vx_8145_egc(aug_i));...
    (vx_8145_egc(sep_i)); (vx_8145_egc(oct_i)); (vx_8145_egc(nov_i));...
    (vx_8145_egc(dec_i))];

vx_8145_mat = [(vx_FW_8145(jan_i)); (vx_FW_8145(feb_i));...
    (vx_FW_8145(mar_i)); (vx_FW_8145(apr_i)); (vx_FW_8145(may_i));...
    (vx_FW_8145(jun_i)); (vx_FW_8145(jul_i)); (vx_FW_8145(aug_i));...
    (vx_FW_8145(sep_i)); (vx_FW_8145(oct_i)); (vx_FW_8145(nov_i));...
    (vx_FW_8145(dec_i))];

% calculating the SD of each of the months for the 3 currents

for jj = 1:12
vx_8145_egc_SD(jj,:) = std(vx_8145_egc_mat(jj,:),'omitnan');
vx_8145_egcc_SD(jj,:) = std(vx_8145_egcc_mat(jj,:),'omitnan');
vx_8145_egccc_SD(jj,:) = std(vx_8145_egccc_mat(jj,:),'omitnan');
vx_8145_SD(jj,:) = std(vx_8145_mat(jj,:),'omitnan');
end


figure
hold on
shadedErrorBar([1:12],vx_8145_egccc_sea,vx_8145_egccc_SD,'lineprops','c')
shadedErrorBar([1:12],vx_8145_egcc_sea,vx_8145_egcc_SD,'lineprops','b')
shadedErrorBar([1:12],vx_8145_egc_sea,vx_8145_egc_SD,'lineprops','r')
shadedErrorBar([1:12],vx_8145_sea,vx_8145_SD,'lineprops','k')

plot([1:12],vx_8145_sea,'-k','LineWidth',2)
plot([1:12],vx_8145_egccc_sea,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:12],vx_8145_egcc_sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],vx_8145_egc_sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([-2 1])

legend('','','','','Total','EGCCC','EGCC','EGC','location','southeast')
ylabel('Volume Transport [Sv]')
title('Seasonal Cycle of Fresh Water Volume Transport at 81.45 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/8145.jpg']);

%% @ 80.05
% establishing overall averages
vx_8005_egccc_sea = [nanmean(vx_8005_egccc(jan_i)), nanmean(vx_8005_egccc(feb_i)),...
    nanmean(vx_8005_egccc(mar_i)), nanmean(vx_8005_egccc(apr_i)), nanmean(vx_8005_egccc(may_i)),...
    nanmean(vx_8005_egccc(jun_i)), nanmean(vx_8005_egccc(jul_i)), nanmean(vx_8005_egccc(aug_i)),...
    nanmean(vx_8005_egccc(sep_i)), nanmean(vx_8005_egccc(oct_i)), nanmean(vx_8005_egccc(nov_i)),...
    nanmean(vx_8005_egccc(dec_i))];

vx_8005_egcc_sea = [nanmean(vx_8005_egcc(jan_i)), nanmean(vx_8005_egcc(feb_i)),...
    nanmean(vx_8005_egcc(mar_i)), nanmean(vx_8005_egcc(apr_i)), nanmean(vx_8005_egcc(may_i)),...
    nanmean(vx_8005_egcc(jun_i)), nanmean(vx_8005_egcc(jul_i)), nanmean(vx_8005_egcc(aug_i)),...
    nanmean(vx_8005_egcc(sep_i)), nanmean(vx_8005_egcc(oct_i)), nanmean(vx_8005_egcc(nov_i)),...
    nanmean(vx_8005_egcc(dec_i))];

vx_8005_egc_sea = [nanmean(vx_8005_egc(jan_i)), nanmean(vx_8005_egc(feb_i)),...
    nanmean(vx_8005_egc(mar_i)), nanmean(vx_8005_egc(apr_i)), nanmean(vx_8005_egc(may_i)),...
    nanmean(vx_8005_egc(jun_i)), nanmean(vx_8005_egc(jul_i)), nanmean(vx_8005_egc(aug_i)),...
    nanmean(vx_8005_egc(sep_i)), nanmean(vx_8005_egc(oct_i)), nanmean(vx_8005_egc(nov_i)),...
    nanmean(vx_8005_egc(dec_i))];

vx_8005_sea = [nanmean(vx_FW_8005(jan_i)), nanmean(vx_FW_8005(feb_i)),...
    nanmean(vx_FW_8005(mar_i)), nanmean(vx_FW_8005(apr_i)), nanmean(vx_FW_8005(may_i)),...
    nanmean(vx_FW_8005(jun_i)), nanmean(vx_FW_8005(jul_i)), nanmean(vx_FW_8005(aug_i)),...
    nanmean(vx_FW_8005(sep_i)), nanmean(vx_FW_8005(oct_i)), nanmean(vx_FW_8005(nov_i)),...
    nanmean(vx_FW_8005(dec_i))];

% now calculating the variation in the seasonal cycle

vx_8005_egccc_mat = [(vx_8005_egccc(jan_i)); (vx_8005_egccc(feb_i));...
    (vx_8005_egccc(mar_i)); (vx_8005_egccc(apr_i)); (vx_8005_egccc(may_i));...
    (vx_8005_egccc(jun_i)); (vx_8005_egccc(jul_i)); (vx_8005_egccc(aug_i));...
    (vx_8005_egccc(sep_i)); (vx_8005_egccc(oct_i)); (vx_8005_egccc(nov_i));...
    (vx_8005_egccc(dec_i))];

vx_8005_egcc_mat = [(vx_8005_egcc(jan_i)); (vx_8005_egcc(feb_i));...
    (vx_8005_egcc(mar_i)); (vx_8005_egcc(apr_i)); (vx_8005_egcc(may_i));...
    (vx_8005_egcc(jun_i)); (vx_8005_egcc(jul_i)); (vx_8005_egcc(aug_i));...
    (vx_8005_egcc(sep_i)); (vx_8005_egcc(oct_i)); (vx_8005_egcc(nov_i));...
    (vx_8005_egcc(dec_i))];

vx_8005_egc_mat = [(vx_8005_egc(jan_i)); (vx_8005_egc(feb_i));...
    (vx_8005_egc(mar_i)); (vx_8005_egc(apr_i)); (vx_8005_egc(may_i));...
    (vx_8005_egc(jun_i)); (vx_8005_egc(jul_i)); (vx_8005_egc(aug_i));...
    (vx_8005_egc(sep_i)); (vx_8005_egc(oct_i)); (vx_8005_egc(nov_i));...
    (vx_8005_egc(dec_i))];

vx_8005_mat = [(vx_FW_8005(jan_i)); (vx_FW_8005(feb_i));...
    (vx_FW_8005(mar_i)); (vx_FW_8005(apr_i)); (vx_FW_8005(may_i));...
    (vx_FW_8005(jun_i)); (vx_FW_8005(jul_i)); (vx_FW_8005(aug_i));...
    (vx_FW_8005(sep_i)); (vx_FW_8005(oct_i)); (vx_FW_8005(nov_i));...
    (vx_FW_8005(dec_i))];

% calculating the SD of each of the months for the 3 currents

for jj = 1:12
vx_8005_egc_SD(jj,:) = std(vx_8005_egc_mat(jj,:),'omitnan');
vx_8005_egcc_SD(jj,:) = std(vx_8005_egcc_mat(jj,:),'omitnan');
vx_8005_egccc_SD(jj,:) = std(vx_8005_egccc_mat(jj,:),'omitnan');
vx_8005_SD(jj,:) = std(vx_8005_mat(jj,:),'omitnan');
end


figure
hold on
shadedErrorBar([1:12],vx_8005_egccc_sea,vx_8005_egccc_SD,'lineprops','c')
shadedErrorBar([1:12],vx_8005_egcc_sea,vx_8005_egcc_SD,'lineprops','b')
shadedErrorBar([1:12],vx_8005_egc_sea,vx_8005_egc_SD,'lineprops','r')
shadedErrorBar([1:12],vx_8005_sea,vx_8005_SD,'lineprops','k')

plot([1:12],vx_8005_sea,'-k','LineWidth',2)
plot([1:12],vx_8005_egccc_sea,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:12],vx_8005_egcc_sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],vx_8005_egc_sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([-2 1])

legend('','','','','Total','EGCCC','EGCC','EGC','location','southeast')
ylabel('Volume Transport [Sv]')
title('Seasonal Cycle of Fresh Water Volume Transport at 80.05 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/8005.jpg']);

%% @ 78.85

% establishing overall averages
vx_7885_egccc_sea = [nanmean(vx_7885_egccc(jan_i)), nanmean(vx_7885_egccc(feb_i)),...
    nanmean(vx_7885_egccc(mar_i)), nanmean(vx_7885_egccc(apr_i)), nanmean(vx_7885_egccc(may_i)),...
    nanmean(vx_7885_egccc(jun_i)), nanmean(vx_7885_egccc(jul_i)), nanmean(vx_7885_egccc(aug_i)),...
    nanmean(vx_7885_egccc(sep_i)), nanmean(vx_7885_egccc(oct_i)), nanmean(vx_7885_egccc(nov_i)),...
    nanmean(vx_7885_egccc(dec_i))];

vx_7885_egcc_sea = [nanmean(vx_7885_egcc(jan_i)), nanmean(vx_7885_egcc(feb_i)),...
    nanmean(vx_7885_egcc(mar_i)), nanmean(vx_7885_egcc(apr_i)), nanmean(vx_7885_egcc(may_i)),...
    nanmean(vx_7885_egcc(jun_i)), nanmean(vx_7885_egcc(jul_i)), nanmean(vx_7885_egcc(aug_i)),...
    nanmean(vx_7885_egcc(sep_i)), nanmean(vx_7885_egcc(oct_i)), nanmean(vx_7885_egcc(nov_i)),...
    nanmean(vx_7885_egcc(dec_i))];

vx_7885_egc_sea = [nanmean(vx_7885_egc(jan_i)), nanmean(vx_7885_egc(feb_i)),...
    nanmean(vx_7885_egc(mar_i)), nanmean(vx_7885_egc(apr_i)), nanmean(vx_7885_egc(may_i)),...
    nanmean(vx_7885_egc(jun_i)), nanmean(vx_7885_egc(jul_i)), nanmean(vx_7885_egc(aug_i)),...
    nanmean(vx_7885_egc(sep_i)), nanmean(vx_7885_egc(oct_i)), nanmean(vx_7885_egc(nov_i)),...
    nanmean(vx_7885_egc(dec_i))];

vx_7885_sea = [nanmean(vx_FW_7885(jan_i)), nanmean(vx_FW_7885(feb_i)),...
    nanmean(vx_FW_7885(mar_i)), nanmean(vx_FW_7885(apr_i)), nanmean(vx_FW_7885(may_i)),...
    nanmean(vx_FW_7885(jun_i)), nanmean(vx_FW_7885(jul_i)), nanmean(vx_FW_7885(aug_i)),...
    nanmean(vx_FW_7885(sep_i)), nanmean(vx_FW_7885(oct_i)), nanmean(vx_FW_7885(nov_i)),...
    nanmean(vx_FW_7885(dec_i))];

% now calculating the variation in the seasonal cycle

vx_7885_egccc_mat = [(vx_7885_egccc(jan_i)); (vx_7885_egccc(feb_i));...
    (vx_7885_egccc(mar_i)); (vx_7885_egccc(apr_i)); (vx_7885_egccc(may_i));...
    (vx_7885_egccc(jun_i)); (vx_7885_egccc(jul_i)); (vx_7885_egccc(aug_i));...
    (vx_7885_egccc(sep_i)); (vx_7885_egccc(oct_i)); (vx_7885_egccc(nov_i));...
    (vx_7885_egccc(dec_i))];

vx_7885_egcc_mat = [(vx_7885_egcc(jan_i)); (vx_7885_egcc(feb_i));...
    (vx_7885_egcc(mar_i)); (vx_7885_egcc(apr_i)); (vx_7885_egcc(may_i));...
    (vx_7885_egcc(jun_i)); (vx_7885_egcc(jul_i)); (vx_7885_egcc(aug_i));...
    (vx_7885_egcc(sep_i)); (vx_7885_egcc(oct_i)); (vx_7885_egcc(nov_i));...
    (vx_7885_egcc(dec_i))];

vx_7885_egc_mat = [(vx_7885_egc(jan_i)); (vx_7885_egc(feb_i));...
    (vx_7885_egc(mar_i)); (vx_7885_egc(apr_i)); (vx_7885_egc(may_i));...
    (vx_7885_egc(jun_i)); (vx_7885_egc(jul_i)); (vx_7885_egc(aug_i));...
    (vx_7885_egc(sep_i)); (vx_7885_egc(oct_i)); (vx_7885_egc(nov_i));...
    (vx_7885_egc(dec_i))];

vx_7885_mat = [(vx_FW_7885(jan_i)); (vx_FW_7885(feb_i));...
    (vx_FW_7885(mar_i)); (vx_FW_7885(apr_i)); (vx_FW_7885(may_i));...
    (vx_FW_7885(jun_i)); (vx_FW_7885(jul_i)); (vx_FW_7885(aug_i));...
    (vx_FW_7885(sep_i)); (vx_FW_7885(oct_i)); (vx_FW_7885(nov_i));...
    (vx_FW_7885(dec_i))];

% calculating the SD of each of the months for the 3 currents

for jj = 1:12
vx_7885_egc_SD(jj,:) = std(vx_7885_egc_mat(jj,:),'omitnan');
vx_7885_egcc_SD(jj,:) = std(vx_7885_egcc_mat(jj,:),'omitnan');
vx_7885_egccc_SD(jj,:) = std(vx_7885_egccc_mat(jj,:),'omitnan');
vx_7885_SD(jj,:) = std(vx_7885_mat(jj,:),'omitnan');
end


figure
hold on
shadedErrorBar([1:12],vx_7885_egccc_sea,vx_7885_egccc_SD,'lineprops','c')
shadedErrorBar([1:12],vx_7885_egcc_sea,vx_7885_egcc_SD,'lineprops','b')
shadedErrorBar([1:12],vx_7885_egc_sea,vx_7885_egc_SD,'lineprops','r')
shadedErrorBar([1:12],vx_7885_sea,vx_7885_SD,'lineprops','k')

plot([1:12],vx_7885_sea,'-k','LineWidth',2)
plot([1:12],vx_7885_egccc_sea,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:12],vx_7885_egcc_sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],vx_7885_egc_sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([-2 1])

legend('','','','','Total','EGCCC','EGCC','EGC','location','southeast')
ylabel('Volume Transport [Sv]')
title('Seasonal Cycle of Fresh Water Volume Transport at 78.85 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/7885.jpg']);

%% @ 76.95

% establishing overall averages
vx_7695_egccc_sea = [nanmean(vx_7695_egccc(jan_i)), nanmean(vx_7695_egccc(feb_i)),...
    nanmean(vx_7695_egccc(mar_i)), nanmean(vx_7695_egccc(apr_i)), nanmean(vx_7695_egccc(may_i)),...
    nanmean(vx_7695_egccc(jun_i)), nanmean(vx_7695_egccc(jul_i)), nanmean(vx_7695_egccc(aug_i)),...
    nanmean(vx_7695_egccc(sep_i)), nanmean(vx_7695_egccc(oct_i)), nanmean(vx_7695_egccc(nov_i)),...
    nanmean(vx_7695_egccc(dec_i))];

vx_7695_egcc_sea = [nanmean(vx_7695_egcc(jan_i)), nanmean(vx_7695_egcc(feb_i)),...
    nanmean(vx_7695_egcc(mar_i)), nanmean(vx_7695_egcc(apr_i)), nanmean(vx_7695_egcc(may_i)),...
    nanmean(vx_7695_egcc(jun_i)), nanmean(vx_7695_egcc(jul_i)), nanmean(vx_7695_egcc(aug_i)),...
    nanmean(vx_7695_egcc(sep_i)), nanmean(vx_7695_egcc(oct_i)), nanmean(vx_7695_egcc(nov_i)),...
    nanmean(vx_7695_egcc(dec_i))];

vx_7695_egc_sea = [nanmean(vx_7695_egc(jan_i)), nanmean(vx_7695_egc(feb_i)),...
    nanmean(vx_7695_egc(mar_i)), nanmean(vx_7695_egc(apr_i)), nanmean(vx_7695_egc(may_i)),...
    nanmean(vx_7695_egc(jun_i)), nanmean(vx_7695_egc(jul_i)), nanmean(vx_7695_egc(aug_i)),...
    nanmean(vx_7695_egc(sep_i)), nanmean(vx_7695_egc(oct_i)), nanmean(vx_7695_egc(nov_i)),...
    nanmean(vx_7695_egc(dec_i))];

vx_7695_sea = [nanmean(vx_FW_7695(jan_i)), nanmean(vx_FW_7695(feb_i)),...
    nanmean(vx_FW_7695(mar_i)), nanmean(vx_FW_7695(apr_i)), nanmean(vx_FW_7695(may_i)),...
    nanmean(vx_FW_7695(jun_i)), nanmean(vx_FW_7695(jul_i)), nanmean(vx_FW_7695(aug_i)),...
    nanmean(vx_FW_7695(sep_i)), nanmean(vx_FW_7695(oct_i)), nanmean(vx_FW_7695(nov_i)),...
    nanmean(vx_FW_7695(dec_i))];

% now calculating the variation in the seasonal cycle

vx_7695_egccc_mat = [(vx_7695_egccc(jan_i)); (vx_7695_egccc(feb_i));...
    (vx_7695_egccc(mar_i)); (vx_7695_egccc(apr_i)); (vx_7695_egccc(may_i));...
    (vx_7695_egccc(jun_i)); (vx_7695_egccc(jul_i)); (vx_7695_egccc(aug_i));...
    (vx_7695_egccc(sep_i)); (vx_7695_egccc(oct_i)); (vx_7695_egccc(nov_i));...
    (vx_7695_egccc(dec_i))];

vx_7695_egcc_mat = [(vx_7695_egcc(jan_i)); (vx_7695_egcc(feb_i));...
    (vx_7695_egcc(mar_i)); (vx_7695_egcc(apr_i)); (vx_7695_egcc(may_i));...
    (vx_7695_egcc(jun_i)); (vx_7695_egcc(jul_i)); (vx_7695_egcc(aug_i));...
    (vx_7695_egcc(sep_i)); (vx_7695_egcc(oct_i)); (vx_7695_egcc(nov_i));...
    (vx_7695_egcc(dec_i))];

vx_7695_egc_mat = [(vx_7695_egc(jan_i)); (vx_7695_egc(feb_i));...
    (vx_7695_egc(mar_i)); (vx_7695_egc(apr_i)); (vx_7695_egc(may_i));...
    (vx_7695_egc(jun_i)); (vx_7695_egc(jul_i)); (vx_7695_egc(aug_i));...
    (vx_7695_egc(sep_i)); (vx_7695_egc(oct_i)); (vx_7695_egc(nov_i));...
    (vx_7695_egc(dec_i))];

vx_7695_mat = [(vx_FW_7695(jan_i)); (vx_FW_7695(feb_i));...
    (vx_FW_7695(mar_i)); (vx_FW_7695(apr_i)); (vx_FW_7695(may_i));...
    (vx_FW_7695(jun_i)); (vx_FW_7695(jul_i)); (vx_FW_7695(aug_i));...
    (vx_FW_7695(sep_i)); (vx_FW_7695(oct_i)); (vx_FW_7695(nov_i));...
    (vx_FW_7695(dec_i))];

% calculating the SD of each of the months for the 3 currents

for jj = 1:12
vx_7695_egc_SD(jj,:) = std(vx_7695_egc_mat(jj,:),'omitnan');
vx_7695_egcc_SD(jj,:) = std(vx_7695_egcc_mat(jj,:),'omitnan');
vx_7695_egccc_SD(jj,:) = std(vx_7695_egccc_mat(jj,:),'omitnan');
vx_7695_SD(jj,:) = std(vx_7695_mat(jj,:),'omitnan');
end


figure
hold on
shadedErrorBar([1:12],vx_7695_egccc_sea,vx_7695_egccc_SD,'lineprops','c')
shadedErrorBar([1:12],vx_7695_egcc_sea,vx_7695_egcc_SD,'lineprops','b')
shadedErrorBar([1:12],vx_7695_egc_sea,vx_7695_egc_SD,'lineprops','r')
shadedErrorBar([1:12],vx_7695_sea,vx_7695_SD,'lineprops','k')

plot([1:12],vx_7695_sea,'-k','LineWidth',2)
plot([1:12],vx_7695_egccc_sea,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:12],vx_7695_egcc_sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],vx_7695_egc_sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([-2 1])

legend('','','','','Total','EGCCC','EGCC','EGC','location','southeast')
ylabel('Volume Transport [Sv]')
title('Seasonal Cycle of Fresh Water Volume Transport at 76.95 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/7695.jpg']);

%% @ 73.95

% establishing overall averages
vx_7395_egccc_sea = [nanmean(vx_7395_egccc(jan_i)), nanmean(vx_7395_egccc(feb_i)),...
    nanmean(vx_7395_egccc(mar_i)), nanmean(vx_7395_egccc(apr_i)), nanmean(vx_7395_egccc(may_i)),...
    nanmean(vx_7395_egccc(jun_i)), nanmean(vx_7395_egccc(jul_i)), nanmean(vx_7395_egccc(aug_i)),...
    nanmean(vx_7395_egccc(sep_i)), nanmean(vx_7395_egccc(oct_i)), nanmean(vx_7395_egccc(nov_i)),...
    nanmean(vx_7395_egccc(dec_i))];

vx_7395_egcc_sea = [nanmean(vx_7395_egcc(jan_i)), nanmean(vx_7395_egcc(feb_i)),...
    nanmean(vx_7395_egcc(mar_i)), nanmean(vx_7395_egcc(apr_i)), nanmean(vx_7395_egcc(may_i)),...
    nanmean(vx_7395_egcc(jun_i)), nanmean(vx_7395_egcc(jul_i)), nanmean(vx_7395_egcc(aug_i)),...
    nanmean(vx_7395_egcc(sep_i)), nanmean(vx_7395_egcc(oct_i)), nanmean(vx_7395_egcc(nov_i)),...
    nanmean(vx_7395_egcc(dec_i))];

vx_7395_egc_sea = [nanmean(vx_7395_egc(jan_i)), nanmean(vx_7395_egc(feb_i)),...
    nanmean(vx_7395_egc(mar_i)), nanmean(vx_7395_egc(apr_i)), nanmean(vx_7395_egc(may_i)),...
    nanmean(vx_7395_egc(jun_i)), nanmean(vx_7395_egc(jul_i)), nanmean(vx_7395_egc(aug_i)),...
    nanmean(vx_7395_egc(sep_i)), nanmean(vx_7395_egc(oct_i)), nanmean(vx_7395_egc(nov_i)),...
    nanmean(vx_7395_egc(dec_i))];

vx_7395_sea = [nanmean(vx_FW_7395(jan_i)), nanmean(vx_FW_7395(feb_i)),...
    nanmean(vx_FW_7395(mar_i)), nanmean(vx_FW_7395(apr_i)), nanmean(vx_FW_7395(may_i)),...
    nanmean(vx_FW_7395(jun_i)), nanmean(vx_FW_7395(jul_i)), nanmean(vx_FW_7395(aug_i)),...
    nanmean(vx_FW_7395(sep_i)), nanmean(vx_FW_7395(oct_i)), nanmean(vx_FW_7395(nov_i)),...
    nanmean(vx_FW_7395(dec_i))];

% now calculating the variation in the seasonal cycle

vx_7395_egccc_mat = [(vx_7395_egccc(jan_i)); (vx_7395_egccc(feb_i));...
    (vx_7395_egccc(mar_i)); (vx_7395_egccc(apr_i)); (vx_7395_egccc(may_i));...
    (vx_7395_egccc(jun_i)); (vx_7395_egccc(jul_i)); (vx_7395_egccc(aug_i));...
    (vx_7395_egccc(sep_i)); (vx_7395_egccc(oct_i)); (vx_7395_egccc(nov_i));...
    (vx_7395_egccc(dec_i))];

vx_7395_egcc_mat = [(vx_7395_egcc(jan_i)); (vx_7395_egcc(feb_i));...
    (vx_7395_egcc(mar_i)); (vx_7395_egcc(apr_i)); (vx_7395_egcc(may_i));...
    (vx_7395_egcc(jun_i)); (vx_7395_egcc(jul_i)); (vx_7395_egcc(aug_i));...
    (vx_7395_egcc(sep_i)); (vx_7395_egcc(oct_i)); (vx_7395_egcc(nov_i));...
    (vx_7395_egcc(dec_i))];

vx_7395_egc_mat = [(vx_7395_egc(jan_i)); (vx_7395_egc(feb_i));...
    (vx_7395_egc(mar_i)); (vx_7395_egc(apr_i)); (vx_7395_egc(may_i));...
    (vx_7395_egc(jun_i)); (vx_7395_egc(jul_i)); (vx_7395_egc(aug_i));...
    (vx_7395_egc(sep_i)); (vx_7395_egc(oct_i)); (vx_7395_egc(nov_i));...
    (vx_7395_egc(dec_i))];

vx_7395_mat = [(vx_FW_7395(jan_i)); (vx_FW_7395(feb_i));...
    (vx_FW_7395(mar_i)); (vx_FW_7395(apr_i)); (vx_FW_7395(may_i));...
    (vx_FW_7395(jun_i)); (vx_FW_7395(jul_i)); (vx_FW_7395(aug_i));...
    (vx_FW_7395(sep_i)); (vx_FW_7395(oct_i)); (vx_FW_7395(nov_i));...
    (vx_FW_7395(dec_i))];

% calculating the SD of each of the months for the 3 currents

for jj = 1:12
vx_7395_egc_SD(jj,:) = std(vx_7395_egc_mat(jj,:),'omitnan');
vx_7395_egcc_SD(jj,:) = std(vx_7395_egcc_mat(jj,:),'omitnan');
vx_7395_egccc_SD(jj,:) = std(vx_7395_egccc_mat(jj,:),'omitnan');
vx_7395_SD(jj,:) = std(vx_7395_mat(jj,:),'omitnan');
end


figure
hold on
shadedErrorBar([1:12],vx_7395_egccc_sea,vx_7395_egccc_SD,'lineprops','c')
shadedErrorBar([1:12],vx_7395_egcc_sea,vx_7395_egcc_SD,'lineprops','b')
shadedErrorBar([1:12],vx_7395_egc_sea,vx_7395_egc_SD,'lineprops','r')
shadedErrorBar([1:12],vx_7395_sea,vx_7395_SD,'lineprops','k')

plot([1:12],vx_7395_sea,'-k','LineWidth',2)
plot([1:12],vx_7395_egccc_sea,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:12],vx_7395_egcc_sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],vx_7395_egc_sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([-2 1])

legend('','','','','Total','EGCCC','EGCC','EGC','location','southeast')
ylabel('Volume Transport [Sv]')
title('Seasonal Cycle of Fresh Water Volume Transport at 73.95 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/7395.jpg']);

%% @ 70.95

% establishing overall averages
% vx_7095_egccc_sea = [nanmean(vx_7095_egccc(jan_i)), nanmean(vx_7095_egccc(feb_i)),...
%     nanmean(vx_7095_egccc(mar_i)), nanmean(vx_7095_egccc(apr_i)), nanmean(vx_7095_egccc(may_i)),...
%     nanmean(vx_7095_egccc(jun_i)), nanmean(vx_7095_egccc(jul_i)), nanmean(vx_7095_egccc(aug_i)),...
%     nanmean(vx_7095_egccc(sep_i)), nanmean(vx_7095_egccc(oct_i)), nanmean(vx_7095_egccc(nov_i)),...
%     nanmean(vx_7095_egccc(dec_i))];

vx_7095_egcc_sea = [nanmean(vx_7095_egcc(jan_i)), nanmean(vx_7095_egcc(feb_i)),...
    nanmean(vx_7095_egcc(mar_i)), nanmean(vx_7095_egcc(apr_i)), nanmean(vx_7095_egcc(may_i)),...
    nanmean(vx_7095_egcc(jun_i)), nanmean(vx_7095_egcc(jul_i)), nanmean(vx_7095_egcc(aug_i)),...
    nanmean(vx_7095_egcc(sep_i)), nanmean(vx_7095_egcc(oct_i)), nanmean(vx_7095_egcc(nov_i)),...
    nanmean(vx_7095_egcc(dec_i))];

vx_7095_egc_sea = [nanmean(vx_7095_egc(jan_i)), nanmean(vx_7095_egc(feb_i)),...
    nanmean(vx_7095_egc(mar_i)), nanmean(vx_7095_egc(apr_i)), nanmean(vx_7095_egc(may_i)),...
    nanmean(vx_7095_egc(jun_i)), nanmean(vx_7095_egc(jul_i)), nanmean(vx_7095_egc(aug_i)),...
    nanmean(vx_7095_egc(sep_i)), nanmean(vx_7095_egc(oct_i)), nanmean(vx_7095_egc(nov_i)),...
    nanmean(vx_7095_egc(dec_i))];

vx_7095_sea = [nanmean(vx_FW_7095(jan_i)), nanmean(vx_FW_7095(feb_i)),...
    nanmean(vx_FW_7095(mar_i)), nanmean(vx_FW_7095(apr_i)), nanmean(vx_FW_7095(may_i)),...
    nanmean(vx_FW_7095(jun_i)), nanmean(vx_FW_7095(jul_i)), nanmean(vx_FW_7095(aug_i)),...
    nanmean(vx_FW_7095(sep_i)), nanmean(vx_FW_7095(oct_i)), nanmean(vx_FW_7095(nov_i)),...
    nanmean(vx_FW_7095(dec_i))];

% now calculating the variation in the seasonal cycle

% vx_7095_egccc_mat = [(vx_7095_egccc(jan_i)); (vx_7095_egccc(feb_i));...
%     (vx_7095_egccc(mar_i)); (vx_7095_egccc(apr_i)); (vx_7095_egccc(may_i));...
%     (vx_7095_egccc(jun_i)); (vx_7095_egccc(jul_i)); (vx_7095_egccc(aug_i));...
%     (vx_7095_egccc(sep_i)); (vx_7095_egccc(oct_i)); (vx_7095_egccc(nov_i));...
%     (vx_7095_egccc(dec_i))];

vx_7095_egcc_mat = [(vx_7095_egcc(jan_i)); (vx_7095_egcc(feb_i));...
    (vx_7095_egcc(mar_i)); (vx_7095_egcc(apr_i)); (vx_7095_egcc(may_i));...
    (vx_7095_egcc(jun_i)); (vx_7095_egcc(jul_i)); (vx_7095_egcc(aug_i));...
    (vx_7095_egcc(sep_i)); (vx_7095_egcc(oct_i)); (vx_7095_egcc(nov_i));...
    (vx_7095_egcc(dec_i))];

vx_7095_egc_mat = [(vx_7095_egc(jan_i)); (vx_7095_egc(feb_i));...
    (vx_7095_egc(mar_i)); (vx_7095_egc(apr_i)); (vx_7095_egc(may_i));...
    (vx_7095_egc(jun_i)); (vx_7095_egc(jul_i)); (vx_7095_egc(aug_i));...
    (vx_7095_egc(sep_i)); (vx_7095_egc(oct_i)); (vx_7095_egc(nov_i));...
    (vx_7095_egc(dec_i))];

vx_7095_mat = [(vx_FW_7095(jan_i)); (vx_FW_7095(feb_i));...
    (vx_FW_7095(mar_i)); (vx_FW_7095(apr_i)); (vx_FW_7095(may_i));...
    (vx_FW_7095(jun_i)); (vx_FW_7095(jul_i)); (vx_FW_7095(aug_i));...
    (vx_FW_7095(sep_i)); (vx_FW_7095(oct_i)); (vx_FW_7095(nov_i));...
    (vx_FW_7095(dec_i))];

% calculating the SD of each of the months for the 3 currents

for jj = 1:12
vx_7095_egc_SD(jj,:) = std(vx_7095_egc_mat(jj,:),'omitnan');
vx_7095_egcc_SD(jj,:) = std(vx_7095_egcc_mat(jj,:),'omitnan');
% vx_7095_egccc_SD(jj,:) = std(vx_7095_egccc_mat(jj,:),'omitnan');
vx_7095_SD(jj,:) = std(vx_7095_mat(jj,:),'omitnan');
end


figure
hold on
% shadedErrorBar([1:12],vx_7095_egccc_sea,vx_7095_egccc_SD,'lineprops','c')
shadedErrorBar([1:12],vx_7095_egcc_sea,vx_7095_egcc_SD,'lineprops','b')
shadedErrorBar([1:12],vx_7095_egc_sea,vx_7095_egc_SD,'lineprops','r')
shadedErrorBar([1:12],vx_7095_sea,vx_7095_SD,'lineprops','k')

plot([1:12],vx_7095_sea,'-k','LineWidth',2)
% plot([1:12],vx_7095_egccc_sea,'Color',[0.5 0.9 0.8],'LineWidth',2)
plot([1:12],vx_7095_egcc_sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],vx_7095_egc_sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([-2 1])

% adding in the Havik and Nilsson values
scatter(8,-0.54,50,'blue','filled')
errorbar(8,-0.54,0.28,'-b','LineWidth',1)

scatter(5,-0.7,50,'blue','v','LineWidth',1)


legend('','','','Total','EGCC','EGC','Håvik et al. (2017)','','Nilsson et al.(2008)','location','southeast')
ylabel('Volume Transport [Sv]')
title('Seasonal Cycle of Fresh Water Volume Transport at 70.95 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/7095.jpg']);

%% Looking at the Seasonality of the FWP of the EGCC

%% @ 81.45
egcc_fwp_8145sea = [nanmean(egcc_fwp_8145(jan_i)), nanmean(egcc_fwp_8145(feb_i)),...
    nanmean(egcc_fwp_8145(mar_i)), nanmean(egcc_fwp_8145(apr_i)), nanmean(egcc_fwp_8145(may_i)),...
    nanmean(egcc_fwp_8145(jun_i)), nanmean(egcc_fwp_8145(jul_i)), nanmean(egcc_fwp_8145(aug_i)),...
    nanmean(egcc_fwp_8145(sep_i)), nanmean(egcc_fwp_8145(oct_i)), nanmean(egcc_fwp_8145(nov_i)),...
    nanmean(egcc_fwp_8145(dec_i))];

egc_fwp_8145sea = [nanmean(egc_fwp_8145(jan_i)), nanmean(egc_fwp_8145(feb_i)),...
    nanmean(egc_fwp_8145(mar_i)), nanmean(egc_fwp_8145(apr_i)), nanmean(egc_fwp_8145(may_i)),...
    nanmean(egc_fwp_8145(jun_i)), nanmean(egc_fwp_8145(jul_i)), nanmean(egc_fwp_8145(aug_i)),...
    nanmean(egc_fwp_8145(sep_i)), nanmean(egc_fwp_8145(oct_i)), nanmean(egc_fwp_8145(nov_i)),...
    nanmean(egc_fwp_8145(dec_i))];

figure
hold on
plot([1:12],egcc_fwp_8145sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],egc_fwp_8145sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
% ylim([-1.5 0])

legend('EGCC','EGC','location','northeast')
ylabel('Volume Transport Percentage')
title('Seasonal Cycle of Fresh Water Volume Transport Percentage at 81.45 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/fwvx_per_8145.jpg']);

%% @ 80.05
egcc_fwp_8005sea = [nanmean(egcc_fwp_8005(jan_i)), nanmean(egcc_fwp_8005(feb_i)),...
    nanmean(egcc_fwp_8005(mar_i)), nanmean(egcc_fwp_8005(apr_i)), nanmean(egcc_fwp_8005(may_i)),...
    nanmean(egcc_fwp_8005(jun_i)), nanmean(egcc_fwp_8005(jul_i)), nanmean(egcc_fwp_8005(aug_i)),...
    nanmean(egcc_fwp_8005(sep_i)), nanmean(egcc_fwp_8005(oct_i)), nanmean(egcc_fwp_8005(nov_i)),...
    nanmean(egcc_fwp_8005(dec_i))];

egc_fwp_8005sea = [nanmean(egc_fwp_8005(jan_i)), nanmean(egc_fwp_8005(feb_i)),...
    nanmean(egc_fwp_8005(mar_i)), nanmean(egc_fwp_8005(apr_i)), nanmean(egc_fwp_8005(may_i)),...
    nanmean(egc_fwp_8005(jun_i)), nanmean(egc_fwp_8005(jul_i)), nanmean(egc_fwp_8005(aug_i)),...
    nanmean(egc_fwp_8005(sep_i)), nanmean(egc_fwp_8005(oct_i)), nanmean(egc_fwp_8005(nov_i)),...
    nanmean(egc_fwp_8005(dec_i))];

figure
hold on
plot([1:12],egcc_fwp_8005sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],egc_fwp_8005sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([0 1])

legend('EGCC','EGC','location','northeast')
ylabel('Volume Transport Percentage')
title('Seasonal Cycle of Fresh Water Volume Transport Percentage at 80.05 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/fwvx_per_8005.jpg']);

%% @ 78.85
egcc_fwp_7885sea = [nanmean(egcc_fwp_7885(jan_i)), nanmean(egcc_fwp_7885(feb_i)),...
    nanmean(egcc_fwp_7885(mar_i)), nanmean(egcc_fwp_7885(apr_i)), nanmean(egcc_fwp_7885(may_i)),...
    nanmean(egcc_fwp_7885(jun_i)), nanmean(egcc_fwp_7885(jul_i)), nanmean(egcc_fwp_7885(aug_i)),...
    nanmean(egcc_fwp_7885(sep_i)), nanmean(egcc_fwp_7885(oct_i)), nanmean(egcc_fwp_7885(nov_i)),...
    nanmean(egcc_fwp_7885(dec_i))];

egc_fwp_7885sea = [nanmean(egc_fwp_7885(jan_i)), nanmean(egc_fwp_7885(feb_i)),...
    nanmean(egc_fwp_7885(mar_i)), nanmean(egc_fwp_7885(apr_i)), nanmean(egc_fwp_7885(may_i)),...
    nanmean(egc_fwp_7885(jun_i)), nanmean(egc_fwp_7885(jul_i)), nanmean(egc_fwp_7885(aug_i)),...
    nanmean(egc_fwp_7885(sep_i)), nanmean(egc_fwp_7885(oct_i)), nanmean(egc_fwp_7885(nov_i)),...
    nanmean(egc_fwp_7885(dec_i))];

figure
hold on
plot([1:12],egcc_fwp_7885sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],egc_fwp_7885sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([0 1])

legend('EGCC','EGC','location','northeast')
ylabel('Volume Transport Percentage')
title('Seasonal Cycle of Fresh Water Volume Transport Percentage at 78.85 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/fwvx_per_7885.jpg']);

%% @ 76.95
egcc_fwp_7695sea = [nanmean(egcc_fwp_7695(jan_i)), nanmean(egcc_fwp_7695(feb_i)),...
    nanmean(egcc_fwp_7695(mar_i)), nanmean(egcc_fwp_7695(apr_i)), nanmean(egcc_fwp_7695(may_i)),...
    nanmean(egcc_fwp_7695(jun_i)), nanmean(egcc_fwp_7695(jul_i)), nanmean(egcc_fwp_7695(aug_i)),...
    nanmean(egcc_fwp_7695(sep_i)), nanmean(egcc_fwp_7695(oct_i)), nanmean(egcc_fwp_7695(nov_i)),...
    nanmean(egcc_fwp_7695(dec_i))];

egc_fwp_7695sea = [nanmean(egc_fwp_7695(jan_i)), nanmean(egc_fwp_7695(feb_i)),...
    nanmean(egc_fwp_7695(mar_i)), nanmean(egc_fwp_7695(apr_i)), nanmean(egc_fwp_7695(may_i)),...
    nanmean(egc_fwp_7695(jun_i)), nanmean(egc_fwp_7695(jul_i)), nanmean(egc_fwp_7695(aug_i)),...
    nanmean(egc_fwp_7695(sep_i)), nanmean(egc_fwp_7695(oct_i)), nanmean(egc_fwp_7695(nov_i)),...
    nanmean(egc_fwp_7695(dec_i))];

figure
hold on
plot([1:12],egcc_fwp_7695sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],egc_fwp_7695sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([0 1])

legend('EGCC','EGC','location','northeast')
ylabel('Volume Transport Percentage')
title('Seasonal Cycle of Fresh Water Volume Transport Percentage at 76.95 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/fwvx_per_7695.jpg']);

%% @ 73.95
egcc_fwp_7395sea = [nanmean(egcc_fwp_7395(jan_i)), nanmean(egcc_fwp_7395(feb_i)),...
    nanmean(egcc_fwp_7395(mar_i)), nanmean(egcc_fwp_7395(apr_i)), nanmean(egcc_fwp_7395(may_i)),...
    nanmean(egcc_fwp_7395(jun_i)), nanmean(egcc_fwp_7395(jul_i)), nanmean(egcc_fwp_7395(aug_i)),...
    nanmean(egcc_fwp_7395(sep_i)), nanmean(egcc_fwp_7395(oct_i)), nanmean(egcc_fwp_7395(nov_i)),...
    nanmean(egcc_fwp_7395(dec_i))];

egc_fwp_7395sea = [nanmean(egc_fwp_7395(jan_i)), nanmean(egc_fwp_7395(feb_i)),...
    nanmean(egc_fwp_7395(mar_i)), nanmean(egc_fwp_7395(apr_i)), nanmean(egc_fwp_7395(may_i)),...
    nanmean(egc_fwp_7395(jun_i)), nanmean(egc_fwp_7395(jul_i)), nanmean(egc_fwp_7395(aug_i)),...
    nanmean(egc_fwp_7395(sep_i)), nanmean(egc_fwp_7395(oct_i)), nanmean(egc_fwp_7395(nov_i)),...
    nanmean(egc_fwp_7395(dec_i))];

figure
hold on
plot([1:12],egcc_fwp_7395sea,'Color',[0.5 0.5 1],'LineWidth',2)
plot([1:12],egc_fwp_7395sea,'Color',[0.7 0 0.4],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([0 1])

legend('EGCC','EGC','location','northeast')
ylabel('Volume Transport Percentage')
title('Seasonal Cycle of Fresh Water Volume Transport Percentage at 73.95 ^oN (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/fwvx_per_7395.jpg']);

%% Putting all of the volume transport percentages on the same figure

figure
hold on

% plot([1:12],egcc_fwp_8145sea,'Color',[0 0 0.4],'LineWidth',2)
plot([1:12],egcc_fwp_8005sea,'Color',[0 0 0],'LineWidth',2)
plot([1:12],egcc_fwp_7885sea,'Color',[0 0 1],'LineWidth',2)
plot([1:12],egcc_fwp_7695sea,'Color',[0 0.6 1],'LineWidth',2)
plot([1:12],egcc_fwp_7395sea,'Color',[0 0.9 1],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([0 1])

legend('80.05','78.85','76.95','73.95','location','northeast')
ylabel('Volume Transport Percentage')
title('Seasonal Cycle of EGCC Fresh Water Volume Transport Percentage (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/FwVx_egcc_perc.jpg']);

%% Putting all of the Fresh Water volume transport on the same figure

figure
hold on

plot([1:12],vx_8145_egcc_sea,'Color',[0 0 0],'LineWidth',2)
plot([1:12],vx_8005_egcc_sea,'Color',[0 0 1],'LineWidth',2)
plot([1:12],vx_7885_egcc_sea,'Color',[0 0.6 1],'LineWidth',2)
plot([1:12],vx_7695_egcc_sea,'Color',[0 0.9 1],'LineWidth',2)
plot([1:12],vx_7395_egcc_sea,'Color',[0 1 0.7],'LineWidth',2)

xticks(1:12);
xlim([1 12])
xticklabels(months)
ylim([-1.5 0])

legend('8145','80.05','78.85','76.95','73.95','location','southeast')
ylabel('Volume Transport [Sv]')
title('Seasonal Cycle of EGCC Fresh Water Volume Transport (RARE)')
subtitle('2003-2019')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_VxFW/Seasonal Cycle/FwVx_egcc.jpg']);
