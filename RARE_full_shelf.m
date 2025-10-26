% RARE full shelf
% Created by: Sara Vianco
% Masters Thesis
% 15 June 2023

%% Preliminaries
close all
clear all

cd('/Volumes/TERABYTE/RARE/matching_theo/')
addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/')

%% initial info/variable establishment
info = ncinfo("rare1.15.2_mn_ocean_reg_2003.nc");
lon = ncread('rare1.15.2_mn_ocean_reg_2003.nc','xt_ocean');
lat = ncread('rare1.15.2_mn_ocean_reg_2003.nc','yt_ocean');

%% specify lat and lon limits
full_lat_i = find(lat == 78.85);
full_lon_i = find(lon>-18 & lon<0);
start_full = [full_lon_i(1) full_lat_i 1 1];
count_full = [length(full_lon_i) 1 38 12];


%% reading in multiple files

ncvars =  {'time', 'xt_ocean', 'yt_ocean', 'sw_ocean', 'v','temp','salt'};
projectdir = '/Volumes/TERABYTE/RARE/matching_theo/';
dinfo = dir( fullfile(projectdir, '*.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );

time = cell(num_files, 1);
xt_ocean = cell(num_files, 1); 
yt_ocean = cell(num_files, 1);
sw_ocean = cell(num_files, 1);
v = cell(num_files, 1);
ml_rho = cell(num_files,1);
v_ml = cell(num_files,1);
temp = cell(num_files,1);
salt = cell(num_files,1);


for K = 1:num_files
  this_file = filenames{K};
  time{K} = ncread(this_file, ncvars{1});
  xt_ocean{K} = ncread(this_file, ncvars{2});
  yt_ocean{K} = ncread(this_file, ncvars{3});
  sw_ocean{K} = ncread(this_file, ncvars{4});
  v{K} = ncread(this_file, ncvars{5},start_full,count_full);
  temp_pot{K} = ncread(this_file, ncvars{6},start_full,count_full);
  sal_pr{K} = ncread(this_file, ncvars{7},start_full,count_full);
end

lon_cell = (xt_ocean{1});
lon_full = lon_cell(full_lon_i);
lat_full = lat(full_lat_i);
lat_cell = (yt_ocean{1});

z_cell = sw_ocean{1};


%% adjusting the time

% making the date number into one that is readable by MATLAB, forming a
% date string afterwards

base = datenum(1979,12,31,0,0,0);
formatIn = 'mm/dd/yy';

for K = 1:num_files
   
date_num{K} = datenum(time{K} + base);
date_string{K} = datestr(date_num{K},'mm/dd/yy');
dv{K} = datevec(date_string{K},formatIn);

end

%% Separating into seperate months

% January
for K = 1:num_files
jan_v_cells{K} = squeeze(mean(v{K}(:,:,:,1),4,'omitnan'));
jan_v_mean = cell2mat(jan_v_cells);
end

[r,c] = size(jan_v_mean);
nlay  = 17;
jan_v_mat  = permute(reshape(jan_v_mean,[r,c/nlay,nlay]),[1,2,3]);
jan_v_mean = mean(jan_v_mat,3,'omitnan');

% February
for K = 1:num_files
feb_v_cells{K} = squeeze(mean(v{K}(:,:,:,2),4,'omitnan'));
feb_v_mean = cell2mat(feb_v_cells);
end

[r,c] = size(feb_v_mean);
nlay  = 17;
feb_v_mat  = permute(reshape(feb_v_mean,[r,c/nlay,nlay]),[1,2,3]);
feb_v_mean = mean(feb_v_mat,3,'omitnan');

% March
for K = 1:num_files
mar_v_cells{K} = squeeze(mean(v{K}(:,:,:,3),4,'omitnan'));
mar_v_mean = cell2mat(mar_v_cells);
end

[r,c] = size(mar_v_mean);
nlay  = 17;
mar_v_mat   = permute(reshape(mar_v_mean,[r,c/nlay,nlay]),[1,2,3]);
mar_v_mean = mean(mar_v_mat,3,'omitnan');

% April
for K = 1:num_files
apr_v_cells{K} = squeeze(mean(v{K}(:,:,:,4),4,'omitnan'));
apr_v_mean = cell2mat(apr_v_cells);
end

[r,c] = size(apr_v_mean);
nlay  = 17;
apr_v_mat   = permute(reshape(apr_v_mean,[r,c/nlay,nlay]),[1,2,3]);
apr_v_mean = mean(apr_v_mat,3,'omitnan');

% May
for K = 1:num_files
may_v_cells{K} = squeeze(mean(v{K}(:,:,:,5),4,'omitnan'));
may_v_mean = cell2mat(may_v_cells);
end

[r,c] = size(may_v_mean);
nlay  = 17;
may_v_mat   = permute(reshape(may_v_mean,[r,c/nlay,nlay]),[1,2,3]);
may_v_mean = mean(may_v_mat,3,'omitnan');

% June
for K = 1:num_files
jun_v_cells{K} = squeeze(mean(v{K}(:,:,:,6),4,'omitnan'));
jun_v_mean = cell2mat(jun_v_cells);
end

[r,c] = size(jun_v_mean);
nlay  = 17;
jun_v_mat   = permute(reshape(jun_v_mean,[r,c/nlay,nlay]),[1,2,3]);
jun_v_mean = mean(jun_v_mat,3,'omitnan');

% July
for K = 1:num_files
jul_v_cells{K} = squeeze(mean(v{K}(:,:,:,7),4,'omitnan'));
jul_v_mean = cell2mat(jul_v_cells);
end

[r,c] = size(jul_v_mean);
nlay  = 17;
jul_v_mat   = permute(reshape(jul_v_mean,[r,c/nlay,nlay]),[1,2,3]);
jul_v_mean = mean(jul_v_mat,3,'omitnan');

% August
for K = 1:num_files
aug_v_cells{K} = squeeze(mean(v{K}(:,:,:,8),4,'omitnan'));
aug_v_mean = cell2mat(aug_v_cells);
end

[r,c] = size(aug_v_mean);
nlay  = 17;
aug_v_mat   = permute(reshape(aug_v_mean,[r,c/nlay,nlay]),[1,2,3]);
aug_v_mean = mean(aug_v_mat,3,'omitnan');

% September
for K = 1:num_files
sep_v_cells{K} = squeeze(mean(v{K}(:,:,:,9),4,'omitnan'));
sep_v_mean = cell2mat(sep_v_cells);
end

[r,c] = size(sep_v_mean);
nlay  = 17;
sep_v_mat   = permute(reshape(sep_v_mean,[r,c/nlay,nlay]),[1,2,3]);
sep_v_mean = mean(sep_v_mat,3,'omitnan');

% October
for K = 1:num_files
oct_v_cells{K} = squeeze(mean(v{K}(:,:,:,10),4,'omitnan'));
oct_v_mean = cell2mat(oct_v_cells);
end

[r,c] = size(oct_v_mean);
nlay  = 17;
oct_v_mat   = permute(reshape(oct_v_mean,[r,c/nlay,nlay]),[1,2,3]);
oct_v_mean = mean(oct_v_mat,3,'omitnan');

% November
for K = 1:num_files
nov_v_cells{K} = squeeze(mean(v{K}(:,:,:,11),4,'omitnan'));
nov_v_mean = cell2mat(nov_v_cells);
end

[r,c] = size(nov_v_mean);
nlay  = 17;
nov_v_mat   = permute(reshape(nov_v_mean,[r,c/nlay,nlay]),[1,2,3]);
nov_v_mean = mean(nov_v_mat,3,'omitnan');

% December
for K = 1:num_files
dec_v_cells{K} = squeeze(mean(v{K}(:,:,:,12),4,'omitnan'));
dec_v_mean = cell2mat(dec_v_cells);
end

[r,c] = size(dec_v_mean);
nlay  = 17;
dec_v_mat   = permute(reshape(dec_v_mean,[r,c/nlay,nlay]),[1,2,3]);
dec_v_mean = mean(dec_v_mat,3,'omitnan');

%% Initial Imaging: Monthly Averages

[Z,LON] = meshgrid(-z_cell,lon_full);

% % January
% figure
% hold on
% contourf(LON,Z,jan_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,jan_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% 
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(jan_v_mean))=0;
% imagesc(jan_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean January EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/01_RARE.jpg']);
% 
% % February
% figure
% hold on
% contourf(LON,Z,feb_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,feb_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(feb_v_mean))=0;
% imagesc(feb_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean February EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/02_RARE.jpg']);
% 
% % March
% figure
% hold on
% contourf(LON,Z,mar_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(mar_v_mean))=0;
% imagesc(mar_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,mar_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean March EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/03_RARE.jpg']);
% 
% % April
% figure
% hold on
% contourf(LON,Z,apr_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(apr_v_mean))=0;
% imagesc(apr_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,apr_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean April EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/04_RARE.jpg']);
% 
% % May
% figure
% hold on
% contourf(LON,Z,may_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(may_v_mean))=0;
% imagesc(may_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,may_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean May EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/05_RARE.jpg']);
% 
% % June
% figure
% hold on
% contourf(LON,Z,jun_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(jun_v_mean))=0;
% imagesc(jun_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% 
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,jun_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean June EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/06_RARE.jpg']);
% 
% % July
% figure
% hold on
% contourf(LON,Z,jul_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(jul_v_mean))=0;
% imagesc(jul_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% 
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,jul_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean July EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/07_RARE.jpg']);
% 
% % August
% figure
% hold on
% contourf(LON,Z,aug_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(aug_v_mean))=0;
% imagesc(aug_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% 
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,aug_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean August EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/08_RARE.jpg']);
% 
% 
% % September
% figure
% hold on
% contourf(LON,Z,sep_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(sep_v_mean))=0;
% imagesc(sep_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% 
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,sep_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean September EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/09_RARE.jpg']);
% 
% % October
% figure
% hold on
% contourf(LON,Z,oct_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(oct_v_mean))=0;
% imagesc(oct_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,oct_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean October EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/10_RARE.jpg']);
% 
% % November
% figure
% hold on
% contourf(LON,Z,nov_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(nov_v_mean))=0;
% imagesc(nov_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% xlim([-18 0]);
% 
% ylim([-1000 0]);
% %ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,nov_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean November EGS Flow (RARE)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'RARE/RARE_full_shelf/11_RARE.jpg']);
% 
% % December
% figure
% hold on
% 
% imAlpha=ones(size(dec_v_mean));
% imAlpha(isnan(dec_v_mean))=0;
% imagesc(dec_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% contourf(LON,Z,dec_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% 
% xlim([-18 0]);
% ylim([-1000 0]);
% %ylim([-2650 0]);
% 
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,dec_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean December EGS Flow (RARE)','FontSize',16);
% 
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%      'RARE/RARE_full_shelf/12_RARE.jpg']);

%% Establishing Density Contours: Full Shelf %%

% adding gibbs toolbox to path
addpath('/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16');
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/html/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/library/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/thermodynamics_from_t/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/pdf/");

% calculating pressure variable from depth variable (z_cell)
p = gsw_p_from_z(-z_cell,lat_full);
p = p';
p_mat = repmat(p,length(lon_full),1);

%% Determining practical salinity (SP) variable
% January
for K = 1:num_files
jan_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,1),4,'omitnan'));
jan_s_mean = cell2mat(jan_s_cells);
end

[r,c] = size(jan_s_mean);
nlay  = 17;
jan_v_mat  = permute(reshape(jan_s_mean,[r,c/nlay,nlay]),[1,2,3]);
jan_s_mean = mean(jan_v_mat,3,'omitnan');

% February
for K = 1:num_files
feb_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,2),4,'omitnan'));
feb_s_mean = cell2mat(feb_s_cells);
end

[r,c] = size(feb_s_mean);
nlay  = 17;
feb_v_mat  = permute(reshape(feb_s_mean,[r,c/nlay,nlay]),[1,2,3]);
feb_s_mean = mean(feb_v_mat,3,'omitnan');

% March
for K = 1:num_files
mar_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,3),4,'omitnan'));
mar_s_mean = cell2mat(mar_s_cells);
end

[r,c] = size(mar_s_mean);
nlay  = 17;
mar_v_mat   = permute(reshape(mar_s_mean,[r,c/nlay,nlay]),[1,2,3]);
mar_s_mean = mean(mar_v_mat,3,'omitnan');

% April
for K = 1:num_files
apr_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,4),4,'omitnan'));
apr_s_mean = cell2mat(apr_s_cells);
end

[r,c] = size(apr_s_mean);
nlay  = 17;
apr_v_mat   = permute(reshape(apr_s_mean,[r,c/nlay,nlay]),[1,2,3]);
apr_s_mean = mean(apr_v_mat,3,'omitnan');

% May
for K = 1:num_files
may_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,5),4,'omitnan'));
may_s_mean = cell2mat(may_s_cells);
end

[r,c] = size(may_s_mean);
nlay  = 17;
may_v_mat   = permute(reshape(may_s_mean,[r,c/nlay,nlay]),[1,2,3]);
may_s_mean = mean(may_v_mat,3,'omitnan');

% June
for K = 1:num_files
jun_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,6),4,'omitnan'));
jun_s_mean = cell2mat(jun_s_cells);
end

[r,c] = size(jun_s_mean);
nlay  = 17;
jun_v_mat   = permute(reshape(jun_s_mean,[r,c/nlay,nlay]),[1,2,3]);
jun_s_mean = mean(jun_v_mat,3,'omitnan');

% July
for K = 1:num_files
jul_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,7),4,'omitnan'));
jul_s_mean = cell2mat(jul_s_cells);
end

[r,c] = size(jul_s_mean);
nlay  = 17;
jul_v_mat   = permute(reshape(jul_s_mean,[r,c/nlay,nlay]),[1,2,3]);
jul_s_mean = mean(jul_v_mat,3,'omitnan');

% August
for K = 1:num_files
aug_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,8),4,'omitnan'));
aug_s_mean = cell2mat(aug_s_cells);
end

[r,c] = size(aug_s_mean);
nlay  = 17;
aug_v_mat   = permute(reshape(aug_s_mean,[r,c/nlay,nlay]),[1,2,3]);
aug_s_mean = mean(aug_v_mat,3,'omitnan');

% September
for K = 1:num_files
sep_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,9),4,'omitnan'));
sep_s_mean = cell2mat(sep_s_cells);
end

[r,c] = size(sep_s_mean);
nlay  = 17;
sep_v_mat   = permute(reshape(sep_s_mean,[r,c/nlay,nlay]),[1,2,3]);
sep_s_mean = mean(sep_v_mat,3,'omitnan');

% October
for K = 1:num_files
oct_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,10),4,'omitnan'));
oct_s_mean = cell2mat(oct_s_cells);
end

[r,c] = size(oct_s_mean);
nlay  = 17;
oct_v_mat   = permute(reshape(oct_s_mean,[r,c/nlay,nlay]),[1,2,3]);
oct_s_mean = mean(oct_v_mat,3,'omitnan');

% November
for K = 1:num_files
nov_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,11),4,'omitnan'));
nov_s_mean = cell2mat(nov_s_cells);
end

[r,c] = size(nov_s_mean);
nlay  = 17;
nov_v_mat   = permute(reshape(nov_s_mean,[r,c/nlay,nlay]),[1,2,3]);
nov_s_mean = mean(nov_v_mat,3,'omitnan');

% December
for K = 1:num_files
dec_s_cells{K} = squeeze(mean(sal_pr{K}(:,:,:,12),4,'omitnan'));
dec_s_mean = cell2mat(dec_s_cells);
end

[r,c] = size(dec_s_mean);
nlay  = 17;
dec_v_mat   = permute(reshape(dec_s_mean,[r,c/nlay,nlay]),[1,2,3]);
dec_s_mean = mean(dec_v_mat,3,'omitnan');

%% Establishing Potential Temperature variables (by month)
% January
for K = 1:num_files
jan_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,1),4,'omitnan'));
jan_t_mean = cell2mat(jan_t_cells);
end

[r,c] = size(jan_t_mean);
nlay  = 17;
jan_v_mat  = permute(reshape(jan_t_mean,[r,c/nlay,nlay]),[1,2,3]);
jan_t_mean = mean(jan_v_mat,3,'omitnan');

% February
for K = 1:num_files
feb_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,2),4,'omitnan'));
feb_t_mean = cell2mat(feb_t_cells);
end

[r,c] = size(feb_t_mean);
nlay  = 17;
feb_v_mat  = permute(reshape(feb_t_mean,[r,c/nlay,nlay]),[1,2,3]);
feb_t_mean = mean(feb_v_mat,3,'omitnan');

% March
for K = 1:num_files
mar_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,3),4,'omitnan'));
mar_t_mean = cell2mat(mar_t_cells);
end

[r,c] = size(mar_t_mean);
nlay  = 17;
mar_v_mat   = permute(reshape(mar_t_mean,[r,c/nlay,nlay]),[1,2,3]);
mar_t_mean = mean(mar_v_mat,3,'omitnan');

% April
for K = 1:num_files
apr_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,4),4,'omitnan'));
apr_t_mean = cell2mat(apr_t_cells);
end

[r,c] = size(apr_t_mean);
nlay  = 17;
apr_v_mat   = permute(reshape(apr_t_mean,[r,c/nlay,nlay]),[1,2,3]);
apr_t_mean = mean(apr_v_mat,3,'omitnan');

% May
for K = 1:num_files
may_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,5),4,'omitnan'));
may_t_mean = cell2mat(may_t_cells);
end

[r,c] = size(may_t_mean);
nlay  = 17;
may_v_mat   = permute(reshape(may_t_mean,[r,c/nlay,nlay]),[1,2,3]);
may_t_mean = mean(may_v_mat,3,'omitnan');

% June
for K = 1:num_files
jun_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,6),4,'omitnan'));
jun_t_mean = cell2mat(jun_t_cells);
end

[r,c] = size(jun_t_mean);
nlay  = 17;
jun_v_mat   = permute(reshape(jun_t_mean,[r,c/nlay,nlay]),[1,2,3]);
jun_t_mean = mean(jun_v_mat,3,'omitnan');

% July
for K = 1:num_files
jul_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,7),4,'omitnan'));
jul_t_mean = cell2mat(jul_t_cells);
end

[r,c] = size(jul_t_mean);
nlay  = 17;
jul_v_mat   = permute(reshape(jul_t_mean,[r,c/nlay,nlay]),[1,2,3]);
jul_t_mean = mean(jul_v_mat,3,'omitnan');

% August
for K = 1:num_files
aug_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,8),4,'omitnan'));
aug_t_mean = cell2mat(aug_t_cells);
end

[r,c] = size(aug_t_mean);
nlay  = 17;
aug_v_mat   = permute(reshape(aug_t_mean,[r,c/nlay,nlay]),[1,2,3]);
aug_t_mean = mean(aug_v_mat,3,'omitnan');

% September
for K = 1:num_files
sep_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,9),4,'omitnan'));
sep_t_mean = cell2mat(sep_t_cells);
end

[r,c] = size(sep_t_mean);
nlay  = 17;
sep_v_mat   = permute(reshape(sep_t_mean,[r,c/nlay,nlay]),[1,2,3]);
sep_t_mean = mean(sep_v_mat,3,'omitnan');

% October
for K = 1:num_files
oct_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,10),4,'omitnan'));
oct_t_mean = cell2mat(oct_t_cells);
end

[r,c] = size(oct_t_mean);
nlay  = 17;
oct_v_mat   = permute(reshape(oct_t_mean,[r,c/nlay,nlay]),[1,2,3]);
oct_t_mean = mean(oct_v_mat,3,'omitnan');

% November
for K = 1:num_files
nov_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,11),4,'omitnan'));
nov_t_mean = cell2mat(nov_t_cells);
end

[r,c] = size(nov_t_mean);
nlay  = 17;
nov_v_mat   = permute(reshape(nov_t_mean,[r,c/nlay,nlay]),[1,2,3]);
nov_t_mean = mean(nov_v_mat,3,'omitnan');

% December
for K = 1:num_files
dec_t_cells{K} = squeeze(mean(temp_pot{K}(:,:,:,12),4,'omitnan'));
dec_t_mean = cell2mat(dec_t_cells);
end

[r,c] = size(dec_t_mean);
nlay  = 17;
dec_v_mat   = permute(reshape(dec_t_mean,[r,c/nlay,nlay]),[1,2,3]);
dec_t_mean = mean(dec_v_mat,3,'omitnan');

%% combining 2D monthly means into arrays
s_mean = cat(3,jan_s_mean, feb_s_mean, mar_s_mean, apr_s_mean, may_s_mean, jun_s_mean,... 
    jul_s_mean, aug_s_mean, sep_s_mean, oct_s_mean, nov_s_mean, dec_s_mean);

t_mean = cat(3,jan_t_mean, feb_t_mean, mar_t_mean, apr_t_mean, may_t_mean, jun_t_mean,... 
    jul_t_mean, aug_t_mean, sep_t_mean, oct_t_mean, nov_t_mean, dec_t_mean);

%% Individual variables for density calculation

SA = NaN(size(s_mean));
for i = 1:12
[SA(:,:,i), in_ocean] = gsw_SA_from_SP(s_mean(:,:,i),p_mat,lon_full,lat_full);
end

CT = NaN(size(t_mean));
for i = 1:12
CT(:,:,i) = gsw_CT_from_pt(SA(:,:,i),t_mean(:,:,i));
end


%% Density Calculation

rho_CT_exact = NaN(size(t_mean));
for i = 1:12
rho_CT_exact(:,:,i) = gsw_rho_CT_exact(SA(:,:,i),CT(:,:,i),p_mat);
end

rho_anom = rho_CT_exact-1000;

%% Adding Mooring Array from FSAOO for Comparison
addpath('/Users/saravianco/Documents/Research/MATLAB/theo/')
load('mooring_xy.mat')

%% Imaging Density Contours

% January
figure
hold on
contourf(LON,Z,jan_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,1),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(jan_v_mean));
imAlpha(isnan(jan_v_mean))=0;
imagesc(jan_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean January EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/01_rho.jpg']);

% February
figure
hold on
contourf(LON,Z,feb_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,2),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(feb_v_mean));
imAlpha(isnan(feb_v_mean))=0;
imagesc(feb_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean February EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/02_rho.jpg']);

% March
figure
hold on
contourf(LON,Z,mar_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,3),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(mar_v_mean));
imAlpha(isnan(mar_v_mean))=0;
imagesc(mar_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean March EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/03_rho.jpg']);

% April
figure
hold on
contourf(LON,Z,apr_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,4),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(apr_v_mean));
imAlpha(isnan(apr_v_mean))=0;
imagesc(apr_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean April EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/04_rho.jpg']);

% May
figure
hold on
contourf(LON,Z,may_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,5),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(may_v_mean));
imAlpha(isnan(may_v_mean))=0;
imagesc(may_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean May EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/05_rho.jpg']);

% June
figure
hold on
contourf(LON,Z,jun_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,6),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(jun_v_mean));
imAlpha(isnan(jun_v_mean))=0;
imagesc(jun_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean June EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/06_rho.jpg']);

% July
figure
hold on
contourf(LON,Z,jul_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,7),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(jul_v_mean));
imAlpha(isnan(jul_v_mean))=0;
imagesc(jul_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean July EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/07_rho.jpg']);

% August
figure
hold on
contourf(LON,Z,aug_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,8),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(aug_v_mean));
imAlpha(isnan(aug_v_mean))=0;
imagesc(aug_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean August EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/08_rho.jpg']);

% September
figure
hold on
contourf(LON,Z,sep_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,9),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(sep_v_mean));
imAlpha(isnan(sep_v_mean))=0;
imagesc(sep_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean September EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/09_rho.jpg']);

% October
figure
hold on
contourf(LON,Z,oct_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,10),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(oct_v_mean));
imAlpha(isnan(oct_v_mean))=0;
imagesc(oct_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean October EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/10_rho.jpg']);

% November
figure
hold on
contourf(LON,Z,nov_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,11),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(nov_v_mean));
imAlpha(isnan(nov_v_mean))=0;
imagesc(nov_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean November EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/11_rho.jpg']);

% December
figure
hold on
contourf(LON,Z,dec_v_mean,8,'linestyle','none');
shading interp
x = colorbar;
[C,d] = contour(LON,Z,rho_anom(:,:,12),[25:0.5:32.5],'-k');
clabel(C,d);

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

imAlpha=ones(size(dec_v_mean));
imAlpha(isnan(dec_v_mean))=0;
imagesc(dec_v_mean,'AlphaData',imAlpha);
set(gca,'color',[0.5 0.5 0.5]);

xlim([-8 -2]);
ylim([-1000 0]);
clim([-0.25 0.11]);
cmocean('balance','pivot',0)

xlabel('Longitude [^oE]','FontSize',14)
ylabel('Depth [m]','FontSize',14)
ylabel(x,'V [m/s]','FontSize',14)
title('Mean December EGS Velocity and Density (RARE)','FontSize',16);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'RARE/RARE_monthly_density/12_rho.jpg']);

%% Volumetric Transport: EGCC and EGCCC

% Establishing fixed box for EGCC transport
% lon_1 = -14.05;
% lon_2 = -6.05;
% 
% lon1_idx = find(lon_full<=-14 & lon_full>=-14.07);
% lon2_idx = find(lon_full<=-6 & lon_full>=-6.07);
% z_idx = find(z_cell>=250,1,'first');
% 
% gridx_diff = lon2_idx-lon1_idx;
% 
% lon_km = haversine([lat_full,lon_1],[lat_full,lon_2]);
% lon_m = lon_km*1000;
% dx = lon_m/gridx_diff; % horizontal distance of a gridcell
% 
% dz = NaN(z_idx,1);
% for i = 1:z_idx
%     dz(i) = z_cell(i+1)-z_cell(i);
% end
% 
% % January
% jan_egcc = jan_v_mean(lon1_idx:lon2_idx,1:z_idx);
% jan_egcc(isnan(jan_egcc)) = 0;
% tx_jan = (sum(dx*jan_egcc*dz))/10^6;
% 
% % February
% feb_egcc = feb_v_mean(lon1_idx:lon2_idx,1:z_idx);
% feb_egcc(isnan(feb_egcc)) = 0;
% tx_feb = (sum(dx*feb_egcc*dz))/10^6;
% 
% % March
% mar_egcc = mar_v_mean(lon1_idx:lon2_idx,1:z_idx);
% mar_egcc(isnan(mar_egcc)) = 0;
% tx_mar = (sum(dx*mar_egcc*dz))/10^6;
% 
% % April
% apr_egcc = apr_v_mean(lon1_idx:lon2_idx,1:z_idx);
% apr_egcc(isnan(apr_egcc)) = 0;
% tx_apr = (sum(dx*apr_egcc*dz))/10^6;
% 
% % May
% may_egcc = may_v_mean(lon1_idx:lon2_idx,1:z_idx);
% may_egcc(isnan(may_egcc)) = 0;
% tx_may = (sum(dx*may_egcc*dz))/10^6;
% 
% % June
% jun_egcc = jun_v_mean(lon1_idx:lon2_idx,1:z_idx);
% jun_egcc(isnan(jun_egcc)) = 0;
% tx_jun = (sum(dx*jun_egcc*dz))/10^6;
% 
% % July
% jul_egcc = jul_v_mean(lon1_idx:lon2_idx,1:z_idx);
% jul_egcc(isnan(jul_egcc)) = 0;
% tx_jul = (sum(dx*jul_egcc*dz))/10^6;
% 
% % August
% aug_egcc = aug_v_mean(lon1_idx:lon2_idx,1:z_idx);
% aug_egcc(isnan(aug_egcc)) = 0;
% tx_aug = (sum(dx*aug_egcc*dz))/10^6;
% 
% % September
% sep_egcc = sep_v_mean(lon1_idx:lon2_idx,1:z_idx);
% sep_egcc(isnan(sep_egcc)) = 0;
% tx_sep = (sum(dx*sep_egcc*dz))/10^6;
% 
% % October
% oct_egcc = oct_v_mean(lon1_idx:lon2_idx,1:z_idx);
% oct_egcc(isnan(oct_egcc)) = 0;
% tx_oct = (sum(dx*oct_egcc*dz))/10^6;
% 
% % November
% nov_egcc = nov_v_mean(lon1_idx:lon2_idx,1:z_idx);
% nov_egcc(isnan(nov_egcc)) = 0;
% tx_nov = (sum(dx*nov_egcc*dz))/10^6;
% 
% % December
% dec_egcc = dec_v_mean(lon1_idx:lon2_idx,1:z_idx);
% dec_egcc(isnan(dec_egcc)) = 0;
% tx_dec = (sum(dx*dec_egcc*dz))/10^6;
% 
% % concatenating monthly volumetric transport values into a single variable
% monthly_tx = [tx_jan,tx_feb,tx_mar,tx_apr,tx_may,tx_jun,tx_jul,tx_aug,...
%     tx_sep,tx_oct,tx_nov,tx_dec];
% month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% 
% figure
% plot(monthly_tx,'-k','linewi',2)
% set(gca,'Ydir','reverse')
% %ylim([-5.3 -2.7])
% xlim([1 12])
% xticks(1:12);
% xticklabels(month)
% 
% ylabel('Transport [Sv]')
% title('EGCC Average Monthly Volumetric Transport (RARE)')
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_vol_tx/EGCC_vol_tx.jpg']);
% 
% %%
% % Establishing fixed box for EGCCC transport
% lon_1 = lon_full(1);
% lon_2 = -14.05;
% 
% lon1_idx = 1;
% lon2_idx = find(lon_full<=-14 & lon_full>=-14.07);
% z_idx = find(z_cell>=250,1,'first');
% 
% gridx_diff = lon2_idx-lon1_idx;
% 
% lon_km = haversine([lat_full,lon_1],[lat_full,lon_2]);
% lon_m = lon_km*1000;
% dx = lon_m/gridx_diff; % horizontal distance of a gridcell
% 
% dz = NaN(z_idx,1);
% for i = 1:z_idx
%     dz(i) = z_cell(i+1)-z_cell(i);
% end
% 
% % January
% jan_egccc = jan_v_mean(lon1_idx:lon2_idx,1:z_idx);
% jan_egccc(isnan(jan_egccc)) = 0;
% tx_jan = (sum(dx*jan_egccc*dz))/10^6;
% 
% % February
% feb_egccc = feb_v_mean(lon1_idx:lon2_idx,1:z_idx);
% feb_egccc(isnan(feb_egccc)) = 0;
% tx_feb = (sum(dx*feb_egccc*dz))/10^6;
% 
% % March
% mar_egccc = mar_v_mean(lon1_idx:lon2_idx,1:z_idx);
% mar_egccc(isnan(mar_egccc)) = 0;
% tx_mar = (sum(dx*mar_egccc*dz))/10^6;
% 
% % April
% apr_egccc = apr_v_mean(lon1_idx:lon2_idx,1:z_idx);
% apr_egccc(isnan(apr_egccc)) = 0;
% tx_apr = (sum(dx*apr_egccc*dz))/10^6;
% 
% % May
% may_egccc = may_v_mean(lon1_idx:lon2_idx,1:z_idx);
% may_egccc(isnan(may_egccc)) = 0;
% tx_may = (sum(dx*may_egccc*dz))/10^6;
% 
% % June
% jun_egccc = jun_v_mean(lon1_idx:lon2_idx,1:z_idx);
% jun_egccc(isnan(jun_egccc)) = 0;
% tx_jun = (sum(dx*jun_egccc*dz))/10^6;
% 
% % July
% jul_egccc = jul_v_mean(lon1_idx:lon2_idx,1:z_idx);
% jul_egccc(isnan(jul_egccc)) = 0;
% tx_jul = (sum(dx*jul_egccc*dz))/10^6;
% 
% % August
% aug_egccc = aug_v_mean(lon1_idx:lon2_idx,1:z_idx);
% aug_egccc(isnan(aug_egccc)) = 0;
% tx_aug = (sum(dx*aug_egccc*dz))/10^6;
% 
% % September
% sep_egccc = sep_v_mean(lon1_idx:lon2_idx,1:z_idx);
% sep_egccc(isnan(sep_egccc)) = 0;
% tx_sep = (sum(dx*sep_egccc*dz))/10^6;
% 
% % October
% oct_egccc = oct_v_mean(lon1_idx:lon2_idx,1:z_idx);
% oct_egccc(isnan(oct_egccc)) = 0;
% tx_oct = (sum(dx*oct_egccc*dz))/10^6;
% 
% % November
% nov_egccc = nov_v_mean(lon1_idx:lon2_idx,1:z_idx);
% nov_egccc(isnan(nov_egccc)) = 0;
% tx_nov = (sum(dx*nov_egccc*dz))/10^6;
% 
% % December
% dec_egccc = dec_v_mean(lon1_idx:lon2_idx,1:z_idx);
% dec_egccc(isnan(dec_egccc)) = 0;
% tx_dec = (sum(dx*dec_egccc*dz))/10^6;
% 
% % concatenating monthly volumetric transport values into a single variable
% monthly_tx = [tx_jan,tx_feb,tx_mar,tx_apr,tx_may,tx_jun,tx_jul,tx_aug,...
%     tx_sep,tx_oct,tx_nov,tx_dec];
% month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% 
% figure
% plot(monthly_tx,'-k','linewi',2)
% set(gca,'Ydir','reverse')
% %ylim([-5.3 -2.7])
% xlim([1 12])
% xticks(1:12);
% xticklabels(month)
% 
% ylabel('Transport [Sv]')
% title('EGCCC Average Monthly Volumetric Transport (RARE)')
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_vol_tx/EGCCC_vol_tx.jpg']);
