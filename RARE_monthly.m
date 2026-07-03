% RARE (monthly averages)
% Created by: Sara Vianco
% Masters Thesis
% 24 MAY 2023

%% Preliminaries
close all
clear all

cd('/Volumes/One Touch/Research Datasets/RARE/matching_theo/')
addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/')

%% initial info/variable establishment
info = ncinfo("rare1.15.2_mn_ocean_reg_2003.nc");
info2 = ncinfo("rare1.15.2_mn_ocean_reg_2019.nc");

time = ncread('rare1.15.2_mn_ocean_reg_2003.nc','time');
lon = ncread('rare1.15.2_mn_ocean_reg_2003.nc','xt_ocean');
lat = ncread('rare1.15.2_mn_ocean_reg_2003.nc','yt_ocean');
z_star = ncread('rare1.15.2_mn_ocean_reg_2003.nc','sw_ocean');
[LON,LAT] = meshgrid(lon,lat);

%% specify FSAOO location
fs_lat_i = find(lat == 78.85);
fs_lon_i = find(lon>-8 & lon<-2);
start_fs = [fs_lon_i(1) fs_lat_i 1 1];
count_fs = [length(fs_lon_i) 1 38 12];


%% reading in multiple files

ncvars =  {'time', 'xt_ocean', 'yt_ocean', 'sw_ocean', 'v','temp','salt'};
projectdir = '/Volumes/One Touch/Research Datasets/RARE/matching_theo/';
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
  v{K} = ncread(this_file, ncvars{5},start_fs,count_fs);
end

for K = 1:num_files
  ptemp{K} = ncread(this_file, ncvars{6},start_fs,count_fs);
  Sp{K} = ncread(this_file, ncvars{7},start_fs,count_fs);
end

lon_cell = (xt_ocean{1});
lon_fs = lon_cell(fs_lon_i);
lat_fs = lat(fs_lat_i);
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

%% separating into seperate months

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



%% initial imaging

[Z,LON] = meshgrid(-z_cell,lon_fs);

% % January
% figure
% hold on
% contourf(LON,Z,jan_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,jan_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% 
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean January EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/01_RARE.jpg');
% 
% % February
% figure
% hold on
% contourf(LON,Z,feb_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,feb_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean February EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/02_RARE.jpg');
% 
% % March
% figure
% hold on
% contourf(LON,Z,mar_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,mar_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean March EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/03_RARE.jpg');
% 
% 
% % April
% figure
% hold on
% contourf(LON,Z,apr_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,apr_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean April EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/04_RARE.jpg');
% 
% % May
% figure
% hold on
% contourf(LON,Z,may_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,may_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean May EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/05_RARE.jpg');
% 
% % June
% figure
% hold on
% contourf(LON,Z,jun_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,jun_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean June EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/06_RARE.jpg');
% 
% % July
% figure
% hold on
% contourf(LON,Z,jul_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,jul_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean July EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/07_RARE.jpg');
% 
% % August
% figure
% hold on
% contourf(LON,Z,aug_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,aug_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean August EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/08_RARE.jpg');
% 
% 
% % September
% figure
% hold on
% contourf(LON,Z,sep_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,sep_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean September EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/09_RARE.jpg');
% 
% % October
% figure
% hold on
% contourf(LON,Z,oct_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,oct_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean October EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/10_RARE.jpg');
% 
% % November
% figure
% hold on
% contourf(LON,Z,nov_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,nov_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean November EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/11_RARE.jpg');
% 
% % December
% figure
% hold on
% contourf(LON,Z,dec_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,dec_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean December EGS Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/12_RARE.jpg');
% 
%% Overall mean

% for K = 1:num_files
% tot_v_cells{K} = squeeze(mean(v{K}(:,:,:,:),4,'omitnan'));
% tot_v_mean = cell2mat(tot_v_cells);
% end
% 
% [r,c] = size(tot_v_mean);
% nlay  = 17;
% tot_v_mat   = permute(reshape(tot_v_mean,[r,c/nlay,nlay]),[1,2,3]);
% tot_v_mean = mean(tot_v_mat,3,'omitnan');
% 
% figure
% hold on
% contourf(LON,Z,tot_v_mean,8,'linestyle','none');
% shading interp
% %set(gca, 'XDir','reverse')
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% [C,d] = contour(LON,Z,tot_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('EGS Mean Flow (RARE)','FontSize',16);
% 
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_mean_figs/total_avg_RARE.jpg');

% %% calculating monthly anomalies
% 
% % January
% jan_v_anom = NaN(size(jan_v_mat));
% for K = 1:num_files
% jan_v_anom(:,:,K) = jan_v_mat(:,:,K)-tot_v_mean;
% end
% jan_v_anom = mean(jan_v_anom,3,'omitnan');
% 
% % February
% feb_v_anom = NaN(size(feb_v_mat));
% for K = 1:num_files
% feb_v_anom(:,:,K) = feb_v_mat(:,:,K)-tot_v_mean;
% end
% feb_v_anom = mean(feb_v_anom,3,'omitnan');
% 
% % March
% mar_v_anom = NaN(size(mar_v_mat));
% for K = 1:num_files
% mar_v_anom(:,:,K) = mar_v_mat(:,:,K)-tot_v_mean;
% end
% mar_v_anom = mean(mar_v_anom,3,'omitnan');
% 
% % April
% apr_v_anom = NaN(size(apr_v_mat));
% for K = 1:num_files
% apr_v_anom(:,:,K) = apr_v_mat(:,:,K)-tot_v_mean;
% end
% apr_v_anom = mean(apr_v_anom,3,'omitnan');
% 
% % May
% may_v_anom = NaN(size(may_v_mat));
% for K = 1:num_files
% may_v_anom(:,:,K) = may_v_mat(:,:,K)-tot_v_mean;
% end
% may_v_anom = mean(may_v_anom,3,'omitnan');
% 
% % June
% jun_v_anom = NaN(size(jun_v_mat));
% for K = 1:num_files
% jun_v_anom(:,:,K) = jun_v_mat(:,:,K)-tot_v_mean;
% end
% jun_v_anom = mean(jun_v_anom,3,'omitnan');
% 
% % July
% jul_v_anom = NaN(size(jul_v_mat));
% for K = 1:num_files
% jul_v_anom(:,:,K) = jul_v_mat(:,:,K)-tot_v_mean;
% end
% jul_v_anom = mean(jul_v_anom,3,'omitnan');
% 
% % August
% aug_v_anom = NaN(size(aug_v_mat));
% for K = 1:num_files
% aug_v_anom(:,:,K) = aug_v_mat(:,:,K)-tot_v_mean;
% end
% aug_v_anom = mean(aug_v_anom,3,'omitnan');
% 
% % September
% sep_v_anom = NaN(size(sep_v_mat));
% for K = 1:num_files
% sep_v_anom(:,:,K) = sep_v_mat(:,:,K)-tot_v_mean;
% end
% sep_v_anom = mean(sep_v_anom,3,'omitnan');
% 
% % October
% oct_v_anom = NaN(size(oct_v_mat));
% for K = 1:num_files
% oct_v_anom(:,:,K) = oct_v_mat(:,:,K)-tot_v_mean;
% end
% oct_v_anom = mean(oct_v_anom,3,'omitnan');
% 
% % November
% nov_v_anom = NaN(size(nov_v_mat));
% for K = 1:num_files
% nov_v_anom(:,:,K) = nov_v_mat(:,:,K)-tot_v_mean;
% end
% nov_v_anom = mean(nov_v_anom,3,'omitnan');
% 
% % December
% dec_v_anom = NaN(size(dec_v_mat));
% for K = 1:num_files
% dec_v_anom(:,:,K) = dec_v_mat(:,:,K)-tot_v_mean;
% end
% dec_v_anom = mean(dec_v_anom,3,'omitnan');


%% Imaging the monthly anomalies

% % January
% figure
% hold on
% contourf(LON,Z,jan_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,jan_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous January EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/01a_RARE.jpg');
% 
% 
% % February
% figure
% hold on
% contourf(LON,Z,feb_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,feb_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous February EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/02a_RARE.jpg');
% 
% % March
% figure
% hold on
% contourf(LON,Z,mar_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,mar_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous March EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/03a_RARE.jpg');
% 
% % April
% figure
% hold on
% contourf(LON,Z,apr_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,apr_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous April EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/04a_RARE.jpg');
% 
% % May
% figure
% hold on
% contourf(LON,Z,may_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,may_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous May EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/05a_RARE.jpg');
% 
% % June
% figure
% hold on
% contourf(LON,Z,jun_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,jun_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous June EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/06a_RARE.jpg');
% 
% % July
% figure
% hold on
% contourf(LON,Z,jul_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,jul_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous July EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/07a_RARE.jpg');
% 
% % August
% figure
% hold on
% contourf(LON,Z,aug_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,aug_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous August EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/08a_RARE.jpg');
% 
% % September
% figure
% hold on
% contourf(LON,Z,sep_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,sep_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous September EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/09a_RARE.jpg');
% 
% % October
% figure
% hold on
% contourf(LON,Z,oct_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,oct_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous October EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/10a_RARE.jpg');
% 
% % November
% figure
% hold on
% contourf(LON,Z,nov_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,nov_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous November EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/11a_RARE.jpg');
% 
% % December
% figure
% hold on
% contourf(LON,Z,dec_v_anom,8,'LineStyle','none');
% shading interp
% x = colorbar;
% ylim([-2650 0]);
% clim([-0.04 0.04]);
% [C,h] = contour(LON,Z,dec_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous December EGS Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_anom_figs/12a_RARE.jpg');


%% Creating a 2D spatial plot of EGS Flow

% Establishing 2D variables
addpath('/Users/saravianco/Documents/Research/MATLAB/m_map/');

% determining lat and lon limits for figure
latS_map = 59;
latN_map = 85;
lonW_map = -45;
lonE_map = 15;

% finding these lat/lon ranges in lat/lon variables
lat2D_i = find(lat>59 & lat<85);
lon2D_i = find(lon>-45 & lon<15);
start2D = [lon2D_i(1) lat2D_i(1) 1];
count2D = [length(lon2D_i) length(lat2D_i) 12];

start3D = [lon2D_i(1) lat2D_i(1) 1 1];
count3D = [length(lon2D_i) length(lat2D_i) length(z_star) 12];

%% Loading in bathymetry data
% % Load 1-minute ETOPO data
% 
% addpath('/Users/saravianco/Documents/Research/MATLAB/ETOPO/')
% info = ncinfo('ETOPO_2022_v1_60s_N90W180_bed.nc');
% 
% lonbath = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','lon');
% latbath = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','lat');
% 
% idxSy = find(abs(latbath-latS_map)==nanmin(abs(latbath-latS_map)),1,'first');
% idxNy = find(abs(latbath-latN_map)==nanmin(abs(latbath-latN_map)),1,'first');
% idxWx = find(abs(lonbath-lonW_map)==nanmin(abs(lonbath-lonW_map)),1,'first');
% idxEx = find(abs(lonbath-lonE_map)==nanmin(abs(lonbath-lonE_map)),1,'first');
% 
% % setting the bounds of the bathymetry data to the area that we're
% % interested in
% lonbath = lonbath(idxWx:idxEx);
% latbath = latbath(idxSy:idxNy);
% nlons = length(lonbath);
% nlats = length(latbath);
% bath = -ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','z',[idxWx idxSy],[nlons nlats]);

%% reading in ML depths for all months/years

% for K = 1:num_files
%   this_file = filenames{K};
%   ml_rho{K} = ncread(this_file,'mlp',start2D,count2D);
% end
% 
% % creating a new v variable that's 600x260x38x12 (for the mixed layer)
% for K = 1:num_files
%   this_file = filenames{K};
%   v_ml{K} = ncread(this_file,'v',start3D,count3D);
% end

 %% z 3D
% 
% z_a = repmat(z_star,1,360,60); % use repmat to create a 3x10x10 copy
% z_3D = permute(z_a,[3 2 1]);

%% separating into seperate months and averaging
[LON2D,LAT2D] = meshgrid(lon_cell(lon2D_i),lat_cell(lat2D_i));

% % January 
% % calculating the mixed layer
% for K = 1:num_files
% jan_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,1),3,'omitnan'));
% end
% jan_ml_mean = cell2mat(jan_ml_cells);
% 
% [r,c] = size(jan_ml_mean);
% nlay  = 17;
% jan_ml_mat   = permute(reshape(jan_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% jan_ml_mean = mean(jan_ml_mat,3,'omitnan');
% jan_ml_overall = mean(jan_ml_mat,'all','omitnan');
% jan_ml_mean(isnan(jan_ml_mean)) = jan_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% jan_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,1),4,'omitnan'));
% end
% 
% jan_vml_mean = cell2mat(jan_vml_cells);
% [r,c,l] = size(jan_vml_mean);
% nlay  = 17;
% jan_vml_mat = permute(reshape(jan_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% jan_vml_3D = mean(jan_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(jan_ml_mean(:,1));
% nlats = length(jan_ml_mean(1,:));
% jan_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% jan_ml_idx(i,j) = find(z_star>jan_ml_mean(i,j),1,'first');
%     end
% end
% 
% jan_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         jan_vml_mean(i,j) = mean(jan_vml_3D(i,j,1:jan_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % February
% % calculating the mixed layer
% for K = 1:num_files
% feb_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,2),3,'omitnan'));
% end
% feb_ml_mean = cell2mat(feb_ml_cells);
% 
% [r,c] = size(feb_ml_mean);
% nlay  = 17;
% feb_ml_mat   = permute(reshape(feb_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% feb_ml_mean = mean(feb_ml_mat,3,'omitnan');
% feb_ml_overall = mean(feb_ml_mat,'all','omitnan');
% feb_ml_mean(isnan(feb_ml_mean)) = feb_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% feb_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,2),4,'omitnan'));
% end
% 
% feb_vml_mean = cell2mat(feb_vml_cells);
% [r,c,l] = size(feb_vml_mean);
% nlay  = 17;
% feb_vml_mat = permute(reshape(feb_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% feb_vml_3D = mean(feb_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(feb_ml_mean(:,1));
% nlats = length(feb_ml_mean(1,:));
% feb_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% feb_ml_idx(i,j) = find(z_star>feb_ml_mean(i,j),1,'first');
%     end
% end
% 
% feb_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         feb_vml_mean(i,j) = mean(feb_vml_3D(i,j,1:feb_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % March
% % calculating the mixed layer
% for K = 1:num_files
% mar_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,3),3,'omitnan'));
% end
% mar_ml_mean = cell2mat(mar_ml_cells);
% 
% [r,c] = size(mar_ml_mean);
% nlay  = 17;
% mar_ml_mat   = permute(reshape(mar_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% mar_ml_mean = mean(mar_ml_mat,3,'omitnan');
% mar_ml_overall = mean(mar_ml_mat,'all','omitnan');
% mar_ml_mean(isnan(mar_ml_mean)) = mar_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% mar_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,3),4,'omitnan'));
% end
% 
% mar_vml_mean = cell2mat(mar_vml_cells);
% [r,c,l] = size(mar_vml_mean);
% nlay  = 17;
% mar_vml_mat = permute(reshape(mar_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% mar_vml_3D = mean(mar_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(mar_ml_mean(:,1));
% nlats = length(mar_ml_mean(1,:));
% mar_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% mar_ml_idx(i,j) = find(z_star>mar_ml_mean(i,j),1,'first');
%     end
% end
% 
% mar_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         mar_vml_mean(i,j) = mean(mar_vml_3D(i,j,1:mar_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % April
% % calculating the mixed layer
% for K = 1:num_files
% apr_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,4),3,'omitnan'));
% end
% apr_ml_mean = cell2mat(apr_ml_cells);
% 
% [r,c] = size(apr_ml_mean);
% nlay  = 17;
% apr_ml_mat   = permute(reshape(apr_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% apr_ml_mean = mean(apr_ml_mat,3,'omitnan');
% apr_ml_overall = mean(apr_ml_mat,'all','omitnan');
% apr_ml_mean(isnan(apr_ml_mean)) = apr_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% apr_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,4),4,'omitnan'));
% end
% 
% apr_vml_mean = cell2mat(apr_vml_cells);
% [r,c,l] = size(apr_vml_mean);
% nlay  = 17;
% apr_vml_mat = permute(reshape(apr_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% apr_vml_3D = mean(apr_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(apr_ml_mean(:,1));
% nlats = length(apr_ml_mean(1,:));
% apr_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% apr_ml_idx(i,j) = find(z_star>apr_ml_mean(i,j),1,'first');
%     end
% end
% 
% apr_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         apr_vml_mean(i,j) = mean(apr_vml_3D(i,j,1:apr_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % May
% % calculating the mixed layer
% for K = 1:num_files
% may_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,5),3,'omitnan'));
% end
% may_ml_mean = cell2mat(may_ml_cells);
% 
% [r,c] = size(may_ml_mean);
% nlay  = 17;
% may_ml_mat   = permute(reshape(may_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% may_ml_mean = mean(may_ml_mat,3,'omitnan');
% may_ml_overall = mean(may_ml_mat,'all','omitnan');
% may_ml_mean(isnan(may_ml_mean)) = may_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% may_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,5),4,'omitnan'));
% end
% 
% may_vml_mean = cell2mat(may_vml_cells);
% [r,c,l] = size(may_vml_mean);
% nlay  = 17;
% may_vml_mat = permute(reshape(may_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% may_vml_3D = mean(may_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(may_ml_mean(:,1));
% nlats = length(may_ml_mean(1,:));
% may_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% may_ml_idx(i,j) = find(z_star>may_ml_mean(i,j),1,'first');
%     end
% end
% 
% may_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         may_vml_mean(i,j) = mean(may_vml_3D(i,j,1:may_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % June
% % calculating the mixed layer
% for K = 1:num_files
% jun_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,6),3,'omitnan'));
% end
% jun_ml_mean = cell2mat(jun_ml_cells);
% 
% [r,c] = size(jun_ml_mean);
% nlay  = 17;
% jun_ml_mat   = permute(reshape(jun_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% jun_ml_mean = mean(jun_ml_mat,3,'omitnan');
% jun_ml_overall = mean(jun_ml_mat,'all','omitnan');
% jun_ml_mean(isnan(jun_ml_mean)) = jun_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% jun_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,6),4,'omitnan'));
% end
% 
% jun_vml_mean = cell2mat(jun_vml_cells);
% [r,c,l] = size(jun_vml_mean);
% nlay  = 17;
% jun_vml_mat = permute(reshape(jun_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% jun_vml_3D = mean(jun_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(jun_ml_mean(:,1));
% nlats = length(jun_ml_mean(1,:));
% jun_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% jun_ml_idx(i,j) = find(z_star>jun_ml_mean(i,j),1,'first');
%     end
% end
% 
% jun_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         jun_vml_mean(i,j) = mean(jun_vml_3D(i,j,1:jun_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % July
% % calculating the mixed layer
% for K = 1:num_files
% jul_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,7),3,'omitnan'));
% end
% jul_ml_mean = cell2mat(jul_ml_cells);
% 
% [r,c] = size(jul_ml_mean);
% nlay  = 17;
% jul_ml_mat   = permute(reshape(jul_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% jul_ml_mean = mean(jul_ml_mat,3,'omitnan');
% jul_ml_overall = mean(jul_ml_mat,'all','omitnan');
% jul_ml_mean(isnan(jul_ml_mean)) = jul_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% jul_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,7),4,'omitnan'));
% end
% 
% jul_vml_mean = cell2mat(jul_vml_cells);
% [r,c,l] = size(jul_vml_mean);
% nlay  = 17;
% jul_vml_mat = permute(reshape(jul_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% jul_vml_3D = mean(jul_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(jul_ml_mean(:,1));
% nlats = length(jul_ml_mean(1,:));
% jul_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% jul_ml_idx(i,j) = find(z_star>jul_ml_mean(i,j),1,'first');
%     end
% end
% 
% jul_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         jul_vml_mean(i,j) = mean(jul_vml_3D(i,j,1:jul_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % August
% % calculating the mixed layer
% for K = 1:num_files
% aug_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,8),3,'omitnan'));
% end
% aug_ml_mean = cell2mat(aug_ml_cells);
% 
% [r,c] = size(aug_ml_mean);
% nlay  = 17;
% aug_ml_mat   = permute(reshape(aug_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% aug_ml_mean = mean(aug_ml_mat,3,'omitnan');
% aug_ml_overall = mean(aug_ml_mat,'all','omitnan');
% aug_ml_mean(isnan(aug_ml_mean)) = aug_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% aug_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,8),4,'omitnan'));
% end
% 
% aug_vml_mean = cell2mat(aug_vml_cells);
% [r,c,l] = size(aug_vml_mean);
% nlay  = 17;
% aug_vml_mat = permute(reshape(aug_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% aug_vml_3D = mean(aug_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(aug_ml_mean(:,1));
% nlats = length(aug_ml_mean(1,:));
% aug_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% aug_ml_idx(i,j) = find(z_star>aug_ml_mean(i,j),1,'first');
%     end
% end
% 
% aug_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         aug_vml_mean(i,j) = mean(aug_vml_3D(i,j,1:aug_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % September
% % calculating the mixed layer
% for K = 1:num_files
% sep_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,9),3,'omitnan'));
% end
% sep_ml_mean = cell2mat(sep_ml_cells);
% 
% [r,c] = size(sep_ml_mean);
% nlay  = 17;
% sep_ml_mat   = permute(reshape(sep_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% sep_ml_mean = mean(sep_ml_mat,3,'omitnan');
% sep_ml_overall = mean(sep_ml_mat,'all','omitnan');
% sep_ml_mean(isnan(sep_ml_mean)) = sep_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% sep_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,9),4,'omitnan'));
% end
% 
% sep_vml_mean = cell2mat(sep_vml_cells);
% [r,c,l] = size(sep_vml_mean);
% nlay  = 17;
% sep_vml_mat = permute(reshape(sep_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% sep_vml_3D = mean(sep_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(sep_ml_mean(:,1));
% nlats = length(sep_ml_mean(1,:));
% sep_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% sep_ml_idx(i,j) = find(z_star>sep_ml_mean(i,j),1,'first');
%     end
% end
% 
% sep_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         sep_vml_mean(i,j) = mean(sep_vml_3D(i,j,1:sep_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % October
% % calculating the mixed layer
% for K = 1:num_files
% oct_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,10),3,'omitnan'));
% end
% oct_ml_mean = cell2mat(oct_ml_cells);
% 
% [r,c] = size(oct_ml_mean);
% nlay  = 17;
% oct_ml_mat   = permute(reshape(oct_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% oct_ml_mean = mean(oct_ml_mat,3,'omitnan');
% oct_ml_overall = mean(oct_ml_mat,'all','omitnan');
% oct_ml_mean(isnan(oct_ml_mean)) = oct_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% oct_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,10),4,'omitnan'));
% end
% 
% oct_vml_mean = cell2mat(oct_vml_cells);
% [r,c,l] = size(oct_vml_mean);
% nlay  = 17;
% oct_vml_mat = permute(reshape(oct_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% oct_vml_3D = mean(oct_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(oct_ml_mean(:,1));
% nlats = length(oct_ml_mean(1,:));
% oct_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% oct_ml_idx(i,j) = find(z_star>oct_ml_mean(i,j),1,'first');
%     end
% end
% 
% oct_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         oct_vml_mean(i,j) = mean(oct_vml_3D(i,j,1:oct_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % November
% % calculating the mixed layer
% for K = 1:num_files
% nov_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,11),3,'omitnan'));
% end
% nov_ml_mean = cell2mat(nov_ml_cells);
% 
% [r,c] = size(nov_ml_mean);
% nlay  = 17;
% nov_ml_mat   = permute(reshape(nov_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% nov_ml_mean = mean(nov_ml_mat,3,'omitnan');
% nov_ml_overall = mean(nov_ml_mat,'all','omitnan');
% nov_ml_mean(isnan(nov_ml_mean)) = nov_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% nov_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,11),4,'omitnan'));
% end
% 
% nov_vml_mean = cell2mat(nov_vml_cells);
% [r,c,l] = size(nov_vml_mean);
% nlay  = 17;
% nov_vml_mat = permute(reshape(nov_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% nov_vml_3D = mean(nov_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(nov_ml_mean(:,1));
% nlats = length(nov_ml_mean(1,:));
% nov_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% nov_ml_idx(i,j) = find(z_star>nov_ml_mean(i,j),1,'first');
%     end
% end
% 
% nov_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         nov_vml_mean(i,j) = mean(nov_vml_3D(i,j,1:nov_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% % December
% % calculating the mixed layer
% for K = 1:num_files
% dec_ml_cells{K} = squeeze(mean(ml_rho{K}(:,:,12),3,'omitnan'));
% end
% dec_ml_mean = cell2mat(dec_ml_cells);
% 
% [r,c] = size(dec_ml_mean);
% nlay  = 17;
% dec_ml_mat   = permute(reshape(dec_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% dec_ml_mean = mean(dec_ml_mat,3,'omitnan');
% dec_ml_overall = mean(dec_ml_mat,'all','omitnan');
% dec_ml_mean(isnan(dec_ml_mean)) = dec_ml_overall;
% 
% % surface velocity
% for K = 1:num_files
% dec_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,12),4,'omitnan'));
% end
% 
% dec_vml_mean = cell2mat(dec_vml_cells);
% [r,c,l] = size(dec_vml_mean);
% nlay  = 17;
% dec_vml_mat = permute(reshape(dec_vml_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
% 
% dec_vml_3D = mean(dec_vml_mat,4,'omitnan');
% 
% % determining the values of v within the mixed layer (for each month)
% nlons = length(dec_ml_mean(:,1));
% nlats = length(dec_ml_mean(1,:));
% dec_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% dec_ml_idx(i,j) = find(z_star>dec_ml_mean(i,j),1,'first');
%     end
% end
% 
% dec_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         dec_vml_mean(i,j) = mean(dec_vml_3D(i,j,1:dec_ml_idx(i,j)-1),3,'omitnan');
%     end
% end

 %% imaging the spatial flow
% addpath('/Users/saravianco/Documents/Research/MATLAB/m_map/data/gshhg-bin-2.3.7/')
% 
% % January
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,jan_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% title('EGS January Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_spatial/01_spatial_RARE.jpg');
% 
% % February
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,feb_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% title('EGS February Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_spatial/02_spatial_RARE.jpg');
% 
% 
% % March
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,mar_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% title('EGS March Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_monthly_spatial/03_spatial_RARE.jpg');
% 
% 
% % April
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,apr_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% 
% title('EGS April Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/04_spatial_RARE.jpg']);
% 
% % May
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,may_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% 
% title('EGS May Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/05_spatial_RARE.jpg']);
% 
% % June
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,jun_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% 
% title('EGS June Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/06_spatial_RARE.jpg']);
% 
% % July
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,jul_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% 
% title('EGS July Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/07_spatial_RARE.jpg']);
% 
% % August
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,aug_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% 
% title('EGS August Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/08_spatial_RARE.jpg']);
% 
% % September
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,sep_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% 
% title('EGS September Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/09_spatial_RARE.jpg']);
% 
% % October
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,oct_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% 
% title('EGS October Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/10_spatial_RARE.jpg']);
% 
% % November
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,nov_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% 
% title('EGS November Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/11_spatial_RARE.jpg']);
% 
% % December
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,dec_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V [m/s]");
% 
% title('EGS December Mean ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/12_spatial_RARE.jpg']);
% 
% %% Overall Mean Flow
% 
% % determining the mean flow across all time at each depth
% for K = 1:num_files
% tot_vml_cells{K} = squeeze(mean(v_ml{K}(:,:,:,:),4,'omitnan'));
% end
% tot_vml_3D = cell2mat(tot_vml_cells);
% 
% [r,c,l] = size(tot_vml_3D);
% nlay  = 17;
% tot_vml_mat   = permute(reshape(tot_vml_3D,[r,c/nlay,l,nlay]),[1,2,3,4]);
% tot_vml_3D = mean(tot_vml_mat,4,'omitnan');
% 
% % determining the average mixed layer across all time
% for K = 1:num_files
% tot_ml_cells{K} = squeeze(mean(ml_rho{K},3,'omitnan'));
% end
% tot_ml_mean = cell2mat(tot_ml_cells);
% 
% [r,c] = size(tot_ml_mean);
% nlay  = 17;
% tot_ml_mat   = permute(reshape(tot_ml_mean,[r,c/nlay,nlay]),[1,2,3]);
% tot_ml_mean = mean(tot_ml_mat,3,'omitnan');
% tot_ml_overall = mean(tot_ml_mat,'all','omitnan');
% tot_ml_mean(isnan(tot_ml_mean)) = tot_ml_overall;
% 
% % determining the ml average velocity across all years
% nlons = length(tot_ml_mean(:,1));
% nlats = length(tot_ml_mean(1,:));
% tot_ml_idx = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
% tot_ml_idx(i,j) = find(z_star>tot_ml_mean(i,j),1,'first');
%     end
% end
% 
% tot_vml_mean = NaN(nlons,nlats);
% 
% for i = 1:nlons
%     for j = 1:nlats
%         tot_vml_mean(i,j) = mean(tot_vml_3D(i,j,1:tot_ml_idx(i,j)-1),3,'omitnan');
%     end
% end
% 
% %% Mean Flow figure
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,tot_vml_mean');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','k');
% m_plot([-13 0],[81.45 81.45],'linewi',2,'color','r');
% m_plot([-18 -0],[80.05 80.05],'linewi',2,'color','r');
% m_plot([-18 0],[78.85 78.85],'linewi',2,'color','r');
% m_plot([-19 0],[76.95 76.95],'linewi',2,'color','r');
% m_plot([-20 -8],[73.95 73.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2);
% clim([-0.3 0.2]);
% ylabel(x,"V' [m/s]");
% title('EGS ML Flow (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_spatial/tot_flow_RARE.jpg']);


%% Monthly Spatial Anomalies

% % January
% jan_vml_year = NaN(size(jan_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         jan_vml_year(i,j,K) = mean(jan_vml_mat(i,j,1:jan_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     jan_vml_anom_3D(:,:,K) = jan_vml_year(:,:,K)-tot_vml_mean;
% end
% jan_vml_anom = mean(jan_vml_anom_3D,3,'omitnan');
% 
% % February
% feb_vml_year = NaN(size(feb_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         feb_vml_year(i,j,K) = mean(feb_vml_mat(i,j,1:feb_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     feb_vml_anom_3D(:,:,K) = feb_vml_year(:,:,K)-tot_vml_mean;
% end
% feb_vml_anom = mean(feb_vml_anom_3D,3,'omitnan');
% 
% % March
% mar_vml_year = NaN(size(mar_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         mar_vml_year(i,j,K) = mean(mar_vml_mat(i,j,1:mar_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     mar_vml_anom_3D(:,:,K) = mar_vml_year(:,:,K)-tot_vml_mean;
% end
% mar_vml_anom = mean(mar_vml_anom_3D,3,'omitnan');
% 
% % April
% apr_vml_year = NaN(size(apr_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         apr_vml_year(i,j,K) = mean(apr_vml_mat(i,j,1:apr_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     apr_vml_anom_3D(:,:,K) = apr_vml_year(:,:,K)-tot_vml_mean;
% end
% apr_vml_anom = mean(apr_vml_anom_3D,3,'omitnan');
% 
% % May
% may_vml_year = NaN(size(may_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         may_vml_year(i,j,K) = mean(may_vml_mat(i,j,1:may_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     may_vml_anom_3D(:,:,K) = may_vml_year(:,:,K)-tot_vml_mean;
% end
% may_vml_anom = mean(may_vml_anom_3D,3,'omitnan');
% 
% % June
% jun_vml_year = NaN(size(jun_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         jun_vml_year(i,j,K) = mean(jun_vml_mat(i,j,1:jun_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     jun_vml_anom_3D(:,:,K) = jun_vml_year(:,:,K)-tot_vml_mean;
% end
% jun_vml_anom = mean(jun_vml_anom_3D,3,'omitnan');
% 
% % July
% jul_vml_year = NaN(size(jul_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         jul_vml_year(i,j,K) = mean(jul_vml_mat(i,j,1:jul_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     jul_vml_anom_3D(:,:,K) = jul_vml_year(:,:,K)-tot_vml_mean;
% end
% jul_vml_anom = mean(jul_vml_anom_3D,3,'omitnan');
% 
% % August
% aug_vml_year = NaN(size(aug_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         aug_vml_year(i,j,K) = mean(aug_vml_mat(i,j,1:aug_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     aug_vml_anom_3D(:,:,K) = aug_vml_year(:,:,K)-tot_vml_mean;
% end
% aug_vml_anom = mean(aug_vml_anom_3D,3,'omitnan');
% 
% % September
% sep_vml_year = NaN(size(sep_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         sep_vml_year(i,j,K) = mean(sep_vml_mat(i,j,1:sep_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     sep_vml_anom_3D(:,:,K) = sep_vml_year(:,:,K)-tot_vml_mean;
% end
% sep_vml_anom = mean(sep_vml_anom_3D,3,'omitnan');
% 
% % October
% oct_vml_year = NaN(size(oct_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         oct_vml_year(i,j,K) = mean(oct_vml_mat(i,j,1:oct_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     oct_vml_anom_3D(:,:,K) = oct_vml_year(:,:,K)-tot_vml_mean;
% end
% oct_vml_anom = mean(oct_vml_anom_3D,3,'omitnan');
% 
% % November
% nov_vml_year = NaN(size(nov_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         nov_vml_year(i,j,K) = mean(nov_vml_mat(i,j,1:nov_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     nov_vml_anom_3D(:,:,K) = nov_vml_year(:,:,K)-tot_vml_mean;
% end
% nov_vml_anom = mean(nov_vml_anom_3D,3,'omitnan');
% 
% % December
% dec_vml_year = NaN(size(dec_v_mat));
% for K = 1:num_files
%     for i = 1:nlons
%         for j = 1:nlats
%         dec_vml_year(i,j,K) = mean(dec_vml_mat(i,j,1:dec_ml_idx(i,j)-1,K),3,'omitnan');
%         end
%     end
% end
% % determining the mean anomalies for the month
% for K = 1:num_files
%     dec_vml_anom_3D(:,:,K) = dec_vml_year(:,:,K)-tot_vml_mean;
% end
% dec_vml_anom = mean(dec_vml_anom_3D,3,'omitnan');


%% spatial anomaly plots

% % January
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,jan_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS January ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/01_sp_anom_RARE.jpg']);
% 
% 
% % February
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,feb_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS February ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/02_sp_anom_RARE.jpg']);
% 
% % March
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,mar_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS March ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/03_sp_anom_RARE.jpg']);
% 
% % April
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,apr_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS April ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/04_sp_anom_RARE.jpg']);
% 
% % May
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,may_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS May ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/05_sp_anom_RARE.jpg']);
% 
% % June
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,jun_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS June ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/06_sp_anom_RARE.jpg']);
% 
% % July
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,jul_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS July ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/07_sp_anom_RARE.jpg']);
% 
% % August
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,aug_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS August ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/08_sp_anom_RARE.jpg']);
% 
% % September
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,sep_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS September ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/09_sp_anom_RARE.jpg']);
% 
% % October
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,oct_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS October ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/10_sp_anom_RARE.jpg']);
% 
% % November
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,nov_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS November ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/11_sp_anom_RARE.jpg']);
% 
% % December
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%  hold on
% m_pcolor(LON2D,LAT2D,dec_vml_anom');
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% x = colorbar;
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     'ytick',60:5:80,'xtick',-45:15:15); m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
%     'ytick',60:5:80,'xtick',-45:15:15);
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');
% 
% set(x,'fontsize',10,'fontname','helvetica',...     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.3:0.1:0.2); clim([-0.3 0.2]);
%     'position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',-0.1:0.05:0.05);
% clim([-0.1 0.05]);
% ylabel(h,"V' [m/s]");
% title('EGS December ML Flow Anomaly (RARE)','FontSize',16);
% saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
%     'RARE_monthly_anom_spatial/12_sp_anom_RARE.jpg']);

%% pressure considerations when defining the EGC
% 
% % Converting epth (z*) to pressure
% addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/");
% p = (gsw_p_from_z(-z_star,lat_fs));
% [SA,in_ocean] = gsw_SA_from_SP(SP,p,long,lat);
% CT = gsw_CT_from_pt(SA,pt);
% rho_CT_exact = gsw_rho_CT_exact(SA,CT,p);

%% Locating the EGC based on literature definitions %%
% in general we are going to be looking at the top portion of the water
% column (maybe around 150 m and shallower), identifying the maximum
% velocity in this area, and then finding where this velocity decays to 20%
% of that maximum value (following Havik et al., 2017)

% defining the 'top layer' => roughly 150m
% ul_idx = find(z_star>150,1,'first');
% 
% % January
% up_150_mean = mean(jan_v_mean(:,1:ul_idx),2);
% core_jan = max(jan_v_mean(:,1:ul_idx),[],'all','ComparisonMethod',"abs");
% 
% outer_lim_jan = 0.3*core_jan;
% out_lim1_idx = find(up_150_mean<=outer_lim_jan,1,'first');
% out_lim2_idx = find(up_150_mean<=outer_lim_jan,1,'last');

% figure
% hold on
% contourf(LON,Z,jan_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,jan_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% clabel(C,d);
% [CC,dd] = contour(LON(out_lim1_idx:end,:),Z(out_lim1_idx:end,:),jan_v_mean(out_lim1_idx:end,:),[core_jan,outer_lim_jan],'-r','linewidth',2);
% clabel(CC,dd);
% 
% xline(-5.25);
% xline(-2.05);
% yline(0);
% yline(-1000);
% % xline(lon_fs(out_lim1_idx),'linewidth',2);
% % xline(lon_fs(60),'linewidth',2);
% 
% ylim([-2650 0]);
% clim([-0.25 0.11]);
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean January EGC (RARE)','FontSize',16);

%% A different approach to Vol Tx: Establishing fixed box for EGC transport
lon_1 = -6.05;
lon_2 = -2.05;
z_1 = 0;
z_2 = -1000;

lon1_idx = find(lon_fs<=-6 & lon_fs>=-6.07);
lon2_idx = find(lon_fs<=-2 & lon_fs>=-2.07);
gridx_diff = lon2_idx-lon1_idx;

lon_km = haversine([lat_fs,lon_1],[lat_fs,lon_2]);
lon_m = lon_km*1000;
dx = lon_m/gridx_diff; % horizontal distance of a gridcell

dz = NaN(23,1);
for i = 1:23
    dz(i) = z_star(i+1)-z_star(i);
end

% [Z_EGC,LON_EGC] = meshgrid(-z_star(1:24),lon_fs(28:60));

% January
jan_egc = jan_v_mean(lon1_idx:lon2_idx,1:23);
% specifying negative values
[row,col] = find(jan_egc>0);
jan_egc(row,col) = 0;
jan_egc(isnan(jan_egc)) = 0;
tx_jan = (sum(dx*jan_egc*dz))/10^6;

jan_egc_sd = std(jan_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_jan = (sum((dx*jan_egc_sd*dz),'omitnan'))/10^6;

% February
feb_egc = feb_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(feb_egc>0);
feb_egc(row,col) = 0;
feb_egc(isnan(feb_egc)) = 0;
tx_feb = (sum(dx*feb_egc*dz))/10^6;

feb_egc_sd = std(feb_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_feb = (sum((dx*feb_egc_sd*dz),'omitnan'))/10^6;

% March
mar_egc = mar_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(mar_egc>0);
mar_egc(row,col) = 0;
mar_egc(isnan(mar_egc)) = 0;
tx_mar = (sum(dx*mar_egc*dz))/10^6;

mar_egc_sd = std(mar_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_mar = (sum((dx*mar_egc_sd*dz),'omitnan'))/10^6;

% April
apr_egc = apr_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(apr_egc>0);
apr_egc(row,col) = 0;
apr_egc(isnan(apr_egc)) = 0;
tx_apr = (sum(dx*apr_egc*dz))/10^6;

apr_egc_sd = std(apr_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_apr = (sum((dx*apr_egc_sd*dz),'omitnan'))/10^6;

% May
may_egc = may_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(may_egc>0);
may_egc(row,col) = 0;
may_egc(isnan(may_egc)) = 0;
tx_may = (sum(dx*may_egc*dz))/10^6;

may_egc_sd = std(may_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_may = (sum((dx*may_egc_sd*dz),'omitnan'))/10^6;


% June
jun_egc = jun_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(jun_egc>0);
jun_egc(row,col) = 0;
jun_egc(isnan(jun_egc)) = 0;
tx_jun = (sum(dx*jun_egc*dz))/10^6;

jun_egc_sd = std(jun_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_jun = (sum((dx*jun_egc_sd*dz),'omitnan'))/10^6;

% July
jul_egc = jul_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(jul_egc>0);
jul_egc(row,col) = 0;
jul_egc(isnan(jul_egc)) = 0;
tx_jul = (sum(dx*jul_egc*dz))/10^6;

jul_egc_sd = std(jul_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_jul = (sum((dx*jul_egc_sd*dz),'omitnan'))/10^6;

% August
aug_egc = aug_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(aug_egc>0);
aug_egc(row,col) = 0;
aug_egc(isnan(aug_egc)) = 0;
tx_aug = (sum(dx*aug_egc*dz))/10^6;

aug_egc_sd = std(aug_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_aug = (sum((dx*aug_egc_sd*dz),'omitnan'))/10^6;

% September
sep_egc = sep_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(sep_egc>0);
sep_egc(row,col) = 0;
sep_egc(isnan(sep_egc)) = 0;
tx_sep = (sum(dx*sep_egc*dz))/10^6;

sep_egc_sd = std(sep_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_sep = (sum((dx*sep_egc_sd*dz),'omitnan'))/10^6;

% October
oct_egc = oct_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(oct_egc>0);
oct_egc(row,col) = 0;
oct_egc(isnan(oct_egc)) = 0;
tx_oct = (sum(dx*oct_egc*dz))/10^6;

oct_egc_sd = std(oct_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_oct = (sum((dx*oct_egc_sd*dz),'omitnan'))/10^6;


% November
nov_egc = nov_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(nov_egc>0);
nov_egc(row,col) = 0;
nov_egc(isnan(nov_egc)) = 0;
tx_nov = (sum(dx*nov_egc*dz))/10^6;

nov_egc_sd = std(nov_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_nov = (sum((dx*nov_egc_sd*dz),'omitnan'))/10^6;

% December
dec_egc = dec_v_mean(lon1_idx:lon2_idx,1:23);
[row,col] = find(dec_egc>0);
dec_egc(row,col) = 0;
dec_egc(isnan(dec_egc)) = 0;
tx_dec = (sum(dx*dec_egc*dz))/10^6;

dec_egc_sd = std(dec_v_mat(lon1_idx:lon2_idx,1:23,:),0,3,'omitnan');
err_dec = (sum((dx*dec_egc_sd*dz),'omitnan'))/10^6;

%
monthly_tx = [tx_jan,tx_feb,tx_mar,tx_apr,tx_may,tx_jun,tx_jul,tx_aug,...
    tx_sep,tx_oct,tx_nov,tx_dec];
monthly_err = [err_jan,err_feb,err_mar,err_apr,err_may,err_jun,err_jul,err_aug,...
    err_sep,err_oct,err_nov,err_dec];
month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

figure
hold on
plot(monthly_tx,'-k','linewi',2)
errorbar(monthly_tx,monthly_err,'.r','LineWidth',1)
set(gca,'Ydir','reverse')
%ylim([-5.3 -2.7])
xlim([1 12])
xticks(1:12);
xticklabels(month)

ylabel('Transport [Sv]')
title('EGC Average Monthly Volumetric Transport (RARE)')
saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
    'RARE_vol_tx/EGC_vol_tx.jpg']);

 % Entire section
lon_1 = lon_fs(1);
lon_2 = lon_fs(end);

lon1_idx = 1;
lon2_idx = 60;
gridx_diff = lon2_idx-lon1_idx;

lon_km = haversine([lat_fs,lon_1],[lat_fs,lon_2]);
lon_m = lon_km*1000;
dx = lon_m/gridx_diff; % horizontal distance of a gridcell

dz = NaN(size(z_star));
for i = 1:(length(z_star(:,1))-1)
    dz(i) = z_star(i+1)-z_star(i);
end

dz(isnan(dz)) = 0;

% January
jan_v_mean(isnan(jan_v_mean)) = 0;
[row,col] = find(jan_v_mean>0);
jan_v_mean(row,col) = 0;
tx_jan = (sum(dx*jan_v_mean*dz))/10^6;

% February
feb_v_mean(isnan(feb_v_mean)) = 0;
[row,col] = find(feb_v_mean>0);
feb_v_mean(row,col) = 0;
tx_feb = (sum(dx*feb_v_mean*dz))/10^6;

% March
mar_v_mean(isnan(mar_v_mean)) = 0;
[row,col] = find(mar_v_mean>0);
mar_v_mean(row,col) = 0;
tx_mar = (sum(dx*mar_v_mean*dz))/10^6;

% April
apr_v_mean(isnan(apr_v_mean)) = 0;
[row,col] = find(apr_v_mean>0);
apr_v_mean(row,col) = 0;
tx_apr = (sum(dx*apr_v_mean*dz))/10^6;

% May
may_v_mean(isnan(may_v_mean)) = 0;
[row,col] = find(may_v_mean>0);
may_v_mean(row,col) = 0;
tx_may = (sum(dx*may_v_mean*dz))/10^6;

% June
jun_v_mean(isnan(jun_v_mean)) = 0;
[row,col] = find(jun_v_mean>0);
jun_v_mean(row,col) = 0;
tx_jun = (sum(dx*jun_v_mean*dz))/10^6;

% July
jul_v_mean(isnan(jul_v_mean)) = 0;
[row,col] = find(jul_v_mean>0);
jul_v_mean(row,col) = 0;
tx_jul = (sum(dx*jul_v_mean*dz))/10^6;

% August
aug_v_mean(isnan(aug_v_mean)) = 0;
[row,col] = find(aug_v_mean>0);
aug_v_mean(row,col) = 0;
tx_aug = (sum(dx*aug_v_mean*dz))/10^6;

% September
sep_v_mean(isnan(sep_v_mean)) = 0;
[row,col] = find(sep_v_mean>0);
sep_v_mean(row,col) = 0;
tx_sep = (sum(dx*sep_v_mean*dz))/10^6;

% October
oct_v_mean(isnan(oct_v_mean)) = 0;
[row,col] = find(oct_v_mean>0);
oct_v_mean(row,col) = 0;
tx_oct = (sum(dx*oct_v_mean*dz))/10^6;

% November
nov_v_mean(isnan(nov_v_mean)) = 0;
[row,col] = find(nov_v_mean>0);
nov_v_mean(row,col) = 0;
tx_nov = (sum(dx*nov_v_mean*dz))/10^6;

% December
dec_v_mean(isnan(dec_v_mean)) = 0;
[row,col] = find(dec_v_mean>0);
dec_v_mean(row,col) = 0;
tx_dec = (sum(dx*dec_v_mean*dz))/10^6;

% timeseries
monthly_tx = [tx_jan,tx_feb,tx_mar,tx_apr,tx_may,tx_jun,tx_jul,tx_aug,...
    tx_sep,tx_oct,tx_nov,tx_dec];
month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

figure
plot(monthly_tx,'-k','linewi',2)
set(gca,'Ydir','reverse')
%ylim([-10 -3])
xlim([1 12])
xticks(1:12);
xticklabels(month)

ylabel('Transport [Sv]')
title('Full Section Average Monthly Volumetric Transport:(RARE)')
saveas(gcf,['/Users/saravianco/Documents/Research/MATLAB/RARE/' ...
    'RARE_vol_tx/vol_tx_monthly_fullsct.jpg']);