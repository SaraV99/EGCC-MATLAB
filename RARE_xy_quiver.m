% RARE Velocity Quiver and Salinity Analysis
    % Evaluating the time-dependent flow at various depths in RARE to
    % determine when the EGCC splits from the EGC and hops onto the shelf
% Created by: Sara Vianco
% Masters Thesis
% 20 July 2023

%% Preliminaries
close all
clear all

cd('/Volumes/One Touch/Research Datasets/RARE/matching_theo/')
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

%% initial info/variable establishment
info = ncinfo("rare1.15.2_mn_ocean_reg_2003.nc");
info2 = ncinfo("rare1.15.2_mn_ocean_reg_2019.nc");

time = ncread('rare1.15.2_mn_ocean_reg_2003.nc','time');
lon = ncread('rare1.15.2_mn_ocean_reg_2003.nc','xt_ocean');
lat = ncread('rare1.15.2_mn_ocean_reg_2003.nc','yt_ocean');
z_star = ncread('rare1.15.2_mn_ocean_reg_2003.nc','sw_ocean');

%% Establishing 3D Variables
% determining lat and lon limits for figure
latS_map = 73;
latN_map = 85;
lonW_map = -35;
lonE_map = 5;

% finding these lat/lon ranges in lat/lon variables
lat3D_i = find(lat>latS_map & lat<latN_map);
lon3D_i = find(lon<lonE_map & lon>lonW_map);
% start2D = [lon3D_i(1) lat3D_i(1) 1];
% count2D = [length(lon3D_i) length(lat3D_i) 12];

start3D = [lon3D_i(1) lat3D_i(1) 1 1];
count3D = [length(lon3D_i) length(lat3D_i) length(z_star) 12];

%%
ncvars =  {'time', 'xt_ocean', 'yt_ocean', 'sw_ocean', 'v','temp','salt','u'};
projectdir = '/Volumes/One Touch/Research Datasets/RARE/matching_theo/';
dinfo = dir( fullfile(projectdir, '*.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );

time = cell(num_files, 1);
xt_ocean = cell(num_files, 1); 
yt_ocean = cell(num_files, 1);
sw_ocean = cell(num_files, 1);
v = cell(num_files, 1);
temp = cell(num_files,1);
salt = cell(num_files,1);

for K = 1:num_files
  this_file = filenames{K};
  time{K} = ncread(this_file, ncvars{1});
  xt_ocean{K} = ncread(this_file, ncvars{2});
  yt_ocean{K} = ncread(this_file, ncvars{3});
  sw_ocean{K} = ncread(this_file, ncvars{4});
  v{K} = ncread(this_file, ncvars{5},start3D,count3D);
end

for K = 1:num_files
  ptemp{K} = ncread(this_file, ncvars{6},start3D,count3D);
  Sp{K} = ncread(this_file, ncvars{7},start3D,count3D);
  u{K} = ncread(this_file, ncvars{8},start3D,count3D);

end

lon_cell = (xt_ocean{1});
lon3D = lon_cell(lon3D_i);
lat3D = lat(lat3D_i);
lat_cell = (yt_ocean{1});


z_cell = sw_ocean{1};

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

for K = 1:num_files
   
date_num{K} = datenum(time{K} + base);
date_string{K} = datestr(date_num{K},'mm/dd/yy');
dv{K} = datevec(date_string{K},formatIn);

end

%% separating into seperate months for V

% January
for K = 1:num_files
jan_v_cells{K} = squeeze(mean(v{K}(:,:,:,1),4,'omitnan'));
jan_v_mean = cell2mat(jan_v_cells);
end

[r,c,l] = size(jan_v_mean);
nlay  = 17;
jan_v_mat  = permute(reshape(jan_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
jan_v_mean = mean(jan_v_mat,4,'omitnan');

% February
for K = 1:num_files
feb_v_cells{K} = squeeze(mean(v{K}(:,:,:,2),4,'omitnan'));
feb_v_mean = cell2mat(feb_v_cells);
end

[r,c,l] = size(feb_v_mean);
nlay  = 17;
feb_v_mat  = permute(reshape(feb_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
feb_v_mean = mean(feb_v_mat,4,'omitnan');

% March
for K = 1:num_files
mar_v_cells{K} = squeeze(mean(v{K}(:,:,:,3),4,'omitnan'));
mar_v_mean = cell2mat(mar_v_cells);
end

[r,c,l] = size(mar_v_mean);
nlay  = 17;
mar_v_mat   = permute(reshape(mar_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
mar_v_mean = mean(mar_v_mat,4,'omitnan');

% April
for K = 1:num_files
apr_v_cells{K} = squeeze(mean(v{K}(:,:,:,4),4,'omitnan'));
apr_v_mean = cell2mat(apr_v_cells);
end

[r,c,l] = size(apr_v_mean);
nlay  = 17;
apr_v_mat   = permute(reshape(apr_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
apr_v_mean = mean(apr_v_mat,4,'omitnan');

% May
for K = 1:num_files
may_v_cells{K} = squeeze(mean(v{K}(:,:,:,5),4,'omitnan'));
may_v_mean = cell2mat(may_v_cells);
end

[r,c,l] = size(may_v_mean);
nlay  = 17;
may_v_mat   = permute(reshape(may_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
may_v_mean = mean(may_v_mat,4,'omitnan');

% June
for K = 1:num_files
jun_v_cells{K} = squeeze(mean(v{K}(:,:,:,6),4,'omitnan'));
jun_v_mean = cell2mat(jun_v_cells);
end

[r,c,l] = size(jun_v_mean);
nlay  = 17;
jun_v_mat   = permute(reshape(jun_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
jun_v_mean = mean(jun_v_mat,4,'omitnan');

% July
for K = 1:num_files
jul_v_cells{K} = squeeze(mean(v{K}(:,:,:,7),4,'omitnan'));
jul_v_mean = cell2mat(jul_v_cells);
end

[r,c,l] = size(jul_v_mean);
nlay  = 17;
jul_v_mat   = permute(reshape(jul_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
jul_v_mean = mean(jul_v_mat,4,'omitnan');

% August
for K = 1:num_files
aug_v_cells{K} = squeeze(mean(v{K}(:,:,:,8),4,'omitnan'));
aug_v_mean = cell2mat(aug_v_cells);
end

[r,c,l] = size(aug_v_mean);
nlay  = 17;
aug_v_mat   = permute(reshape(aug_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
aug_v_mean = mean(aug_v_mat,4,'omitnan');

% September
for K = 1:num_files
sep_v_cells{K} = squeeze(mean(v{K}(:,:,:,9),4,'omitnan'));
sep_v_mean = cell2mat(sep_v_cells);
end

[r,c,l] = size(sep_v_mean);
nlay  = 17;
sep_v_mat   = permute(reshape(sep_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
sep_v_mean = mean(sep_v_mat,4,'omitnan');

% October
for K = 1:num_files
oct_v_cells{K} = squeeze(mean(v{K}(:,:,:,10),4,'omitnan'));
oct_v_mean = cell2mat(oct_v_cells);
end

[r,c,l] = size(oct_v_mean);
nlay  = 17;
oct_v_mat   = permute(reshape(oct_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
oct_v_mean = mean(oct_v_mat,4,'omitnan');

% November
for K = 1:num_files
nov_v_cells{K} = squeeze(mean(v{K}(:,:,:,11),4,'omitnan'));
nov_v_mean = cell2mat(nov_v_cells);
end

[r,c,l] = size(nov_v_mean);
nlay  = 17;
nov_v_mat   = permute(reshape(nov_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
nov_v_mean = mean(nov_v_mat,4,'omitnan');

% December
for K = 1:num_files
dec_v_cells{K} = squeeze(mean(v{K}(:,:,:,12),4,'omitnan'));
dec_v_mean = cell2mat(dec_v_cells);
end

[r,c,l] = size(dec_v_mean);
nlay  = 17;
dec_v_mat   = permute(reshape(dec_v_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
dec_v_mean = mean(dec_v_mat,4,'omitnan');

%% separating into seperate months for U

% January
for K = 1:num_files
jan_u_cells{K} = squeeze(mean(u{K}(:,:,:,1),4,'omitnan'));
jan_u_mean = cell2mat(jan_u_cells);
end

[r,c,l] = size(jan_u_mean);
nlay  = 17;
jan_u_mat  = permute(reshape(jan_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
jan_u_mean = mean(jan_u_mat,4,'omitnan');

% February
for K = 1:num_files
feb_u_cells{K} = squeeze(mean(u{K}(:,:,:,2),4,'omitnan'));
feb_u_mean = cell2mat(feb_u_cells);
end

[r,c,l] = size(feb_u_mean);
nlay  = 17;
feb_u_mat  = permute(reshape(feb_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
feb_u_mean = mean(feb_u_mat,4,'omitnan');

% March
for K = 1:num_files
mar_u_cells{K} = squeeze(mean(u{K}(:,:,:,3),4,'omitnan'));
mar_u_mean = cell2mat(mar_u_cells);
end

[r,c,l] = size(mar_u_mean);
nlay  = 17;
mar_u_mat   = permute(reshape(mar_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
mar_u_mean = mean(mar_u_mat,4,'omitnan');

% April
for K = 1:num_files
apr_u_cells{K} = squeeze(mean(u{K}(:,:,:,4),4,'omitnan'));
apr_u_mean = cell2mat(apr_u_cells);
end

[r,c,l] = size(apr_u_mean);
nlay  = 17;
apr_u_mat   = permute(reshape(apr_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
apr_u_mean = mean(apr_u_mat,4,'omitnan');

% May
for K = 1:num_files
may_u_cells{K} = squeeze(mean(u{K}(:,:,:,5),4,'omitnan'));
may_u_mean = cell2mat(may_u_cells);
end

[r,c,l] = size(may_u_mean);
nlay  = 17;
may_u_mat   = permute(reshape(may_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
may_u_mean = mean(may_u_mat,4,'omitnan');

% June
for K = 1:num_files
jun_u_cells{K} = squeeze(mean(u{K}(:,:,:,6),4,'omitnan'));
jun_u_mean = cell2mat(jun_u_cells);
end

[r,c,l] = size(jun_u_mean);
nlay  = 17;
jun_u_mat   = permute(reshape(jun_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
jun_u_mean = mean(jun_u_mat,4,'omitnan');

% July
for K = 1:num_files
jul_u_cells{K} = squeeze(mean(u{K}(:,:,:,7),4,'omitnan'));
jul_u_mean = cell2mat(jul_u_cells);
end

[r,c,l] = size(jul_u_mean);
nlay  = 17;
jul_u_mat   = permute(reshape(jul_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
jul_u_mean = mean(jul_u_mat,4,'omitnan');

% August
for K = 1:num_files
aug_u_cells{K} = squeeze(mean(u{K}(:,:,:,8),4,'omitnan'));
aug_u_mean = cell2mat(aug_u_cells);
end

[r,c,l] = size(aug_u_mean);
nlay  = 17;
aug_u_mat   = permute(reshape(aug_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
aug_u_mean = mean(aug_u_mat,4,'omitnan');

% September
for K = 1:num_files
sep_u_cells{K} = squeeze(mean(u{K}(:,:,:,9),4,'omitnan'));
sep_u_mean = cell2mat(sep_u_cells);
end

[r,c,l] = size(sep_u_mean);
nlay  = 17;
sep_u_mat   = permute(reshape(sep_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
sep_u_mean = mean(sep_u_mat,4,'omitnan');

% October
for K = 1:num_files
oct_u_cells{K} = squeeze(mean(u{K}(:,:,:,10),4,'omitnan'));
oct_u_mean = cell2mat(oct_u_cells);
end

[r,c,l] = size(oct_u_mean);
nlay  = 17;
oct_u_mat   = permute(reshape(oct_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
oct_u_mean = mean(oct_u_mat,4,'omitnan');

% November
for K = 1:num_files
nov_u_cells{K} = squeeze(mean(u{K}(:,:,:,11),4,'omitnan'));
nov_u_mean = cell2mat(nov_u_cells);
end

[r,c,l] = size(nov_u_mean);
nlay  = 17;
nov_u_mat   = permute(reshape(nov_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
nov_u_mean = mean(nov_u_mat,4,'omitnan');

% December
for K = 1:num_files
dec_u_cells{K} = squeeze(mean(u{K}(:,:,:,12),4,'omitnan'));
dec_u_mean = cell2mat(dec_u_cells);
end

[r,c,l] = size(dec_u_mean);
nlay  = 17;
dec_u_mat   = permute(reshape(dec_u_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
dec_u_mean = mean(dec_u_mat,4,'omitnan');

%% separating into seperate months for Salinity

z_a = repmat(z_star,1,120,400); % use repmat to create a 3x10x10 copy
z_3D = permute(z_a,[3 2 1]);

for j = 1:38
p(:,:,j) = gsw_p_from_z(-z_3D(:,:,j),lat3D);
end

% January
for K = 1:num_files
jan_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,1),4,'omitnan'));
jan_s_mean = cell2mat(jan_s_cells);
end

[r,c,l] = size(jan_s_mean);
nlay  = 17;
jan_s_mat  = permute(reshape(jan_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);

for K = 1:num_files
    for j = 1:38
    [jan_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(jan_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
jan_sa_mean = mean(jan_sa_mat,4,'omitnan');
    end
end


% February
for K = 1:num_files
feb_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,2),4,'omitnan'));
feb_s_mean = cell2mat(feb_s_cells);
end

[r,c,l] = size(feb_s_mean);
nlay  = 17;
feb_s_mat  = permute(reshape(feb_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
feb_sa_mean = mean(feb_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [feb_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(feb_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
feb_sa_mean = mean(feb_sa_mat,4,'omitnan');
    end
end

% March
for K = 1:num_files
mar_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,3),4,'omitnan'));
mar_s_mean = cell2mat(mar_s_cells);
end

[r,c,l] = size(mar_s_mean);
nlay  = 17;
mar_s_mat   = permute(reshape(mar_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
mar_sa_mean = mean(mar_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [mar_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(mar_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
mar_sa_mean = mean(mar_sa_mat,4,'omitnan');
    end
end

% April
for K = 1:num_files
apr_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,4),4,'omitnan'));
apr_s_mean = cell2mat(apr_s_cells);
end

[r,c,l] = size(apr_s_mean);
nlay  = 17;
apr_s_mat   = permute(reshape(apr_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
apr_sa_mean = mean(apr_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [apr_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(apr_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
apr_sa_mean = mean(apr_sa_mat,4,'omitnan');
    end
end

% May
for K = 1:num_files
may_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,5),4,'omitnan'));
may_s_mean = cell2mat(may_s_cells);
end

[r,c,l] = size(may_s_mean);
nlay  = 17;
may_s_mat   = permute(reshape(may_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
may_sa_mean = mean(may_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [may_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(may_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
may_sa_mean = mean(may_sa_mat,4,'omitnan');
    end
end

% June
for K = 1:num_files
jun_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,6),4,'omitnan'));
jun_s_mean = cell2mat(jun_s_cells);
end

[r,c,l] = size(jun_s_mean);
nlay  = 17;
jun_s_mat   = permute(reshape(jun_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
jun_sa_mean = mean(jun_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [jun_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(jun_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
jun_sa_mean = mean(jun_sa_mat,4,'omitnan');
    end
end

% July
for K = 1:num_files
jul_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,7),4,'omitnan'));
jul_s_mean = cell2mat(jul_s_cells);
end

[r,c,l] = size(jul_s_mean);
nlay  = 17;
jul_s_mat   = permute(reshape(jul_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
jul_sa_mean = mean(jul_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [jul_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(jul_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
jul_sa_mean = mean(jul_sa_mat,4,'omitnan');
    end
end

% August
for K = 1:num_files
aug_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,8),4,'omitnan'));
aug_s_mean = cell2mat(aug_s_cells);
end

[r,c,l] = size(aug_s_mean);
nlay  = 17;
aug_s_mat   = permute(reshape(aug_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
aug_sa_mean = mean(aug_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [aug_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(aug_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
aug_sa_mean = mean(aug_sa_mat,4,'omitnan');
    end
end

% September
for K = 1:num_files
sep_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,9),4,'omitnan'));
sep_s_mean = cell2mat(sep_s_cells);
end

[r,c,l] = size(sep_s_mean);
nlay  = 17;
sep_s_mat   = permute(reshape(sep_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
sep_sa_mean = mean(sep_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [sep_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(sep_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
sep_sa_mean = mean(sep_sa_mat,4,'omitnan');
    end
end

% October
for K = 1:num_files
oct_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,10),4,'omitnan'));
oct_s_mean = cell2mat(oct_s_cells);
end

[r,c,l] = size(oct_s_mean);
nlay  = 17;
oct_s_mat   = permute(reshape(oct_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
oct_sa_mean = mean(oct_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [oct_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(oct_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
oct_sa_mean = mean(oct_sa_mat,4,'omitnan');
    end
end

% November
for K = 1:num_files
nov_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,11),4,'omitnan'));
nov_s_mean = cell2mat(nov_s_cells);
end

[r,c,l] = size(nov_s_mean);
nlay  = 17;
nov_s_mat   = permute(reshape(nov_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
nov_sa_mean = mean(nov_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [nov_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(nov_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
nov_sa_mean = mean(nov_sa_mat,4,'omitnan');
    end
end

% December
for K = 1:num_files
dec_s_cells{K} = squeeze(mean(Sp{K}(:,:,:,12),4,'omitnan'));
dec_s_mean = cell2mat(dec_s_cells);
end

[r,c,l] = size(dec_s_mean);
nlay  = 17;
dec_s_mat   = permute(reshape(dec_s_mean,[r,c/nlay,l,nlay]),[1,2,3,4]);
dec_sa_mean = mean(dec_s_mat,4,'omitnan');

for K = 1:num_files
    for j = 1:38
    [dec_sa_mat(:,:,j,K), in_ocean] = gsw_SA_from_SP(dec_s_mat(:,:,j,K),p(:,:,j),lon3D,lat3D);
dec_sa_mean = mean(dec_sa_mat,4,'omitnan');
    end
end

%% Flow and Salinity in the Uper 46m (rounding to 50)
[LON,LAT] = meshgrid(lon3D,lat3D);


%%
jan_v_2014 = v{2}(:,:,1,1);
jan_u_2014 = u{2}(:,:,1,1);
[LON,LAT] = meshgrid(lon3D,lat3D);


figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
cmocean('deep')

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
m_quiver(LON,LAT,jan_u_2014',jan_v_2014',2,'k','LineWidth',0.7,'MaxHeadSize',2);
%%

% sa_50 = cat(3,mean(jan_sa_mean(:,:,1:9),3,'omitnan'),mean(feb_sa_mean(:,:,1:9),3,'omitnan'),mean(mar_sa_mean(:,:,1:9),3,'omitnan'),...
%     mean(apr_sa_mean(:,:,1:9),3,'omitnan'),mean(may_sa_mean(:,:,1:9),3,'omitnan'),mean(jun_sa_mean(:,:,1:9),3,'omitnan'),...
%     mean(jul_sa_mean(:,:,1:9),3,'omitnan'),mean(aug_sa_mean(:,:,1:9),3,'omitnan'),mean(sep_sa_mean(:,:,1:9),3,'omitnan'),...
%     mean(oct_sa_mean(:,:,1:9),3,'omitnan'),mean(nov_sa_mean(:,:,1:9),3,'omitnan'),mean(dec_sa_mean(:,:,1:9),3,'omitnan'));
% u_50 = cat(3,mean(jan_u_mean(:,:,1:9),3,'omitnan'),mean(feb_u_mean(:,:,1:9),3,'omitnan'),mean(mar_u_mean(:,:,1:9),3,'omitnan'),...
%     mean(apr_u_mean(:,:,1:9),3,'omitnan'),mean(may_u_mean(:,:,1:9),3,'omitnan'),mean(jun_u_mean(:,:,1:9),3,'omitnan'),...
%     mean(jul_u_mean(:,:,1:9),3,'omitnan'),mean(aug_u_mean(:,:,1:9),3,'omitnan'),mean(sep_u_mean(:,:,1:9),3,'omitnan'),...
%     mean(oct_u_mean(:,:,1:9),3,'omitnan'),mean(nov_u_mean(:,:,1:9),3,'omitnan'),mean(dec_u_mean(:,:,1:9),3,'omitnan'));
% v_50 = cat(3,mean(jan_v_mean(:,:,1:9),3,'omitnan'),mean(feb_v_mean(:,:,1:9),3,'omitnan'),mean(mar_v_mean(:,:,1:9),3,'omitnan'),...
%     mean(apr_v_mean(:,:,1:9),3,'omitnan'),mean(may_v_mean(:,:,1:9),3,'omitnan'),mean(jun_v_mean(:,:,1:9),3,'omitnan'),...
%     mean(jul_v_mean(:,:,1:9),3,'omitnan'),mean(aug_v_mean(:,:,1:9),3,'omitnan'),mean(sep_v_mean(:,:,1:9),3,'omitnan'),...
%     mean(oct_v_mean(:,:,1:9),3,'omitnan'),mean(nov_v_mean(:,:,1:9),3,'omitnan'),mean(dec_v_mean(:,:,1:9),3,'omitnan'));
% 
% %% Flow and Salinity in the Uper 217m (rounding to 200m)
% 
% sa_200 = cat(3,mean(jan_sa_mean(:,:,1:16),3,'omitnan'),mean(feb_sa_mean(:,:,1:16),3,'omitnan'),mean(mar_sa_mean(:,:,1:16),3,'omitnan'),...
%     mean(apr_sa_mean(:,:,1:16),3,'omitnan'),mean(may_sa_mean(:,:,1:16),3,'omitnan'),mean(jun_sa_mean(:,:,1:16),3,'omitnan'),...
%     mean(jul_sa_mean(:,:,1:16),3,'omitnan'),mean(aug_sa_mean(:,:,1:16),3,'omitnan'),mean(sep_sa_mean(:,:,1:16),3,'omitnan'),...
%     mean(oct_sa_mean(:,:,1:16),3,'omitnan'),mean(nov_sa_mean(:,:,1:16),3,'omitnan'),mean(dec_sa_mean(:,:,1:16),3,'omitnan'));
% u_200 = cat(3,mean(jan_u_mean(:,:,1:16),3,'omitnan'),mean(feb_u_mean(:,:,1:16),3,'omitnan'),mean(mar_u_mean(:,:,1:16),3,'omitnan'),...
%     mean(apr_u_mean(:,:,1:16),3,'omitnan'),mean(may_u_mean(:,:,1:16),3,'omitnan'),mean(jun_u_mean(:,:,1:16),3,'omitnan'),...
%     mean(jul_u_mean(:,:,1:16),3,'omitnan'),mean(aug_u_mean(:,:,1:16),3,'omitnan'),mean(sep_u_mean(:,:,1:16),3,'omitnan'),...
%     mean(oct_u_mean(:,:,1:16),3,'omitnan'),mean(nov_u_mean(:,:,1:16),3,'omitnan'),mean(dec_u_mean(:,:,1:16),3,'omitnan'));
% v_200 = cat(3,mean(jan_v_mean(:,:,1:16),3,'omitnan'),mean(feb_v_mean(:,:,1:16),3,'omitnan'),mean(mar_v_mean(:,:,1:16),3,'omitnan'),...
%     mean(apr_v_mean(:,:,1:16),3,'omitnan'),mean(may_v_mean(:,:,1:16),3,'omitnan'),mean(jun_v_mean(:,:,1:16),3,'omitnan'),...
%     mean(jul_v_mean(:,:,1:16),3,'omitnan'),mean(aug_v_mean(:,:,1:16),3,'omitnan'),mean(sep_v_mean(:,:,1:16),3,'omitnan'),...
%     mean(oct_v_mean(:,:,1:16),3,'omitnan'),mean(nov_v_mean(:,:,1:16),3,'omitnan'),mean(dec_v_mean(:,:,1:16),3,'omitnan'));
% 
% %% Flow and Salinity between 150-300m
% 
% sa_150 = cat(3,mean(jan_sa_mean(:,:,14:18),3,'omitnan'),mean(feb_sa_mean(:,:,14:18),3,'omitnan'),mean(mar_sa_mean(:,:,14:18),3,'omitnan'),...
%     mean(apr_sa_mean(:,:,14:18),3,'omitnan'),mean(may_sa_mean(:,:,14:18),3,'omitnan'),mean(jun_sa_mean(:,:,14:18),3,'omitnan'),...
%     mean(jul_sa_mean(:,:,14:18),3,'omitnan'),mean(aug_sa_mean(:,:,14:18),3,'omitnan'),mean(sep_sa_mean(:,:,14:18),3,'omitnan'),...
%     mean(oct_sa_mean(:,:,14:18),3,'omitnan'),mean(nov_sa_mean(:,:,14:18),3,'omitnan'),mean(dec_sa_mean(:,:,14:18),3,'omitnan'));
% u_150 = cat(3,mean(jan_u_mean(:,:,14:18),3,'omitnan'),mean(feb_u_mean(:,:,14:18),3,'omitnan'),mean(mar_u_mean(:,:,14:18),3,'omitnan'),...
%     mean(apr_u_mean(:,:,14:18),3,'omitnan'),mean(may_u_mean(:,:,14:18),3,'omitnan'),mean(jun_u_mean(:,:,14:18),3,'omitnan'),...
%     mean(jul_u_mean(:,:,14:18),3,'omitnan'),mean(aug_u_mean(:,:,14:18),3,'omitnan'),mean(sep_u_mean(:,:,14:18),3,'omitnan'),...
%     mean(oct_u_mean(:,:,14:18),3,'omitnan'),mean(nov_u_mean(:,:,14:18),3,'omitnan'),mean(dec_u_mean(:,:,14:18),3,'omitnan'));
% v_150 = cat(3,mean(jan_v_mean(:,:,14:18),3,'omitnan'),mean(feb_v_mean(:,:,14:18),3,'omitnan'),mean(mar_v_mean(:,:,14:18),3,'omitnan'),...
%     mean(apr_v_mean(:,:,14:18),3,'omitnan'),mean(may_v_mean(:,:,14:18),3,'omitnan'),mean(jun_v_mean(:,:,14:18),3,'omitnan'),...
%     mean(jul_v_mean(:,:,14:18),3,'omitnan'),mean(aug_v_mean(:,:,14:18),3,'omitnan'),mean(sep_v_mean(:,:,14:18),3,'omitnan'),...
%     mean(oct_v_mean(:,:,14:18),3,'omitnan'),mean(nov_v_mean(:,:,14:18),3,'omitnan'),mean(dec_v_mean(:,:,14:18),3,'omitnan'));

%% Flow and Salinity top 600m

sa_600 = cat(3,mean(jan_sa_mean(:,:,1:21),3,'omitnan'),mean(feb_sa_mean(:,:,1:21),3,'omitnan'),mean(mar_sa_mean(:,:,1:21),3,'omitnan'),...
    mean(apr_sa_mean(:,:,1:21),3,'omitnan'),mean(may_sa_mean(:,:,1:21),3,'omitnan'),mean(jun_sa_mean(:,:,1:21),3,'omitnan'),...
    mean(jul_sa_mean(:,:,1:21),3,'omitnan'),mean(aug_sa_mean(:,:,1:21),3,'omitnan'),mean(sep_sa_mean(:,:,1:21),3,'omitnan'),...
    mean(oct_sa_mean(:,:,1:21),3,'omitnan'),mean(nov_sa_mean(:,:,1:21),3,'omitnan'),mean(dec_sa_mean(:,:,1:21),3,'omitnan'));
u_600 = cat(3,mean(jan_u_mean(:,:,1:21),3,'omitnan'),mean(feb_u_mean(:,:,1:21),3,'omitnan'),mean(mar_u_mean(:,:,1:21),3,'omitnan'),...
    mean(apr_u_mean(:,:,1:21),3,'omitnan'),mean(may_u_mean(:,:,1:21),3,'omitnan'),mean(jun_u_mean(:,:,1:21),3,'omitnan'),...
    mean(jul_u_mean(:,:,1:21),3,'omitnan'),mean(aug_u_mean(:,:,1:21),3,'omitnan'),mean(sep_u_mean(:,:,1:21),3,'omitnan'),...
    mean(oct_u_mean(:,:,1:21),3,'omitnan'),mean(nov_u_mean(:,:,1:21),3,'omitnan'),mean(dec_u_mean(:,:,1:21),3,'omitnan'));
v_600 = cat(3,mean(jan_v_mean(:,:,1:21),3,'omitnan'),mean(feb_v_mean(:,:,1:21),3,'omitnan'),mean(mar_v_mean(:,:,1:21),3,'omitnan'),...
    mean(apr_v_mean(:,:,1:21),3,'omitnan'),mean(may_v_mean(:,:,1:21),3,'omitnan'),mean(jun_v_mean(:,:,1:21),3,'omitnan'),...
    mean(jul_v_mean(:,:,1:21),3,'omitnan'),mean(aug_v_mean(:,:,1:21),3,'omitnan'),mean(sep_v_mean(:,:,1:21),3,'omitnan'),...
    mean(oct_v_mean(:,:,1:21),3,'omitnan'),mean(nov_v_mean(:,:,1:21),3,'omitnan'),mean(dec_v_mean(:,:,1:21),3,'omitnan'));

speed_600 = NaN(size(u_600));
for ii = 1:12
speed_600(:,:,ii) = sqrt(u_600(:,:,ii).^2+v_600(:,:,ii).^2);
end

%% Creating an animation: Reusing code for all depth levelsa    

cd('/Users/saravianco/Documents/Research/MATLAB/RARE/RARE_quiver/');
% writerObj = VideoWriter('vs_300_quiver_RARE.avi','Uncompressed AVI');
% writerObj.FrameRate = 1;
% open(writerObj);
% 
u_ref = NaN(size(LON));
v_ref = NaN(size(LON));
u_ref(41,50) = 0.25;
v_ref(41,50) = 0;

% fig1 = figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
% hold on

% for ii = 1:12
% fig1 = figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
% clf
% hold on
% m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
% m_pcolor(LON,LAT,sa_300(:,:,ii)');
% cmocean('haline')
% x = colorbar;
% m_quiver(LON(5:3:end,40:5:390),LAT(5:3:end,40:5:390),u_300(40:5:390,5:3:end,ii)',...
%     v_300(40:5:390,5:3:end,ii)',2.5,'k','LineWidth',0.5,'Clipping','on');
% % set(gca,'Clipping','off')
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r'); % plots the FSAOO mooring locations
% 
% m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
% m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
%     'ytick',60:5:80,'xtick',-35:15:5);
% m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','w');
% ylabel(x,'SA [g/kg]')
% set(x,'fontsize',10,'fontname','helvetica','position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',29:1:35);
% clim([29.5 35.5]);
% m_quiver(LON,LAT,u_ref,v_ref,2.5,'k','LineWidth',0.7);
% 
% 
% dim = [0.37 0.4 .045 .06];
% str = '0.25 m/s';
% annotation('textbox',dim,'String',str,'FontSize',13,'FitBoxToText','off','Rotation',-14.5);
% 
% formatSpec = '%s NEGS Velocity and Salinity in the top 50m (RARE)';
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

%% Salinity and Velocity
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(LON,LAT,sa_600_avg');
cmocean('haline')
x = colorbar;
m_quiver(LON(2:1:end,35:7:395),LAT(2:1:end,35:7:395),u_600_avg(35:7:395,2:1:end)',...
    v_600_avg(35:7:395,2:1:end)',2.5,'k','LineWidth',0.5,'Clipping','on');
m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',0.5,'color','w');
ylabel(x,'SA [g/kg]')
set(x,'fontsize',10,'fontname','helvetica','position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',29:1:35);
clim([29.5 35.5]);
m_quiver(LON,LAT,u_ref,v_ref,1,'k','LineWidth',0.7,'MaxHeadSize',2);

dim = [0.37 0.4 .045 .06];
str = '0.25 m/s';
annotation('textbox',dim,'String',str,'FontSize',13,'FitBoxToText','off','Rotation',-14.5);

title('Time Mean Velocity and Salinity: RARE 2003-2019','FontSize',16)
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_quiver/RARE_tmean_v_sal.jpg']);

%% Velocity and Speed
u_300_avg(41,50) = 1;
v_300_avg(41,50) = 0;

figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(LON,LAT,speed_600_avg');
cmocean('speed')
x = colorbar;
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','w');
m_quiver(LON(2:1:end,35:7:395),LAT(2:1:end,35:7:395),u_600_avg(35:7:395,2:1:end)',...
    v_600_avg(35:7:395,2:1:end)',2.5,'k','LineWidth',0.5,'Clipping','on');
m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);

set(x,'fontsize',10,'fontname','helvetica','position',[0.8 0.6 0.01 0.3],'tickdir','out','ytick',0:0.05:0.2);
clim([0 0.2]);
m_quiver(LON,LAT,u_ref,v_ref,1,'k','LineWidth',0.7,'MaxHeadSize',2);

dim = [0.37 0.4 .045 .06];
str = '0.25 m/s';
annotation('textbox',dim,'String',str,'FontSize',13,'FitBoxToText','off','Rotation',-14.5);

ylabel(x,'Speed [m/s]')
title('Time Mean Flow: RARE 2003-2019','FontSize',16)
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_quiver/RARE_tmean_v_spd.jpg']);  
% 
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

%% Identifying Polygon coordinates: Nordostrundingen latitude
% [x,y] = ginput(1);
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
m_plot([polyNW_lon polyNE_lon],[lat_nord lat_nord],'-g','LineWidth',2)
m_plot([polySW_74 polySE_lon],[polySE_lat polySE_lat],'-r','LineWidth',2)
m_plot(polyNE_lon,polyNE_lat,'*w')
m_plot(polySE_lon,polySE_lat,'*w')

% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
title('Shelf Polygon','FontSize',16)
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/GLORYS quiver/shelf_polygon.jpg']);

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

for i = 1:12
v_600m_res(i,:) = interp2(lat3D,lon3D,v_600(:,:,i),lat_600_poly,lon_600_poly,'nearest');
u_600m_res(i,:) = interp2(lat3D,lon3D,u_600(:,:,i),lat_600_poly,lon_600_poly,'nearest');
end
% it worked!

v_600m_res = v_600m_res';
u_600m_res = u_600m_res';

% finding the high-res shelf angle
[x,y]  = m_ll2xy(lon_600_poly, lat_600_poly);

clear dx
clear dy

dx = NaN(1,1605);
dy = NaN(1,1605);
% calculating a simple slope 
for m = 2:(length(lon_600_poly)-1)
dy(m) = (y(m+1)-y(m-1));
dx(m) = (x(m+1)-x(m-1));
end

ang_rough = atan2d(dy,dx);
vec_ang = -ang_rough;
vec_ang(1,1606) = NaN;

v_600_cont = reshape(v_600m_res,[1,(1606*12)]);
u_600_cont = reshape(u_600m_res,[1,(1606*12)]);

U = cat(1,u_600_cont,v_600_cont);

R = NaN(2,2,(length(u_600m_res)));

vec_ang_yr = repmat(vec_ang,1,12);

for jj = 1:length(vec_ang_yr)
R(:,:,jj)=[cosd(vec_ang_yr(jj)), -sind(vec_ang_yr(jj)) ; sind(vec_ang_yr(jj)), cosd(vec_ang_yr(jj))];
U_r(:,jj) = R(:,:,jj)*U(:,jj);
end

u_600r = U_r(1,:);
v_600r = U_r(2,:);

u_600r = reshape(u_600r,[1606,12]);
v_600r = reshape(v_600r,[1606,12]);

u_600r = u_600r(1:end-1,:);
v_600r = v_600r(1:end-1,:);

% averaging over all months and plotting
u_600c_avg = nanmean(u_600m_res,2);
v_600c_avg = nanmean(v_600m_res,2);

%% testing the coordinates in RARE coordinate system
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
cmocean('deep')
b6 = m_plot(lon_600_poly,lat_600_poly,'-r','LineWidth',1);
m_quiver(lon_600_poly(1:100:end),lat_600_poly(1:100:end),u_600c_avg(1:100:end)',...
    v_600c_avg(1:100:end)','k','LineWidth',1,'Clipping','on');
m_quiver(lon_600_poly(1:100:end-1), lat_600_poly(1:100:end-1), dx(1:100:end), dy(1:100:end),...
    '-r','ShowArrowHead','off') 

% m_plot([lon_nord polyNE_lon],[lat_nord lat_nord],'-r','LineWidth',2)
% m_plot([lon_74 polySE_lon],[polySE_lat polySE_lat],'-r','LineWidth',2)
% m_plot(polyNE_lon,polyNE_lat,'*w')
% m_plot(polySE_lon,polySE_lat,'*w')

% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);

title('Velocity Quivers and Corresponding Shelf Tangents (RARE)')

% %%
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

clear lon_km lon_m 

lon_km = NaN(length(lon_600_poly)-1,1);

for jj = 1:(length(lon_600_poly)-1)
lon_km(jj,1) = haversine([lat_600_poly(jj),lon_600_poly(jj)],[lat_600_poly(jj+1),lon_600_poly(jj+1)]); 
% haversign calculates the geometric distance in km between two coordinates
end

% lon_km(end,1) = lon_km(end-1,1);
lon_m = lon_km*1000;

% to calculate onshelf transport, we're going to find where v is negative and
% multiply that against the size of our grid, as well as by 600m, since our
% velocities are averaged over the top 600m. That way our final units will
% be in Sv
vol_tx_600 = NaN(size(v_600r));

onshelf_idx = find(v_600r<0);
v_onshelf = zeros(size(v_600r));
v_onshelf(onshelf_idx) = v_600r(onshelf_idx);

onshelf_vx = 600*v_onshelf.*lon_m;
vol_tx_600(onshelf_idx) = onshelf_vx((onshelf_idx));

for i = 1:12
onshelf_vx_tot(i) = (nansum(onshelf_vx(:,i),1))/10^6;
end

offshelf_idx = find(v_600r>0);
v_offshelf = zeros(size(v_600r));
v_offshelf(offshelf_idx) = v_600r(offshelf_idx);

offshelf_vx = 600*v_offshelf.*lon_m;
vol_tx_600(offshelf_idx) = offshelf_vx((offshelf_idx));

for i = 1:12
offshelf_vx_tot(i) = (nansum(offshelf_vx(:,i),1))/10^6;
end

% averaging across all years
vx_yr_avg = nanmean(vol_tx_600,2)/10^6;

vol_tx_600 = vol_tx_600/10^6;

vol_tx_600 = vol_tx_600(2:end,:);

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
m_scatter(lon_600_poly(1:end-1),lat_600_poly(1:end-1),100,(1000*vx_yr_avg)','Filled');
clim([-30 5])
cmocean('balance','pivot',0)
h = colorbar;
% m_quiver(lon_600_poly,lat_600_poly,u_600_con,v_600_con,1,'k','LineWidth',0.5,'Clipping','on');
title('Average Cross-Shelf Volumetric Transport in the top 600m (RARE)')
ylabel(h,'Volumetric Transport [mSv]')
% m_quiver(LON,LAT,u_ref,v_ref,2.5,'k','LineWidth',0.7);
% 
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/Vx_600_yearly.jpg']);


%% displaying transport on/off the shelf 
% (same code as above, but running for different months
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
%m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
m_scatter(lon_600_poly(1,2:end-1),lat_600_poly(1,2:end-1),100,1000*vol_tx_600(:,12)','Filled');
clim([-30 5])
cmocean('balance','pivot',0)
h = colorbar;
title('Average December Cross-Shelf Volumetric Transport (RARE)')
ylabel(h,'Volumetric Transport [mSv]')

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/vx_600_12.jpg']);

%% visualizing the volumetric transport by plotting against distance

months = {'Jan','Feb','Mar','Apr','May','Jul','Jul','Aug','Sep','Oct','Nov','Dec'};
% establishing a variable that's distance from the northern tip of our
% polygon

dist_nord = cumsum(lon_m(2:end));

vol_tx_600m = 1000*vol_tx_600;

figure
hold on
imagesc([1:12],dist_nord/1000,vol_tx_600m);
set(gca,'Ydir','reverse');
x = colorbar;
clim([-40 40])
colormap(cmocean('balance','pivot',0))

xlim([0.5 12.5])
xticks(1:12);
xticklabels(months)
ylabel(x,'Volumetric Transport [mSv]','FontSize',10)
title('Volumetric Transport along the NEGS (RARE)')
subtitle('Top 600m')
set(gcf,'color','w');
ylabel('Distance from Nordostrundingen [km]','FontSize',10)
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/monthly_vx_03_19.jpg']);

%% using a line plot to image the temporal (not super helpful)

% c = cmocean('thermal',77);
% shelf_mean_vx = mean(vol_tx_600,1);
% sm_vx_smooth = movmean(shelf_mean_vx,6);
% 
% p = polyfit(1:12,shelf_mean_vx,4);
% y1 = polyval(p,1:12);
% 
% y1_sm = smooth(y1,2);
% 
% figure
% hold on
% % for ii = 1:7:77
% % plot([1:12],vol_tx_600(ii,:),'linewi',1,'Color',c(ii,:));
% % end
% plot(1:12, shelf_mean_vx,'k','linewi',2)
% % plot(dist_nord/1000, sm_vx_smooth,'b','linewi',2);
% plot(1:12,y1_sm,'r','linewi',2)
% 
% xlim([1 12])
% xticks(1:12);
% xticklabels(months)
% legend('Mean VT','Polynomial Fit')
% ylim([-0.18 0])
% ylabel('VT [Sv]')
% 
% title('Seasonal Cycle of Cross-shelf Volumetric Transport')
% set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_quiver/vol tx monthly/sea_cyc_line.jpg']);

% Across the shelf, average monthly volumetric transport is all on-shelf.
% Meaning that there must be some dynamic on the shelf that's exporting
% this water that's coming on

%% Looking at the Volumetric Transport for different seasons

% to plot this value, we can look at different seasons, and further down,
% we ca look at indvidual years if we wanted to 

VT_winter = mean(vol_tx_600m(:,[12,1,2]),2);
VT_spring = mean(vol_tx_600m(:,[3,4,5]),2);
VT_summer = mean(vol_tx_600m(:,[6,7,8]),2);
VT_fall = mean(vol_tx_600m(:,[9,10,11]),2);

figure
hold on
plot(dist_nord/1000,VT_winter,'b')
plot(dist_nord/1000,VT_spring,'g')
plot(dist_nord/1000,VT_summer,'m')
plot(dist_nord/1000,VT_fall,'r')
yline(0,'--k','linewidth',1)
ylim([-30 10])

legend('Winter (DJF)','Spring (MAM)','Summer (JJA)','Autumn (SON)','location','southeast')
xlabel('Distance along 600m Isobath [km]')
ylabel('VT [Sv]')
title('Seasonal Volumetric Transport')
set(gcf,'color','w');
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_quiver/vol tx monthly/sea_vt_line.jpg']);

%% looking at the seasonal cycle across all years and months: U
 for K = 1:17
     u_all_600{K} = squeeze(mean(u{K}(:,:,1:21,:),3,'omitnan'));
 end

u_all_mat = cell2mat(u_all_600);

[r,c,l] = size(u_all_mat);
nlay  = 17;
u_all_600_mat = permute(reshape(u_all_mat,[r,c/nlay,l*nlay]),[1,2,3]);

%% looking at the seasonal cycle across all years and months: V
 for K = 1:17
     v_all_600{K} = squeeze(mean(v{K}(:,:,1:21,:),3,'omitnan'));
 end

v_all_mat = cell2mat(v_all_600);

[r,c,l] = size(v_all_mat);
nlay  = 17;
v_all_600_mat = permute(reshape(v_all_mat,[r,c/nlay,l*nlay]),[1,2,3]);

%%
clear v_600_monthly u_600_monthly
for i = 1:204
v_600_monthly(i,:) = interp2(lat3D,lon3D,v_all_600_mat(:,:,i),lat_600_poly,lon_600_poly,'nearest');
u_600_monthly(i,:) = interp2(lat3D,lon3D,u_all_600_mat(:,:,i),lat_600_poly,lon_600_poly,'nearest');
end

u_600_monthly = u_600_monthly';
v_600_monthly = v_600_monthly';

v_600_mcont = reshape(v_600_monthly,[1,(1606*204)]);
u_600_mcont = reshape(u_600_monthly,[1,(1606*204)]);


U_mo = cat(1,u_600_mcont,v_600_mcont);

R_mo = NaN(2,2,(length(u_600_mcont)));

vec_ang_mo = repmat(vec_ang,1,204);

for jj = 1:length(vec_ang_mo)
R_mo(:,:,jj)=[cosd(vec_ang_mo(jj)), -sind(vec_ang_mo(jj)) ; sind(vec_ang_mo(jj)), cosd(vec_ang_mo(jj))];
U_r_mo(:,jj) = R_mo(:,:,jj)*U_mo(:,jj);
end

u_600r_mo = U_r_mo(1,:);
v_600r_mo = U_r_mo(2,:);

u_600r_mo = reshape(u_600r_mo,[1606,204]);
v_600r_mo = reshape(v_600r_mo,[1606,204]);

u_600r_mo = u_600r_mo(2:end,:);
v_600r_mo = v_600r_mo(2:end,:);

%% determining the onshore/off shore velocities

vx_mo = NaN(size(v_600r_mo));

% finding the indicies where the seasonal rotated v velocities are negative
% to diagnose onshelf transport
onshelf_idx = find(v_600r_mo<0);
v_onshelf = zeros(size(v_600r_mo));
v_onshelf(onshelf_idx) = v_600r_mo(onshelf_idx);

onshelf_vx = 600*v_onshelf.*lon_m;
vx_mo(onshelf_idx) = onshelf_vx((onshelf_idx));

for i = 1:204
onshelf_vx_tot_mo(i) = (nansum(onshelf_vx(:,i),1))/10^6;
end

% finding the indicies where the seasonal rotated v velocities are positive
% to diagnose offshelf transport
offshelf_idx = find(v_600r_mo>0);
v_offshelf = zeros(size(v_600r_mo));
v_offshelf(offshelf_idx) = v_600r_mo(offshelf_idx);

offshelf_vx = 600*v_offshelf.*lon_m;
vx_mo(offshelf_idx) = offshelf_vx((offshelf_idx));

for i = 1:204
offshelf_vx_tot_mo(i) = (nansum(offshelf_vx(:,i),1))/10^6;
end
% converting to sverdrups
vol_tx_monthly = vx_mo/10^6;

%% Reconfiguring to have the VT separated by years -> months

% vol_tx_monthly is already converted to sverdrups
[r,c] = size(vol_tx_monthly);
vx_yr_mo = reshape(vol_tx_monthly,[r,c/nlay,nlay]);

%% Imaging this Yearly Volumetric Transport

% calculating different means
% Mean volumetric transport across our entire section for each month of
% each year
vx_yr_SectionAvg = squeeze(mean(vx_yr_mo,1,'omitnan'));
vx_year_tavg = mean(vx_yr_SectionAvg,2);
vx_year_sd = std(vx_yr_SectionAvg);

c = cmocean('thermal',17);

figure
hold on
for K = 1:17
    plot([1:12],1000*vx_yr_SectionAvg(:,K),'Color',c(K,:))
end
    plot([1:12],1000*vx_year_tavg,'-k','LineWidth',2)
legend('2003','2004','2005','2006','2007','2008','2009','2010','2011',...
    '2012','2013','2014','2015','2016','2017','2018','2019','location','eastoutside')

xlim([1 12])
xticks(1:12);
xticklabels(months)
ylim([-20 5])
ylabel('VT [mSv]')
title('Cross-shelf Volumetric Transport for 2003-2019 (RARE)')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/year_seas.jpg']);

% There's one year that's particularly anomolous: 2006 (2012 is also a bit
% low)

%% Volumetric Transport Anomalies

vx_yr_anom = vx_yr_SectionAvg-vx_year_tavg;

figure
hold on
for K = 1:17
    plot([1:12],1000*vx_yr_anom(:,K),'Color',c(K,:))
end
legend('2003','2004','2005','2006','2007','2008','2009','2010','2011',...
    '2012','2013','2014','2015','2016','2017','2018','2019','location','eastoutside')

xlim([1 12])
xticks(1:12);
xticklabels(months)
ylim([-10 10])
ylabel('VT [mSv]')
title('Cross-shelf Volumetric Transport Anomalies for 2003-2019 (RARE)')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/year_seas_anom.jpg']);

%% Creating a colormap across all months/years
% need to adjust the time and the volumetric transport arrays to both
% resemble 77x204

date_mat = cell2mat(date_num);
[r,c] = size(date_mat);
date_num_all = reshape(date_mat,[1,r*c]);
% date_hyp = datetime(date_num_all, 'ConvertFrom', 'datenum', 'Format', 'MM-dd-yyyy');  
% dv_mat = cell2mat(dv);

%%

figure
hold on
imagesc(date_num_all,dist_nord/1000,1000*vol_tx_monthly);
set(gca,'Ydir','reverse');
dateFormat = 11;
datetick('x',dateFormat)
xlim([date_num_all(1,1) date_num_all(end)])
x = colorbar;
clim([-50 50])
colormap(cmocean('balance','pivot',0))
ylabel('Distance from Nordostrundingen [km]','FontSize',10)
ylabel(x,'VT [mSv]')
title('Cross-shelf Volumetric Transport for 2003-2019 (RARE)')
set(gcf,'color','w');
% 
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/vt_all.jpg']);

%% Volumetric Transport Anomalies (x shelf)
% 
% % first looking at the anomalies across all months/years based on the
% % average at each location
vx_test = mean(vol_tx_monthly,2);
vx_anom_1605 = vol_tx_monthly - vx_test;

figure
hold on
imagesc(date_num_all,dist_nord/1000,1000*vx_anom_1605);
set(gca,'Ydir','reverse');
dateFormat = 11;
datetick('x',dateFormat)
xlim([date_num_all(1,1) date_num_all(end)])
x = colorbar;
clim([-30 30])
colormap(cmocean('balance','pivot',0))
ylabel('Distance from Nordostrundingen [km]','FontSize',10)
ylabel(x,"VT' [mSv]")
title('Cross-shelf Volumetric Transport Anomalies for 2003-2019 (RARE)')
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/vt_all_anom.jpg']);

%% Plotting 100/200/300...km along shelfbreak

dist_nord_km = dist_nord/1000;
xq = [0:100:900];
dist_nord_km_highres = interp1(dist_nord_km,dist_nord_km,xq);

for dd = 1:length(xq)
dist_100s_idx(dd) = find(abs(dist_nord_km-xq(dd))==nanmin(abs(dist_nord_km-xq(dd))));
end

xq_string = string(xq);

%%
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
cmocean('deep')

b6 = m_plot(lon_600_poly,lat_600_poly,'-r','LineWidth',2);
m_plot([polyNW_lon polyNE_lon],[lat_nord lat_nord],'-r','LineWidth',2)
m_plot([polySW_74 polySE_lon],[polySE_lat polySE_lat],'-r','LineWidth',2)
% m_plot(polyNE_lon,polyNE_lat,'*w')
% m_plot(polySE_lon,polySE_lat,'*w')
m_scatter(lon_600_poly(dist_100s_idx),lat_600_poly(dist_100s_idx),100,lat_600_poly(dist_100s_idx)','*w');
% m_plot([-8 -2],[78.95 78.95],'linewi',2,'color','r');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
m_text(lon_600_poly(dist_100s_idx)+1, lat_600_poly(dist_100s_idx),xq_string,'color','w','fontweight','bold')

title('Shelf Polygon with 100km Indications','FontSize',16)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/dist_km_sect.jpg']);

%% Cluster Analysis: Creating a Dendrogram 
% the purpose of this is to try and diagnose how many clusters there are
% within the volumetric transport data

tree = linkage(vol_tx_monthly,'average');

figure()
dendrogram(tree,100,'ColorThreshold','default')
% from this dendrogram, I think I'm going to start with either 7 or 8 nodes

%% K-means clustering

% 1 cluster
[idx1,C] = kmeans(vol_tx_monthly(1:end-1,:),1);
[silh1,h] = silhouette(vol_tx_monthly(1:end-1,:),idx1);

cluster1 = mean(silh1); % NaN

% 2 cluster
[idx2,C] = kmeans(vol_tx_monthly(1:end-1,:),2);
[silh2,h] = silhouette(vol_tx_monthly(1:end-1,:),idx2);

cluster2 = mean(silh2); 

% 3 cluster
[idx3,C] = kmeans(vol_tx_monthly(1:end-1,:),3);
[silh3,h] = silhouette(vol_tx_monthly(1:end-1,:),idx3);

cluster3 = mean(silh3); 

% 4 cluster
[idx4,C] = kmeans(vol_tx_monthly(1:end-1,:),4);
[silh4,h] = silhouette(vol_tx_monthly(1:end-1,:),idx4);

cluster4 = mean(silh4); 

% 5 cluster
[idx5,C] = kmeans(vol_tx_monthly(1:end-1,:),5);
[silh5,h] = silhouette(vol_tx_monthly(1:end-1,:),idx5);

cluster5 = mean(silh5); 

% 6 cluster
[idx6,C] = kmeans(vol_tx_monthly(1:end-1,:),6);
[silh6,h] = silhouette(vol_tx_monthly(1:end-1,:),idx6);

cluster6 = mean(silh6); 

% 7 cluster
[idx7,C] = kmeans(vol_tx_monthly(1:end-1,:),7);
[silh7,h] = silhouette(vol_tx_monthly(1:end-1,:),idx7);

cluster7 = mean(silh7); 

% 8 cluster
[idx8,C] = kmeans(vol_tx_monthly(1:end-1,:),8);
[silh8,h] = silhouette(vol_tx_monthly(1:end-1,:),idx8);

cluster8 = mean(silh8); 

% 9 cluster
[idx9,C] = kmeans(vol_tx_monthly(1:end-1,:),9);
[silh9,h] = silhouette(vol_tx_monthly(1:end-1,:),idx9);

cluster9 = mean(silh9); 

% 10 cluster
[idx10,C] = kmeans(vol_tx_monthly(1:end-1,:),10);
[silh10,h] = silhouette(vol_tx_monthly(1:end-1,:),idx10);

cluster10 = mean(silh10); 

% based on a dendrogram, testing only first 10 clusters. A more visual
% cutoff existed around 7-8 but using the mean silhouette value for the top
% 10 clusters allowed for a more objective approach. 

%%
figure
scatter(idx3,dist_nord_km,'k','Filled');
set(gca,'Ydir','reverse');

title('K-Means Clustering Results (3 clusters, RARE)')
xlabel('Cluster Number','FontSize',10)
ylabel('Distance from Nordostrundingen [km]','FontSize',10)

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/idx3.jpg']);

% it looks like there's three prominant clusters here

%% to best analyze the progression of these

% in a two cluster scenario, finding the upper and lower limits of the
% two clusters

% first cluster
idx2_1_i = find(idx2 == 1,1,'first');
idx2_1_f = find(idx2 == 1,1,'last');

% second cluster
idx2_2_i = find(idx2 == 2,1,'first');
idx2_2_f = find(idx2 == 2,1,'last');

% giving a geographical reference
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(lonbath,latbath,bath')
m_contour(lonbath,latbath,bath',[100 350 500 1000],'linew',1,'color','k');
cmocean('deep')

b6 = m_plot(lon_600_poly,lat_600_poly,'-r','LineWidth',2);
m_plot([polyNW_lon polyNE_lon],[lat_nord lat_nord],'-r','LineWidth',2)
m_plot([polySW_74 polySE_lon],[polySE_lat polySE_lat],'-r','LineWidth',2)
% m_plot(polyNE_lon,polyNE_lat,'*w')
% m_plot(polySE_lon,polySE_lat,'*w')
m_scatter(lon_600_poly(dist_100s_idx),lat_600_poly(dist_100s_idx),100,lat_600_poly(dist_100s_idx)','*w');
m_plot([lon_600_poly(idx2_1_i)-2 lon_600_poly(idx2_1_i)+2],[lat_600_poly(idx2_1_i) lat_600_poly(idx2_1_i)],'linewi',2,'color','y');
m_plot([lon_600_poly(idx2_1_f)-2 lon_600_poly(idx2_1_f)+2],[lat_600_poly(idx2_1_f) lat_600_poly(idx2_1_f)],'linewi',2,'color','y');
m_plot([lon_600_poly(idx2_2_i)-2 lon_600_poly(idx2_2_i)+2],[lat_600_poly(idx2_2_i) lat_600_poly(idx2_2_i)],'linewi',2,'color','g');
m_plot([lon_600_poly(idx2_2_f)-2 lon_600_poly(idx2_2_f)+2],[lat_600_poly(idx2_2_f) lat_600_poly(idx2_2_f)],'linewi',2,'color','g');

m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);

title('Cluster Locations along 600m','FontSize',16)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/geo_kmean_sections.jpg']);

%%
figure
hold on
imagesc(date_num_all,dist_nord/1000,1000*vol_tx_monthly);
plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_i) dist_nord_km(idx2_1_i)],'linewi',2,'color','y');
plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_1_f) dist_nord_km(idx2_1_f)],'linewi',2,'color','y');
plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_2_i) dist_nord_km(idx2_2_i)],'linewi',2,'color','g');
plot([date_num_all(1,1) date_num_all(end)],[dist_nord_km(idx2_2_f) dist_nord_km(idx2_2_f)],'linewi',2,'color','g');

set(gca,'Ydir','reverse');
dateFormat = 11;
datetick('x',dateFormat)
xlim([date_num_all(1,1) date_num_all(end)])
x = colorbar;
clim([-30 30])
colormap(cmocean('balance','pivot',0))
ylabel('Distance from Nordostrundingen [km]','FontSize',10)
ylabel(x,'VT [mSv]')
title('Cross-shelf Volumetric Transport for 2003-2019 (RARE)')
set(gcf,'color','w');

% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_quiver/vol tx monthly/heatmap_kmean_sections.jpg']);

%% At a yearly resolution
figure
hold on
imagesc([1:12],dist_nord/1000,1000*vol_tx_600);

set(gca,'Ydir','reverse');
dateFormat = 11;
datetick('x',dateFormat)
% xlim([date_num_all(1,1) date_num_all(end)])
x = colorbar;
% clim([-30 30])
colormap(cmocean('balance','pivot',0))
ylabel('Distance from Nordostrundingen [km]','FontSize',10)
ylabel(x,'VT [mSv]')
title('Cross-shelf Volumetric Transport for 2003-2019 (RARE)')
set(gcf,'color','w');

% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'RARE/RARE_quiver/vol tx monthly/heatmap_kmean_sections.jpg']);

%% performing an ANOVA test on the different sections

% %
% figure
% histogram(vx_g2) % (vx_g1) % (vx_g3)
% ^ each of the three sections were tested for normality and visually
% passed the test
% [h,p,ksstat,cv] = kstest(vx_g1) < using the Kolmogorov-Smirnoff test
% indicates that none of these distributions are normal...

% vx_g1_mean = (mean(vol_tx_monthly(idx3_1_i:idx3_1_f),1))';
% vx_g2_mean = (mean(vol_tx_monthly(idx3_2_i:idx3_2_f),1))';
% vx_g3_mean = (mean(vol_tx_monthly(idx3_3_i:idx3_3_f),1))';

vx_g1 = vol_tx_monthly(idx2_1_i:idx2_1_f,:)';
vx_g2 = vol_tx_monthly(idx2_2_i:idx2_2_f,:)';

vx_group = cat(2,vx_g1,vx_g2);
group = idx2;

[p,tbl,stats] = anova1(vx_group,group);

[p,tbl,stats] = kruskalwallis(vx_group,group);
% According to the unbalanced ANOVA test the different groups do indeed
% have different means, and are thus, distinct from one another. P =
% 7.66e-26.However, when I average across all the timesteps, the p-value
% reduces to 0.22, which causes us to fail to reject the null hypothesis
% (the different groups are the same).

% A Kruskal-Wallis test, which does not assume normalcy, still returns the
% "reject the null hypothesis" output, which suggests that the different
% groups do not come from the same populations.

%% Looking at Auto-Correlation for the volumetric transport
[r,p] = corr(vol_tx_monthly'); % this looks at the linear correlation between 
% the different areas on the shelf as we move through time (for all 17
% years)

dist_nord_km = dist_nord_km(2:end-1);

figure
hold on
imagesc(dist_nord_km,dist_nord_km,r)
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx2_1_i) dist_nord_km(idx2_1_i)],'linewi',2,'color','y');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx2_1_f) dist_nord_km(idx2_1_f)],'linewi',2,'color','y');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx2_2_i) dist_nord_km(idx2_2_i)],'linewi',2,'color','g');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx2_2_f) dist_nord_km(idx2_2_f)],'linewi',2,'color','g');

set(gca,'Ydir','reverse');
set(gca,'Xdir','reverse');
colormap(cmocean('balance','pivot',0))
colorbar
clim([-0.5 1])
ylabel('Distance from Nordostrundingen [km]','FontSize',12)
xlabel('Distance from Nordostrundingen [km]','FontSize',12)

title('VT Correlation along the 600m Isobath (Monthly, RARE)')

set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/heatmap_vx_corr_yr.jpg']);

%%
[r,p] = corr(vol_tx_600'); % this looks at the linear correlation between 
% the different areas on the shelf as we move through time (for all 17
% years)

figure
hold on
imagesc(dist_nord_km,dist_nord_km,r)
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_1_i) dist_nord_km(idx3_1_i)],'linewi',2,'color','y');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_1_f) dist_nord_km(idx3_1_f)],'linewi',2,'color','y');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_2_i) dist_nord_km(idx3_2_i)],'linewi',2,'color','g');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_2_f) dist_nord_km(idx3_2_f)],'linewi',2,'color','g');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_3_i) dist_nord_km(idx3_3_i)],'linewi',2,'color','b');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_3_f) dist_nord_km(idx3_3_f)],'linewi',2,'color','b');

set(gca,'Ydir','reverse');
set(gca,'Xdir','reverse');
colormap(cmocean('balance','pivot',0))
colorbar
clim([-1 1])
ylabel('Distance from Nordostrundingen [km]','FontSize',12)
xlabel('Distance from Nordostrundingen [km]','FontSize',12)

title('VT Correlation along the 600m Isobath (Yearly, RARE)')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/heatmap_vx_corr_yravg.jpg']);

%%
[r,p] = corr(vx_anom_1605');

figure
hold on
imagesc(dist_nord_km,dist_nord_km,r)
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_1_i) dist_nord_km(idx3_1_i)],'linewi',2,'color','y');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_1_f) dist_nord_km(idx3_1_f)],'linewi',2,'color','y');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_2_i) dist_nord_km(idx3_2_i)],'linewi',2,'color','g');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_2_f) dist_nord_km(idx3_2_f)],'linewi',2,'color','g');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_3_i) dist_nord_km(idx3_3_i)],'linewi',2,'color','b');
% plot([dist_nord_km(1,1) dist_nord_km(end)],[dist_nord_km(idx3_3_f) dist_nord_km(idx3_3_f)],'linewi',2,'color','b');

set(gca,'Ydir','reverse');
set(gca,'Xdir','reverse');
colormap(cmocean('balance','pivot',0))
colorbar
% clim([-1 1])
ylabel('Distance from Nordostrundingen [km]','FontSize',12)
xlabel('Distance from Nordostrundingen [km]','FontSize',12)

title('VT Anomaly Correlation along the 600m Isobath (Monthly, RARE)')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/heatmap_vx_corr_anom.jpg']);

%% Constructing a Volumetric Transport Budget
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

%% testing
figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);
hold on
m_proj('lambert','lat',[latS_map latN_map],'lons',[lonW_map lonE_map]);
m_pcolor(lonbath,latbath,bath')
cmocean('deep')
m_contour(lonbath,latbath,bath',[100 200 250 350 500 1000],'linew',1,'color','k');
% m_plot(lon_600,lat_600,'-k','LineWidth',1)
m_gshhs_i('patch',[0.95 0.95 0.95],'edgecolor','k','linew',0.5);
m_grid('box','fancy','tickdir','out','fontsize',12,'fontname','helvetica',...     
    'ytick',60:5:80,'xtick',-35:15:5);
% m_plot(ll_600_poly(:,1),ll_600_poly(:,2),'-g','LineWidth',1)
m_plot([e_74_lon,w_74_lon,e_nord_lon,w_nord_lon],[s_74_lat,s_74_lat,n_nord_lat,n_nord_lat],'*r')
colorbar
clim([-100 5000])
title('Plotting Polygon Corners (RARE)')

%% Creating a Volume Budget for the Single Mean Across the Entire Timespan
%% Southern Edge of the Polygon
% applplying the polygon limits to the v velocities that are flowing out of
% the polygon and further down the shelf. Assessing the "offshore" tranport
% at this latitude to create more of a volume budget. This is a fairly
% preliminary volume budget...just trying it out

v_74 = v_600_avg(w_74i:e_74i,idxlatS);

lon_74_span = lon3D(w_74i:e_74i);
dxS_km = NaN(size(lon_74_span));
for nn = 2:length(lon_74_span)-1
    dxS_km(nn) = haversine([s_74_lat lon_74_span(nn-1)],[s_74_lat lon_74_span(nn+1)]);
end
dxS_m = 1000*dxS_km;

vx_600_S = NaN(size(v_74));
v_oni = find(v_74>0);  % assuming that we're going straight across the 74
% degree latitude line then all the +v will be considered to be "on" the
% shelf. In reality, this flow is going inside of our polygon
v_74_on = zeros(size(v_74));
v_74_on(v_oni) = v_74(v_oni);
onshelf_74 = 600*v_74_on.*dxS_m;

vx_600_S(v_oni) = onshelf_74(v_oni);

v_offi = find(v_74<0);  
v_74_off = zeros(size(v_74));
v_74_off(v_offi) = v_74(v_offi);
offshelf_74 = 600*v_74_off.*dxS_m;

vx_600_S(v_offi) = offshelf_74(v_offi);

%% Northern Edge of the Polygon
% applplying the polygon limits to the v velocities that are flowing out of
% the polygon and further down the shelf. Assessing the "offshore" tranport
% at this latitude to create more of a volume budget. This is a fairly
% preliminary volume budget...just trying it out

v_nord = v_600_avg(w_nordi:e_nordi,idxlatN);

lon_nord_span = lon3D(w_nordi:e_nordi);

dxN_km = NaN(size(lon_nord_span));
for nn = 2:length(lon_nord_span)-1
    dxN_km(nn) = haversine([n_nord_lat lon_nord_span(nn-1)],[n_nord_lat lon_nord_span(nn+1)]);
end
dxN_m = 1000*dxN_km;

vx_600_N = NaN(size(v_nord));
v_oni = find(v_nord<0);  % assuming that we're going straight across the Nord
% latitude line then all the -v will be considered to be "on" the
% shelf. In reality, this flow is going inside of our polygon
v_nord_on = zeros(size(v_nord));
v_nord_on(v_oni) = v_nord(v_oni);
onshelf_nord = 600*v_nord_on.*dxN_m;

vx_600_N(v_oni) = onshelf_nord(v_oni);

v_offi = find(v_nord>0);  
v_nord_off = zeros(size(v_nord));
v_nord_off(v_offi) = v_nord(v_offi);
offshelf_nord = 600*v_nord_off.*dxN_m;

vx_600_N(v_offi) = offshelf_nord(v_offi);

% summing up all the onshore and offshore transport
onshelf_NS_total = (nansum(abs(onshelf_nord),1)+nansum(abs(onshelf_74),1))/10^6;
offshelf_NS_total = (nansum(abs(offshelf_nord),1)+nansum(abs(offshelf_74),1))/10^6;
% this is now in Sverdrups

onshelf_600_avg = abs(nanmean(onshelf_vx_tot));
offshelf_600_avg = abs(nanmean(offshelf_vx_tot));

onshelf_avg = onshelf_NS_total+onshelf_600_avg;
offshelf_avg = offshelf_NS_total+offshelf_600_avg;

% from this analysis, it looks like our volumetric budget is not
% accurate. We have 21.2316Sv of water flowing onto the shelf and only
% 11.9941Sv exiting the shelf. How can we account for the extra 10Sv?

%% Volume Budget for Individual Months
% Basically doing the same thing as above but at on a monthly average
%% Southern Edge of the Polygon: Monthly Resolution

v_74_mo = squeeze(v_600(w_74i:e_74i,idxlatS,:));

vx_600_S_mo = NaN(size(v_74_mo));
v_oni = find(v_74_mo>0); 
v_74_on = zeros(size(v_74_mo));
v_74_on(v_oni) = v_74_mo(v_oni);

for i = 1:12
onshelf_74_mo(:,i) = 600*v_74_on(:,i).*dxS_m;
end

vx_600_S_mo(v_oni) = onshelf_74_mo(v_oni);

v_offi = find(v_74_mo<0);  
v_74_off = zeros(size(v_74_mo));
v_74_off(v_offi) = v_74_mo(v_offi);

for i = 1:12
offshelf_74_mo(:,i) = 600*v_74_off(:,i).*dxS_m;
end 

vx_600_S_mo(v_offi) = offshelf_74_mo(v_offi);
vx_600_S_mo = squeeze(vx_600_S_mo);

%% Northern Edge of the Polygon: Monthly Resolution
v_nord_mo = squeeze(v_600(w_nordi:e_nordi,idxlatN,:));

lon_nord_span = lon3D(w_nordi:e_nordi);

vx_600_N_mo = NaN(size(v_nord_mo));
v_oni = find(v_nord_mo<0);  
v_nord_on = zeros(size(v_nord_mo));
v_nord_on(v_oni) = v_nord_mo(v_oni);

for i = 1:12
onshelf_nord_mo(:,i) = 600*v_nord_on(:,i).*dxN_m;
end 

vx_600_N_mo(v_oni) = onshelf_nord_mo(v_oni);

v_offi = find(v_nord_mo>0);  
v_nord_off = zeros(size(v_nord_mo));
v_nord_off(v_offi) = v_nord_mo(v_offi);

for i = 1:12
offshelf_nord_mo(:,i) = 600*v_nord_off(:,i).*dxN_m;
end

vx_600_N_mo(v_offi) = offshelf_nord_mo(v_offi);
vx_600_N_mo = squeeze(vx_600_N_mo);

offshelf_NS_mo = (nanmean(abs(offshelf_nord_mo),1)+nanmean(abs(offshelf_74_mo),1))/10^6;
onshelf_NS_mo =  (nanmean(abs(onshelf_nord_mo),1)+nanmean(abs(onshelf_74_mo),1))/10^6;

onshelf_avg_mo = abs(onshelf_vx_tot)+abs(onshelf_NS_mo); % makes sense that the 
% outer 600m isobath VT would dominate this term
offshelf_avg_mo = abs(offshelf_vx_tot)+abs(offshelf_NS_mo);

%% plotting the volumetric budget in the top 600m
offshelf_avg_test = mean(offshelf_avg_mo);
onshelf_avg_test = mean(onshelf_avg_mo);

net_VT_mo = onshelf_avg_mo-offshelf_avg_mo;
net_VT_avg = onshelf_avg_test-offshelf_avg_test;

figure
hold on
plot([1:12],-offshelf_avg_mo,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
yline(-offshelf_avg_test,'--','Color',[0.6350 0.0780 0.1840],'LineWidth',1)
plot([1:12],onshelf_avg_mo,'-','Color',[0 0.4470 0.7410],'LineWidth',2)
yline(onshelf_avg_test,'--','Color',[0 0.4470 0.7410],'LineWidth',1)
plot([1:12],net_VT_mo,'-k','LineWidth',2)
yline(net_VT_avg,'--k','LineWidth',1)
legend('Monthly Average Offshelf VT','2003-2019 Average Offshelf VT',...
    'Monthly Average Onshelf VT','2003-2019 Average Onshelf VT',...
    'Net Average Monthly VT','2003-2009 Average Net VT','location','southeast')
ylim([-17 30])
xlim([1 12])
xticks(1:12);
xticklabels(months)
title('Volumetric Transport Budget in the top 600m (RARE)')
ylabel('VT [Sv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/net_VT_12.jpg']);

%% Southern Edge of the Polygon: All Individual Months 

v_74_204 = squeeze(v_all_600_mat(w_74i:e_74i,idxlatS,:));

vx_600_S_204 = NaN(size(v_74_204));
v_oni = find(v_74_204>0); 
v_74_on_204 = zeros(size(v_74_204));
v_74_on_204(v_oni) = v_74_204(v_oni);

onshelf_74_204 = NaN(size(v_74_on_204));
for i = 1:204
onshelf_74_204(:,i) = 600*v_74_on_204(:,i).*dxS_m;
end

vx_600_S_204(v_oni) = onshelf_74_204(v_oni);

v_offi = find(v_74_204<0);  
v_74_off_204 = zeros(size(v_74_204));
v_74_off_204(v_offi) = v_74_204(v_offi);

for i = 1:204
offshelf_74_204(:,i) = 600*v_74_off_204(:,i).*dxS_m;
end 

vx_600_S_204(v_offi) = offshelf_74_204(v_offi);
vx_600_S_204 = squeeze(vx_600_S_204);

%% Northern Edge of the Polygon: All Individual Months
v_nord_204 = squeeze(v_all_600_mat(w_nordi:e_nordi,idxlatN,:));

vx_600_N_204 = NaN(size(v_nord_204));
v_oni = find(v_nord_204<0);  
v_nord_on_204 = zeros(size(v_nord_204));
v_nord_on_204(v_oni) = v_nord_204(v_oni);

for i = 1:204
onshelf_nord_204(:,i) = 600*v_nord_on_204(:,i).*dxN_m;
end 

vx_600_N_204(v_oni) = onshelf_nord_204(v_oni);

v_offi = find(v_nord_204>0);  
v_nord_off_204 = zeros(size(v_nord_204));
v_nord_off_204(v_offi) = v_nord_204(v_offi);

for i = 1:204
offshelf_nord_204(:,i) = 600*v_nord_off_204(:,i).*dxN_m;
end

vx_600_N_204(v_offi) = offshelf_nord_204(v_offi);
vx_600_N_204 = squeeze(vx_600_N_204);

offshelf_NS_204 = (nanmean(abs(offshelf_nord_204))+nanmean(abs(offshelf_74_204)))/10^6;
onshelf_NS_204 =  (nanmean(abs(onshelf_nord_204))+nanmean(abs(onshelf_74_204)))/10^6;

onshelf_avg_204 = abs(onshelf_vx_tot_mo)+abs(onshelf_NS_204);
offshelf_avg_204 = abs(offshelf_vx_tot_mo)+abs(offshelf_NS_204);

%% Plotting the volumetric budget across all 204 months
offshelf204_avg_test = mean(offshelf_avg_204);
onshelf204_avg_test = mean(onshelf_avg_204);

net_VT_204 = onshelf_avg_204-offshelf_avg_204;
net_VT_204avg = onshelf204_avg_test-offshelf204_avg_test;

years = [{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012',...
    '2013','2014','2015','2016','2017','2018','2019'}];

figure
hold on
plot([1:204],-offshelf_avg_204,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
yline(-offshelf204_avg_test,'--','Color',[0.6350 0.0780 0.1840],'LineWidth',1)
plot([1:204],onshelf_avg_204,'-','Color',[0 0.4470 0.7410],'LineWidth',2)
yline(onshelf204_avg_test,'--','Color',[0 0.4470 0.7410],'LineWidth',1)
plot([1:204],net_VT_204,'-k','LineWidth',2)
yline(net_VT_204avg,'--k','LineWidth',1)
legend('Monthly Average Offshelf VT','2003-2019 Average Offshelf VT',...
    'Monthly Average Onshelf VT','2003-2019 Average Onshelf VT',...
    'Net Average Monthly VT','2003-2009 Average Net VT','location','southoutside')
ylim([-10 20])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)
% dateFormat = 11;
% datetick('x',dateFormat)

title('Volumetric Transport Budget in the top 600m (RARE)')
ylabel('VT [Sv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/net_VT_204.jpg']);

%% Plotting the top and bottom volumetric transport seperately
%% beginning with the northern flow
offshelf_nord_avg = abs((nanmean(offshelf_nord_204,1))/10^6);
onshelf_nord_avg = abs((nanmean(onshelf_nord_204,1))/10^6);

figure
hold on
plot([1:204],-1000*offshelf_nord_avg,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot([1:204],1000*onshelf_nord_avg,'-','Color',[0 0.4470 0.7410],'LineWidth',2)

legend('Average Offshelf Vol Tx','Average Onshelf Vol Tx')
% ylim([])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)

title('Volumetric Transport Budget at the Northern Edge of the Polygon (RARE)')
ylabel('Vol Tx [mSv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/net_VT_Nord.jpg']);

%% continuing with the southern flow
offshelf_74_avg = abs((nanmean(offshelf_74_204,1))/10^6);
onshelf_74_avg = abs((nanmean(onshelf_74_204,1))/10^6);

figure
hold on
plot([1:204],-1000*offshelf_74_avg,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot([1:204],1000*onshelf_74_avg,'-','Color',[0 0.4470 0.7410],'LineWidth',2)

legend('Average Offshelf Vol Tx','Average Onshelf Vol Tx')
% ylim([])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)

title('Volumetric Transport Budget at the Southern Edge of the Polygon (RARE)')
ylabel('Vol Tx [mSv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/net_VT_74.jpg']);

%% finishing with the 600m isobath

figure
hold on
plot([1:204],-offshelf_vx_tot_mo,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot([1:204],-onshelf_vx_tot_mo,'-','Color',[0 0.4470 0.7410],'LineWidth',2)

legend('Average Offshelf Vol Tx','Average Onshelf Vol Tx')
% ylim([])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)

title('Volumetric Transport Budget at the Outer Edge of the Polygon (RARE)')
ylabel('Vol Tx [Sv]')

set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/net_VT_outer.jpg']);

%% just to check
% I'm going to add the off shelf and onshelf values together again to see
% if I get the same plot. Or, at least see what the net result is. 

% onshelf
onshelf_test = abs(onshelf_vx_tot_mo)+abs(onshelf_74_avg)+abs(onshelf_nord_avg);
offshelf_test = abs(offshelf_vx_tot_mo)+abs(offshelf_74_avg)+abs(offshelf_nord_avg);
net_test = onshelf_test-offshelf_test;

figure
hold on
plot([1:204],-offshelf_test,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot([1:204],onshelf_test,'-','Color',[0 0.4470 0.7410],'LineWidth',2)
plot([1:204],net_test,'-k','LineWidth',2)

legend('Offshelf','Onshelf','Net')
% ylim([])
xticks(1:12:204);
xlim([1 204])
xticklabels(years)
ylabel('Vol Tx [Sv]')
title('NEGS Polygon Volume Budget')
set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'RARE/RARE_quiver/vol tx monthly/net_VT_polygon.jpg']);

