

%% Preliminaries
close all
clear all

cd('/Users/saravianco/Documents/Research/MATLAB/GLORYS/')

addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/html/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/library/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/thermodynamics_from_t/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/pdf/");

%%
load('TS Diagrams/ts_vars_GLORYS.mat')

%% Interpolating the depths to avoid averaging

z_interp = [0:5:600];

clear SA_interp
clear CT_interp
% Interpolate the matrix onto the new grid
SA_interp = NaN(479,143,121,204);
CT_interp = NaN(479,143,121,204);

for jj = 1:204
SA_interp(:,:,:,jj) = interp3(lat3D,lon3D,z,SA(:,:,:,jj),lat3D,lon3D,z_interp);
CT_interp(:,:,:,jj) = interp3(lat3D,lon3D,z,CT(:,:,:,jj),lat3D,lon3D,z_interp);
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

xlim([29.5 36])
ylim([-2.3 7])

xlabel('Absolute Salinity [PSU]','FontSize',18)
ylabel('Conservative Temperature [^oC]','FontSize',18)
ylabel(x,'Depth [m]','FontSize',18)
title('2003-2019 GLORYS Averages','FontSize',20)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/TS Diagrams/ts_avg204.jpg']);


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
sigma0_5 = (gsw_rho(SA_ref,CT_ref,500))-1000;

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

xlabel('Absolute Salinity [g/kg]','FontSize',25)
ylabel('Conservative Temperature [^oC]','FontSize',25)
ylabel(x,'Relative Percentage of Sample [%]','FontSize',25)

% title('2003-2019 GLORYS Averages: Ekman Layer (0-50 m) ','FontSize',25)
% title('2003-2019 GLORYS Averages: Sub-Ekman Shelf Layer (50-200 m) ','FontSize',20)
% title('2003-2019 GLORYS Averages: Shelfbreak Layer (200-600 m) ','FontSize',20)

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
'GLORYS/TS Diagrams/ts_ekman.jpg']);
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/TS Diagrams/ts_shelf.jpg']);
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
% 'GLORYS/TS Diagrams/ts_shelfbreak.jpg']);