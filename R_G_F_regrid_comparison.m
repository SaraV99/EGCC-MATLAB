% RARE/GLORYS/FSAOO mooring Regridded Comparison
% Created by: Sara Vianco
% Masters Thesis
% 18 July 2023

%% Preliminaries
close all
clear all

cd('/Users/saravianco/Documents/Research/MATLAB/Regridded Comparison/')
addpath('/Users/saravianco/Documents/Research/MATLAB/theo/regridded plots/')
addpath('/Users/saravianco/Documents/Research/MATLAB/RARE/RARE regridded/')
addpath('/Users/saravianco/Documents/Research/MATLAB/GLORYS/GLORYS_regridded/')
addpath('/Users/saravianco/Documents/Research/MATLAB/github_repo/') % this is for cmocean
addpath('/Users/saravianco/Documents/MATLAB//haversine/') % this is for haversine

% loading all regridded data

load("GLORYS_regrid_vars.mat")
load("RARE_regridded_vars.mat")
load("fsaoo_regrid_vars.mat")
load('bath_xy.mat')

%% Calculating Overall Means

% Velocity
v_2D_FSAOO = mean(v_interp_FSAOO,3,'omitnan');
v_2D_RARE = mean(v_interp_RARE,3,'omitnan');
v_2D_GLORYS = mean(v_interp_GLORYS,3,'omitnan');

% Absolute Salinity
SA_2D_FSAOO = mean(SAi,3,'omitnan');
SA_2D_RARE = mean(SA_interp_RARE,3,'omitnan');
SA_2D_GLORYS = mean(SA_interp_GLORYS,3,'omitnan');

% Conservative Temperature
CT_2D_FSAOO = mean(CTi,3,'omitnan');
CT_2D_RARE = mean(CT_interp_RARE,3,'omitnan');
CT_2D_GLORYS = mean(CT_interp_GLORYS,3,'omitnan');

% Density Anomaly
rhoa_2D_FSAOO = mean(rho_anom_i,3,'omitnan');
rhoa_2D_RARE = mean(rhoa_interp_RARE,3,'omitnan');
rhoa_2D_GLORYS = mean(rhoa_interp_GLORYS,3,'omitnan');

%% Calculating Differences: Overall Mean

% RARE - FSAOO
v_diff_RF = v_2D_RARE-v_2D_FSAOO;
SA_diff_RF = SA_2D_RARE-SA_2D_FSAOO;
CT_diff_RF = CT_2D_RARE-CT_2D_FSAOO;
rhoa_diff_RF = rhoa_2D_RARE-rhoa_2D_FSAOO;

% GLORYS - FSAOO
v_diff_GF = v_2D_GLORYS-v_2D_FSAOO;
SA_diff_GF = SA_2D_GLORYS-SA_2D_FSAOO;
CT_diff_GF = CT_2D_GLORYS-CT_2D_FSAOO;
rhoa_diff_GF = rhoa_2D_GLORYS-rhoa_2D_FSAOO;

% RARE - GLORYS
v_diff_RG = v_2D_RARE-v_2D_GLORYS;
SA_diff_RG = SA_2D_RARE-SA_2D_GLORYS;
CT_diff_RG = CT_2D_RARE-CT_2D_GLORYS;
rhoa_diff_RG = rhoa_2D_RARE-rhoa_2D_GLORYS;

%% Plotting Overall Means %%

% Adding Mooring locations
moorings_x = [-8 -6.5 -5.5 -5 -4 -3 -2];
moorings_y = zeros(1,7);

%%
% FSAOO
figure
hold on
contourf(XQ,YQ,v_2D_FSAOO,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
[C,d] = contour(XQ,YQ,v_2D_FSAOO,[-0.15,-0.1,-0.05,0,0.05],'-k');
clabel(C,d)

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

xlim([-8 -2])
ylim([-1000 0])
clim([-0.15 0.05]);
cmocean('balance','pivot',0)

% title('Overall Mean V: FSAOO','FontSize',16)
title('FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'V [m/s]','FontSize',20);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/velocity/FSAOO.jpg']);

%% GLORYS
figure
hold on
contourf(XQ,YQ,v_2D_GLORYS,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
[C,d] = contour(XQ,YQ,v_2D_GLORYS,[-0.15,-0.1,-0.05,0,0.05],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
clim([-0.15 0.05]);
cmocean('balance','pivot',0)

% title('Overall Mean V: GLORYS','FontSize',16)
title('GLORYS','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'V [m/s]','FontSize',20);
set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/velocity/GLORYS.jpg']);


%% RARE
figure
hold on
contourf(XQ,YQ,v_2D_RARE,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
[C,d] = contour(XQ,YQ,v_2D_RARE,[-0.15,-0.1,-0.05,0,0.05],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
clim([-0.15 0.05]);
cmocean('balance','pivot',0)

% title('Overall Mean V: RARE','FontSize',16)
title('RARE','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'V [m/s]','FontSize',20);
set(gcf,'color','w');
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/velocity/RARE.jpg']);

%% Plotting Velocity Differences %%
%[bathx,bathy] = ginput(20);

%% RARE - FSAOO
figure
hold on
contourf(XQ,YQ,v_diff_RF,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,v_diff_RF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12; % Set the desired font size for both axes

clim([-0.08 0.08]);
cmocean('balance','pivot',0)

% title('Overall Velocity Differences, RARE - FSAOO','FontSize',16)
title('RARE - FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'V [m/s]','FontSize',20);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/velocity_diff/RARE_FSAOO/overall.jpg']);

%% GLORYS - FSAOO
figure
hold on
contourf(XQ,YQ,v_diff_GF,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
[C,d] = contour(XQ,YQ,v_diff_GF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12; % Set the desired font size for both axes

clim([-0.08 0.08]);
cmocean('balance','pivot',0)

title('GLORYS - FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'V [m/s]','FontSize',20);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/velocity_diff/GLORYS_FSAOO/overall.jpg']);


%% RARE - GLORYS
figure
hold on
contourf(XQ,YQ,v_diff_RG,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
[C,d] = contour(XQ,YQ,v_diff_RG,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-0.08 0.08]);
cmocean('balance','pivot',0)

title('RARE - GLORYS','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'V [m/s]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/velocity_diff/RARE_GLORYS/overall.jpg']);

%% Plotting Salinity Differences
%[bathx,bathy] = ginput(15);
%% RARE - FSAOO
figure
hold on
contourf(XQ,YQ,SA_diff_RF,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,SA_diff_RF,[-0.1,0,0.1,0.5,1,1.2],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-0.2 1.2]);
cmocean('balance','pivot',0)

title('RARE - FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'SA [g/kg]','FontSize',20);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/salinity_diff/RARE_FSAOO/overall.jpg']);

%% GLORYS - FSAOO
figure
hold on
contourf(XQ,YQ,SA_diff_GF,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,SA_diff_GF,[-0.1,0,0.1,0.5,1,1.2],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-0.2 1.2]);
cmocean('balance','pivot',0)

title('GLORYS - FSAOO','FontSize',22)

xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'SA [g/kg]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/salinity_diff/GLORYS_FSAOO/overall.jpg']);

%% RARE - GLORYS
figure
hold on
contourf(XQ,YQ,SA_diff_RG,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
[C,d] = contour(XQ,YQ,SA_diff_RG,[-0.1,0,0.1,0.5,1,1.2],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-0.2 1.2]);
cmocean('balance','pivot',0)

title('RARE - GLORYS','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'SA [g/kg]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/salinity_diff/RARE_GLORYS/overall.jpg']);

%% Plotting Overall Salinities
figure
hold on
contourf(XQ,YQ,SA_2D_FSAOO,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,SA_2D_FSAOO,[32,33,34,35,35.1],'-k');
clabel(C,d)

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([32 35]);
cmocean('haline')

title('FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'SA [g/kg]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/salinity/FSAOO.jpg']);


figure
hold on
contourf(XQ,YQ,SA_2D_GLORYS,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,SA_2D_GLORYS,[32,33,34,35,35.1],'-k');
clabel(C,d)

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([32 35]);
cmocean('haline')

title('GLORYS','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'SA [g/kg]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/salinity/GLORYS.jpg']);

figure
hold on
contourf(XQ,YQ,SA_2D_RARE,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,SA_2D_RARE,[32,33,34,35,35.1],'-k');
clabel(C,d)

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([32 35]);
cmocean('haline')

title('RARE','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'SA [g/kg]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/salinity/RARE.jpg']);

%% Calculating where there's Atlantic Water

AWt_FSAOO_idx = find(CT_2D_FSAOO>2 & SA_2D_FSAOO>34.92);
AW_FSAOO = zeros(size(CT_2D_FSAOO));
AW_FSAOO(AWt_FSAOO_idx) = 1;

AWt_RARE_idx = find(CT_2D_RARE>2 & SA_2D_RARE>34.92);
AW_RARE = zeros(size(CT_2D_RARE));
AW_RARE(AWt_RARE_idx) = 1;

AWt_GLORYS_idx = find(CT_2D_GLORYS>2 & SA_2D_GLORYS>34.92);
AW_GLORYS = zeros(size(CT_2D_GLORYS));
AW_GLORYS(AWt_GLORYS_idx) = 1;

%% Polar Surface Water

PWt_FSAOO_idx = find(CT_2D_FSAOO<0 & rhoa_2D_FSAOO<27.7);
% PWt_FSAOO_idx = find(CT_2D_FSAOO<0 & SA_2D_FSAOO<34);
PW_FSAOO = zeros(size(CT_2D_FSAOO));
PW_FSAOO(PWt_FSAOO_idx) = 1;

PWt_RARE_idx = find(CT_2D_RARE<0 & rhoa_2D_RARE<27.7);
PW_RARE = zeros(size(CT_2D_RARE));
PW_RARE(PWt_RARE_idx) = 1;

PWt_GLORYS_idx = find(CT_2D_GLORYS<0 & rhoa_2D_GLORYS<27.7);
PW_GLORYS = zeros(size(CT_2D_GLORYS));
PW_GLORYS(PWt_GLORYS_idx) = 1;

%% Plotting Overall Temperature
figure
hold on
contourf(XQ,YQ,CT_2D_FSAOO,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,CT_2D_FSAOO,[-3,-2,-1,0,1,2,3],'-k');
clabel(C,d)

[E,f] = contour(XQ,YQ,AW_FSAOO,[1 1],'-m','LineWidth',2);
[G,h] = contour(XQ,YQ,PW_FSAOO,[1 1],'-w','LineWidth',1.5);


plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-3 3]);
cmocean('thermal')

title('FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'CT [^oC]','FontSize',20);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/temperature/FSAOO.jpg']);

%% RARE

figure
hold on
contourf(XQ,YQ,CT_2D_RARE,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,CT_2D_RARE,[-3,-2,-1,0,1,2,3],'-k');
clabel(C,d)

[E,f] = contour(XQ,YQ,AW_RARE,[1 1],'-m','LineWidth',2);
[G,h] = contour(XQ,YQ,PW_RARE,[1 1],'-w','LineWidth',1.5);


% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-3 3]);
cmocean('thermal')

title('RARE','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'CT [^oC]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/temperature/RARE.jpg']);

%% GLORYS
figure
hold on
contourf(XQ,YQ,CT_2D_GLORYS,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,CT_2D_GLORYS,[-3,-2,-1,0,1,2,3],'-k');
clabel(C,d)

[E,f] = contour(XQ,YQ,AW_GLORYS,[1 1],'-m','LineWidth',2);
[G,h] = contour(XQ,YQ,PW_GLORYS,[1 1],'-w','LineWidth',1.5);


% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-3 3]);
cmocean('thermal')

title('GLORYS','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'CT [^oC]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/temperature/GLORYS.jpg']);

%% Plotting Temperature Differences

% RARE - FSAOO
figure
hold on
contourf(XQ,YQ,CT_diff_RF,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,CT_diff_RF,[-2,-1.5,-1,0,1,2,3],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-2 3]);
cmocean('balance','pivot',0)

title('RARE - FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'CT [^oC]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/temperature_diff/RARE_FSAOO/overall.jpg']);

% GLORYS - FSAOO
figure
hold on
contourf(XQ,YQ,CT_diff_GF,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,CT_diff_GF,[-2,-1.5,-1,0,1,2,3],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-2 3]);
cmocean('balance','pivot',0)

title('GLORYS - FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'CT [^oC]','FontSize',20);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/temperature_diff/GLORYS_FSAOO/overall.jpg']);

%% RARE - GLORYS
figure
hold on
contourf(XQ,YQ,CT_diff_RG,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,CT_diff_RG,[-2,-1.5,-1,0,1,2,3],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-2 3]);
cmocean('balance','pivot',0)

title('RARE - GLORYS','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'CT [^oC]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/temperature_diff/RARE_GLORYS/overall.jpg']);

%% Plotting Density Differences

% RARE - FSAOO
figure
hold on
contourf(XQ,YQ,rhoa_diff_RF,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,rhoa_diff_RF,[-0.25,0,0.25,0.5,0.75,1],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-0.2 1]);
cmocean('balance','pivot',0)

title('RARE - FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'Alpha [kg/m^3]','FontSize',20);
set(gcf,'color','w');

imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/density_diff/RARE_FSAOO/overall.jpg']);

% GLORYS - FSAOO
figure
hold on
contourf(XQ,YQ,rhoa_diff_GF,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,rhoa_diff_GF,[-0.25,0,0.25,0.5,0.75,1],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-0.2 1]);
cmocean('balance','pivot',0)

title('GLORYS - FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'Alpha [kg/m^3]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/density_diff/GLORYS_FSAOO/overall.jpg']);


% RARE - GLORYS
figure
hold on
contourf(XQ,YQ,rhoa_diff_RG,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,rhoa_diff_RF,[-0.25,0,0.25,0.5,0.75,1],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

clim([-0.2 1]);
cmocean('balance','pivot',0)

title('RARE - GLORYS','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'Alpha [kg/m^3]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/density_diff/RARE_FSAOO/RARE_GLORYS.jpg']);

%% Plotting Overall Densities

% RARE
figure
hold on
contourf(XQ,YQ,rhoa_2D_RARE,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,rhoa_2D_RARE,[26:1:36],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

% clim([-0.2 1]);
cmocean('dense')

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

title('RARE','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'Alpha [kg/m^3]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/density/RARE.jpg']);

% GLORYS
figure
hold on
contourf(XQ,YQ,rhoa_2D_GLORYS,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,rhoa_2D_GLORYS,[26:1:36],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

% clim([-0.2 1]);
cmocean('dense')

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

title('GLORYS','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'Alpha [kg/m^3]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/density/GLORYS.jpg']);

% FSAOO
figure
hold on
contourf(XQ,YQ,rhoa_2D_FSAOO,10,'LineStyle','none')
x = colorbar;
fill(bathx,bathy,[0.5 0.5 0.5])
[C,d] = contour(XQ,YQ,rhoa_2D_FSAOO,[26:1:36],'-k');
clabel(C,d)

xlim([-8 -2])
ylim([-1000 0])
ax = gca;
ax.FontSize = 12;

% clim([-0.2 1]);
cmocean('dense')

plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
xline(moorings_x,'--r')

title('FSAOO','FontSize',22)
xlabel('Longitude [ ^oE]','FontSize',20);
ylabel('Depth [m]','FontSize',20);
x.FontSize = 12;
ylabel(x,'Alpha [kg/m^3]','FontSize',20);
set(gcf,'color','w');


imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/density/FSAOO.jpg']);
%% A Quick Foray into T-tests%%
% manually calculating the t statistic and the critical t values so as to
% compare the different data sources

%% Velocity

% RARE - FSAOO
v_1D_RARE = mean(v_interp_RARE,"all",'omitnan');
v_sd_RARE = std(v_interp_RARE,0,'all','omitnan');

v_1D_FSAOO = mean(v_interp_FSAOO,"all",'omitnan');
v_sd_FSAOO = std(v_interp_FSAOO,0,'all','omitnan');

nm = v_1D_RARE-v_1D_FSAOO;
d1 = (v_sd_RARE^2)/(1001*121);
d2 = (v_sd_FSAOO^2)/(1001*121);

t_stat = nm/sqrt(d1+d2);
t_crit = tinv(0.975,(1001*121)-1);

% GLORYS - FSAOO
v_1D_GLORYS = mean(v_interp_GLORYS,"all",'omitnan');
v_sd_GLORYS = std(v_interp_GLORYS,0,'all','omitnan');

v_1D_FSAOO = mean(v_interp_FSAOO,"all",'omitnan');
v_sd_FSAOO = std(v_interp_FSAOO,0,'all','omitnan');

nm = v_1D_GLORYS-v_1D_FSAOO;
d1 = (v_sd_GLORYS^2)/(1001*121);
d2 = (v_sd_FSAOO^2)/(1001*121);

t_stat = nm/sqrt(d1+d2);
t_crit = tinv(0.975,(1001*121)-1);

%% Volume Transport Comparisons %%
lon_1 = -6;
lon_2 = -2;

lat_fs = 78.9;

lon1_idx = find(XQ(1,:)==lon_1);
lon2_idx = find(XQ(1,:)==lon_2);
gridx_diff = lon2_idx-lon1_idx; % how many grid cells are between the two longitudes

lon_km = haversine([lat_fs,lon_1],[lat_fs,lon_2]); % haversign calculates the geometric distance in km 
% between two coordinates
lon_m = lon_km*1000; % converting to m
dx = lon_m/gridx_diff; % horizontal distance of a gridcell

dz = ones(1001,1); % I already specified the depth spacing to be 1m between cells when I interpolated

%% RARE

% January
jan_v_RARE = v_interp_RARE(:,:,1);
jan_egc_RARE = jan_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(jan_egc_RARE>0);
jan_egc_RARE(row,col) = 0;
jan_egc_RARE(isnan(jan_egc_RARE)) = 0;

jan_vx_RARE = sum(dz'*jan_egc_RARE*dx)/10^6;

% February
feb_v_RARE = v_interp_RARE(:,:,2);
feb_egc_RARE = feb_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(feb_egc_RARE>0);
feb_egc_RARE(row,col) = 0;
feb_egc_RARE(isnan(feb_egc_RARE)) = 0;

feb_vx_RARE = sum(dz'*feb_egc_RARE*dx)/10^6;

% March
mar_v_RARE = v_interp_RARE(:,:,3);
mar_egc_RARE = mar_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(mar_egc_RARE>0);
mar_egc_RARE(row,col) = 0;
mar_egc_RARE(isnan(mar_egc_RARE)) = 0;

mar_vx_RARE = sum(dz'*mar_egc_RARE*dx)/10^6;

% April
apr_v_RARE = v_interp_RARE(:,:,4);
apr_egc_RARE = apr_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(apr_egc_RARE>0);
apr_egc_RARE(row,col) = 0;
apr_egc_RARE(isnan(apr_egc_RARE)) = 0;

apr_vx_RARE = sum(dz'*apr_egc_RARE*dx)/10^6;

% May
may_v_RARE = v_interp_RARE(:,:,5);
may_egc_RARE = may_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(may_egc_RARE>0);
may_egc_RARE(row,col) = 0;
may_egc_RARE(isnan(may_egc_RARE)) = 0;

may_vx_RARE = sum(dz'*may_egc_RARE*dx)/10^6;

% June
jun_v_RARE = v_interp_RARE(:,:,6);
jun_egc_RARE = jun_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(jun_egc_RARE>0);
jun_egc_RARE(row,col) = 0;
jun_egc_RARE(isnan(jun_egc_RARE)) = 0;

jun_vx_RARE = sum(dz'*jun_egc_RARE*dx)/10^6;

% July
jul_v_RARE = v_interp_RARE(:,:,7);
jul_egc_RARE = jul_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(jul_egc_RARE>0);
jul_egc_RARE(row,col) = 0;
jul_egc_RARE(isnan(jul_egc_RARE)) = 0;

jul_vx_RARE = sum(dz'*jul_egc_RARE*dx)/10^6;

% August
aug_v_RARE = v_interp_RARE(:,:,8);
aug_egc_RARE = aug_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(aug_egc_RARE>0);
aug_egc_RARE(row,col) = 0;
aug_egc_RARE(isnan(aug_egc_RARE)) = 0;

aug_vx_RARE = sum(dz'*aug_egc_RARE*dx)/10^6;

% September
sep_v_RARE = v_interp_RARE(:,:,9);
sep_egc_RARE = sep_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(sep_egc_RARE>0);
sep_egc_RARE(row,col) = 0;
sep_egc_RARE(isnan(sep_egc_RARE)) = 0;

sep_vx_RARE = sum(dz'*sep_egc_RARE*dx)/10^6;

% October
oct_v_RARE = v_interp_RARE(:,:,10);
oct_egc_RARE = oct_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(oct_egc_RARE>0);
oct_egc_RARE(row,col) = 0;
oct_egc_RARE(isnan(oct_egc_RARE)) = 0;

oct_vx_RARE = sum(dz'*oct_egc_RARE*dx)/10^6;

% November
nov_v_RARE = v_interp_RARE(:,:,11);
nov_egc_RARE = nov_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(nov_egc_RARE>0);
nov_egc_RARE(row,col) = 0;
nov_egc_RARE(isnan(nov_egc_RARE)) = 0;

nov_vx_RARE = sum(dz'*nov_egc_RARE*dx)/10^6;

% December
dec_v_RARE = v_interp_RARE(:,:,12);
dec_egc_RARE = dec_v_RARE(:,lon1_idx:lon2_idx);
[row,col] = find(dec_egc_RARE>0);
dec_egc_RARE(row,col) = 0;
dec_egc_RARE(isnan(dec_egc_RARE)) = 0;

dec_vx_RARE = sum(dz'*dec_egc_RARE*dx)/10^6;

%% GLORYS
% January
jan_v_GLORYS = v_interp_GLORYS(:,:,1);
jan_egc_GLORYS = jan_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(jan_egc_GLORYS>0);
jan_egc_GLORYS(row,col) = 0;
jan_egc_GLORYS(isnan(jan_egc_GLORYS)) = 0;

jan_vx_GLORYS = sum(dz'*jan_egc_GLORYS*dx)/10^6;

% February
feb_v_GLORYS = v_interp_GLORYS(:,:,2);
feb_egc_GLORYS = feb_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(feb_egc_GLORYS>0);
feb_egc_GLORYS(row,col) = 0;
feb_egc_GLORYS(isnan(feb_egc_GLORYS)) = 0;

feb_vx_GLORYS = sum(dz'*feb_egc_GLORYS*dx)/10^6;

% March
mar_v_GLORYS = v_interp_GLORYS(:,:,3);
mar_egc_GLORYS = mar_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(mar_egc_GLORYS>0);
mar_egc_GLORYS(row,col) = 0;
mar_egc_GLORYS(isnan(mar_egc_GLORYS)) = 0;

mar_vx_GLORYS = sum(dz'*mar_egc_GLORYS*dx)/10^6;

% April
apr_v_GLORYS = v_interp_GLORYS(:,:,4);
apr_egc_GLORYS = apr_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(apr_egc_GLORYS>0);
apr_egc_GLORYS(row,col) = 0;
apr_egc_GLORYS(isnan(apr_egc_GLORYS)) = 0;

apr_vx_GLORYS = sum(dz'*apr_egc_GLORYS*dx)/10^6;

% May
may_v_GLORYS = v_interp_GLORYS(:,:,5);
may_egc_GLORYS = may_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(may_egc_GLORYS>0);
may_egc_GLORYS(row,col) = 0;
may_egc_GLORYS(isnan(may_egc_GLORYS)) = 0;

may_vx_GLORYS = sum(dz'*may_egc_GLORYS*dx)/10^6;

% June
jun_v_GLORYS = v_interp_GLORYS(:,:,6);
jun_egc_GLORYS = jun_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(jun_egc_GLORYS>0);
jun_egc_GLORYS(row,col) = 0;
jun_egc_GLORYS(isnan(jun_egc_GLORYS)) = 0;

jun_vx_GLORYS = sum(dz'*jun_egc_GLORYS*dx)/10^6;

% July
jul_v_GLORYS = v_interp_GLORYS(:,:,7);
jul_egc_GLORYS = jul_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(jul_egc_GLORYS>0);
jul_egc_GLORYS(row,col) = 0;
jul_egc_GLORYS(isnan(jul_egc_GLORYS)) = 0;

jul_vx_GLORYS = sum(dz'*jul_egc_GLORYS*dx)/10^6;

% August
aug_v_GLORYS = v_interp_GLORYS(:,:,8);
aug_egc_GLORYS = aug_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(aug_egc_GLORYS>0);
aug_egc_GLORYS(row,col) = 0;
aug_egc_GLORYS(isnan(aug_egc_GLORYS)) = 0;

aug_vx_GLORYS = sum(dz'*aug_egc_GLORYS*dx)/10^6;

% September
sep_v_GLORYS = v_interp_GLORYS(:,:,9);
sep_egc_GLORYS = sep_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(sep_egc_GLORYS>0);
sep_egc_GLORYS(row,col) = 0;
sep_egc_GLORYS(isnan(sep_egc_GLORYS)) = 0;

sep_vx_GLORYS = sum(dz'*sep_egc_GLORYS*dx)/10^6;

% October
oct_v_GLORYS = v_interp_GLORYS(:,:,10);
oct_egc_GLORYS = oct_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(oct_egc_GLORYS>0);
oct_egc_GLORYS(row,col) = 0;
oct_egc_GLORYS(isnan(oct_egc_GLORYS)) = 0;

oct_vx_GLORYS = sum(dz'*oct_egc_GLORYS*dx)/10^6;

% November
nov_v_GLORYS = v_interp_GLORYS(:,:,11);
nov_egc_GLORYS = nov_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(nov_egc_GLORYS>0);
nov_egc_GLORYS(row,col) = 0;
nov_egc_GLORYS(isnan(nov_egc_GLORYS)) = 0;

nov_vx_GLORYS = sum(dz'*nov_egc_GLORYS*dx)/10^6;

% December
dec_v_GLORYS = v_interp_GLORYS(:,:,12);
dec_egc_GLORYS = dec_v_GLORYS(:,lon1_idx:lon2_idx);
[row,col] = find(dec_egc_GLORYS>0);
dec_egc_GLORYS(row,col) = 0;
dec_egc_GLORYS(isnan(dec_egc_GLORYS)) = 0;

dec_vx_GLORYS = sum(dz'*dec_egc_GLORYS*dx)/10^6;

%% FSAOO
% January
jan_v_FSAOO = v_interp_FSAOO(:,:,1);
jan_egc_FSAOO = jan_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(jan_egc_FSAOO>0);
jan_egc_FSAOO(row,col) = 0;
jan_egc_FSAOO(isnan(jan_egc_FSAOO)) = 0;

jan_vx_FSAOO = sum(dz'*jan_egc_FSAOO*dx)/10^6;

% February
feb_v_FSAOO = v_interp_FSAOO(:,:,2);
feb_egc_FSAOO = feb_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(feb_egc_FSAOO>0);
feb_egc_FSAOO(row,col) = 0;
feb_egc_FSAOO(isnan(feb_egc_FSAOO)) = 0;

feb_vx_FSAOO = sum(dz'*feb_egc_FSAOO*dx)/10^6;

% March
mar_v_FSAOO = v_interp_FSAOO(:,:,3);
mar_egc_FSAOO = mar_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(mar_egc_FSAOO>0);
mar_egc_FSAOO(row,col) = 0;
mar_egc_FSAOO(isnan(mar_egc_FSAOO)) = 0;

mar_vx_FSAOO = sum(dz'*mar_egc_FSAOO*dx)/10^6;

% April
apr_v_FSAOO = v_interp_FSAOO(:,:,4);
apr_egc_FSAOO = apr_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(apr_egc_FSAOO>0);
apr_egc_FSAOO(row,col) = 0;
apr_egc_FSAOO(isnan(apr_egc_FSAOO)) = 0;

apr_vx_FSAOO = sum(dz'*apr_egc_FSAOO*dx)/10^6;

% May
may_v_FSAOO = v_interp_FSAOO(:,:,5);
may_egc_FSAOO = may_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(may_egc_FSAOO>0);
may_egc_FSAOO(row,col) = 0;
may_egc_FSAOO(isnan(may_egc_FSAOO)) = 0;

may_vx_FSAOO = sum(dz'*may_egc_FSAOO*dx)/10^6;

% June
jun_v_FSAOO = v_interp_FSAOO(:,:,6);
jun_egc_FSAOO = jun_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(jun_egc_FSAOO>0);
jun_egc_FSAOO(row,col) = 0;
jun_egc_FSAOO(isnan(jun_egc_FSAOO)) = 0;

jun_vx_FSAOO = sum(dz'*jun_egc_FSAOO*dx)/10^6;

% July
jul_v_FSAOO = v_interp_FSAOO(:,:,7);
jul_egc_FSAOO = jul_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(jul_egc_FSAOO>0);
jul_egc_FSAOO(row,col) = 0;
jul_egc_FSAOO(isnan(jul_egc_FSAOO)) = 0;

jul_vx_FSAOO = sum(dz'*jul_egc_FSAOO*dx)/10^6;

% August
aug_v_FSAOO = v_interp_FSAOO(:,:,8);
aug_egc_FSAOO = aug_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(aug_egc_FSAOO>0);
aug_egc_FSAOO(row,col) = 0;
aug_egc_FSAOO(isnan(aug_egc_FSAOO)) = 0;

aug_vx_FSAOO = sum(dz'*aug_egc_FSAOO*dx)/10^6;

% September
sep_v_FSAOO = v_interp_FSAOO(:,:,9);
sep_egc_FSAOO = sep_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(sep_egc_FSAOO>0);
sep_egc_FSAOO(row,col) = 0;
sep_egc_FSAOO(isnan(sep_egc_FSAOO)) = 0;

sep_vx_FSAOO = sum(dz'*sep_egc_FSAOO*dx)/10^6;

% October
oct_v_FSAOO = v_interp_FSAOO(:,:,10);
oct_egc_FSAOO = oct_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(oct_egc_FSAOO>0);
oct_egc_FSAOO(row,col) = 0;
oct_egc_FSAOO(isnan(oct_egc_FSAOO)) = 0;

oct_vx_FSAOO = sum(dz'*oct_egc_FSAOO*dx)/10^6;

% November
nov_v_FSAOO = v_interp_FSAOO(:,:,11);
nov_egc_FSAOO = nov_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(nov_egc_FSAOO>0);
nov_egc_FSAOO(row,col) = 0;
nov_egc_FSAOO(isnan(nov_egc_FSAOO)) = 0;

nov_vx_FSAOO = sum(dz'*nov_egc_FSAOO*dx)/10^6;

% December
dec_v_FSAOO = v_interp_FSAOO(:,:,12);
dec_egc_FSAOO = dec_v_FSAOO(:,lon1_idx:lon2_idx);
[row,col] = find(dec_egc_FSAOO>0);
dec_egc_FSAOO(row,col) = 0;
dec_egc_FSAOO(isnan(dec_egc_FSAOO)) = 0;

dec_vx_FSAOO = sum(dz'*dec_egc_FSAOO*dx)/10^6;

%% Concatenating 

vx_RARE_total = [jan_vx_RARE,feb_vx_RARE,mar_vx_RARE,apr_vx_RARE,may_vx_RARE,...
    jun_vx_RARE,jul_vx_RARE, aug_vx_RARE,sep_vx_RARE,oct_vx_RARE,nov_vx_RARE,...
    dec_vx_RARE];

vx_GLORYS_total = [jan_vx_GLORYS,feb_vx_GLORYS,mar_vx_GLORYS,apr_vx_GLORYS,may_vx_GLORYS,...
    jun_vx_GLORYS,jul_vx_GLORYS, aug_vx_GLORYS,sep_vx_GLORYS,oct_vx_GLORYS,nov_vx_GLORYS,...
    dec_vx_GLORYS];

vx_FSAOO_total = [jan_vx_FSAOO,feb_vx_FSAOO,mar_vx_FSAOO,apr_vx_FSAOO,may_vx_FSAOO,...
    jun_vx_FSAOO,jul_vx_FSAOO, aug_vx_FSAOO,sep_vx_FSAOO,oct_vx_FSAOO,nov_vx_FSAOO,...
    dec_vx_FSAOO];

% calculating correlations
[R_RG,P_RG] = corrcoef(vx_RARE_total,vx_GLORYS_total);
[R_RF,P_RF] = corrcoef(vx_RARE_total,vx_FSAOO_total);
[R_GF,P_GF] = corrcoef(vx_GLORYS_total,vx_FSAOO_total);

%%

% plotting Volume transports
month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

figure
hold on
plot(vx_RARE_total,'-g','LineWidth',2)
plot(vx_GLORYS_total,'-b','LineWidth',2)
plot(vx_FSAOO_total,'-k','LineWidth',2)

xlim([1 12])
ylim([-7 -1])
xticks(1:12);
xticklabels(month)

str = {'R_G_L_O_R_Y_S_,_F_S_A_O_O = 0.6696','R_R_A_R_E_,_F_S_A_O_O = 0.6267',...
    'R_R_A_R_E_,_G_L_O_R_Y_S = 0.6751'};
annotation('textbox', [0.7 0.25 0.2 0.05], 'String', str,'FitBoxToText','on',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontSize', 10, 'BackgroundColor', 'none');

set(gcf,'color','w');
ylabel('Transport [Sv]')
title('EGC Monthly Volume Transport Comparison: RARE, GLORYS, FSAOO 2004-2019')
legend('RARE','GLORYS','FSAOO')
imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
    'Regridded Comparison/vol_tx/overall.jpg']);

%%
test_RARE = v_2D_RARE(:);
test_GLORYS = v_2D_GLORYS(:);
test_FSA00 = v_2D_FSAOO(:);

[h,p,ci,stats] = ttest2(v_2D_GLORYS,v_2D_FSAOO,"Tail","both","Alpha",0.05,"Dim",1);

% both RARE and GLORYS failed the ttest when I reduced the matricies into
% one column vector

% %% Seasonal Comparisons: Velocity
% 
% % Winter (DJF)
% 
% for i = 1:3
%     for j = [1 2 12]
% w_v_diff_RF(:,:,i) = v_interp_RARE(:,:,j)-v_interp_FSAOO(:,:,j);
%     end
% end
% w_v_diff_RF = mean(w_v_diff_RF,3,'omitnan');
% 
% for i = 1:3
%     for j = [1 2 12]
% w_v_diff_GF(:,:,i) = v_interp_GLORYS(:,:,j)-v_interp_FSAOO(:,:,j);
%     end
% end
% w_v_diff_GF = mean(w_v_diff_GF,3,'omitnan');
% 
% %% Plotting Winter Differences
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% subplot(1,2,1)
% hold on
% contourf(XQ,YQ,w_v_diff_RF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,w_v_diff_RF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.08 0.08]);
% cmocean('balance','pivot',0)
% 
% title('RARE - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'V [m/s]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% subplot(1,2,2)
% hold on
% contourf(XQ,YQ,w_v_diff_GF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,w_v_diff_GF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.08 0.08]);
% cmocean('balance','pivot',0)
% 
% title('GLORYS - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'V [m/s]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% sgtitle('Winter Velocity Differences','FontSize',16,'FontWeight','bold') 
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'Regridded Comparison/velocity_diff/seasonal/winter.jpg']);
% 
% %% Spring (MAM)
% 
% for i = 1:3
%     for j = [3 4 5]
% sp_v_diff_RF(:,:,i) = v_interp_RARE(:,:,j)-v_interp_FSAOO(:,:,j);
%     end
% end
% sp_v_diff_RF = mean(sp_v_diff_RF,3,'omitnan');
% 
% for i = 1:3
%     for j = [3 4 5]
% sp_v_diff_GF(:,:,i) = v_interp_GLORYS(:,:,j)-v_interp_FSAOO(:,:,j);
%     end
% end
% sp_v_diff_GF = mean(sp_v_diff_GF,3,'omitnan');
% 
% %% Plotting Spring Differences
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% subplot(1,2,1)
% hold on
% contourf(XQ,YQ,sp_v_diff_RF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,sp_v_diff_RF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.08 0.08]);
% cmocean('balance','pivot',0)
% 
% title('RARE - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'V [m/s]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% subplot(1,2,2)
% hold on
% contourf(XQ,YQ,sp_v_diff_GF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,sp_v_diff_GF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.08 0.08]);
% cmocean('balance','pivot',0)
% 
% title('GLORYS - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'V [m/s]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% sgtitle('Spring Velocity Differences','FontSize',16,'FontWeight','bold')
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'Regridded Comparison/velocity_diff/seasonal/spring.jpg']);
% 
% %% Summer (JJA)
% 
% for i = 1:3
%     for j = [6 7 8]
% su_v_diff_RF(:,:,i) = v_interp_RARE(:,:,j)-v_interp_FSAOO(:,:,j);
%     end
% end
% su_v_diff_RF = mean(su_v_diff_RF,3,'omitnan');
% 
% for i = 1:3
%     for j = [6 7 8]
% su_v_diff_GF(:,:,i) = v_interp_GLORYS(:,:,j)-v_interp_FSAOO(:,:,j);
%     end
% end
% su_v_diff_GF = mean(su_v_diff_GF,3,'omitnan');
% 
% %% Plotting Summer Differences
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% subplot(1,2,1)
% hold on
% contourf(XQ,YQ,su_v_diff_RF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,su_v_diff_RF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.08 0.08]);
% cmocean('balance','pivot',0)
% 
% title('RARE - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'V [m/s]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% subplot(1,2,2)
% hold on
% contourf(XQ,YQ,su_v_diff_GF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,su_v_diff_GF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.08 0.08]);
% cmocean('balance','pivot',0)
% 
% title('GLORYS - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'V [m/s]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% sgtitle('Summer Velocity Differences','FontSize',16,'FontWeight','bold')
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'Regridded Comparison/velocity_diff/seasonal/summer.jpg']);
% 
% %% Autumn (SON)
% 
% for i = 1:3
%     for j = [9 10 11]
% au_v_diff_RF(:,:,i) = v_interp_RARE(:,:,j)-v_interp_FSAOO(:,:,j);
%     end
% end
% au_v_diff_RF = mean(au_v_diff_RF,3,'omitnan');
% 
% for i = 1:3
%     for j = [9 10 11]
% au_v_diff_GF(:,:,i) = v_interp_GLORYS(:,:,j)-v_interp_FSAOO(:,:,j);
%     end
% end
% au_v_diff_GF = mean(au_v_diff_GF,3,'omitnan');
% 
% %% Plotting Summer Differences
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% subplot(1,2,1)
% hold on
% contourf(XQ,YQ,au_v_diff_RF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,au_v_diff_RF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.08 0.08]);
% cmocean('balance','pivot',0)
% 
% title('RARE - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'V [m/s]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% subplot(1,2,2)
% hold on
% contourf(XQ,YQ,au_v_diff_GF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,au_v_diff_GF,[-0.04,-0.02,0,0.02,0.04,0.06,0.8],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.08 0.08]);
% cmocean('balance','pivot',0)
% 
% title('GLORYS - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'V [m/s]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% sgtitle('Autumn Velocity Differences','FontSize',16,'FontWeight','bold')
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'Regridded Comparison/velocity_diff/seasonal/autumn.jpg']);
% 
% %% Seasonal Comparisons: Salinity
% 
% % Winter (DJF)
% 
% for i = 1:3
%     for j = [1 2 12]
% w_SA_diff_RF(:,:,i) = SA_interp_RARE(:,:,j)-SAi(:,:,j);
%     end
% end
% w_SA_diff_RF = mean(w_SA_diff_RF,3,'omitnan');
% 
% for i = 1:3
%     for j = [1 2 12]
% w_SA_diff_GF(:,:,i) = SA_interp_GLORYS(:,:,j)-SAi(:,:,j);
%     end
% end
% w_SA_diff_GF = mean(w_SA_diff_GF,3,'omitnan');
% 
% %% Plotting Winter Differences
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% subplot(1,2,1)
% hold on
% contourf(XQ,YQ,w_SA_diff_RF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,w_SA_diff_RF,[-0.1,0,0.1,0.5,1,1.2],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.2 1.2]);
% cmocean('balance','pivot',0)
% 
% title('RARE - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'SA [g/kg]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% subplot(1,2,2)
% hold on
% contourf(XQ,YQ,w_SA_diff_GF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,w_SA_diff_GF,[-0.1,0,0.1,0.5,1,1.2],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.2 1.2]);
% cmocean('balance','pivot',0)
% 
% title('GLORYS - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'SA [g/kg]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% sgtitle('Winter Salinity Differences','FontSize',16,'FontWeight','bold') 
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'Regridded Comparison/salinity_diff/seasonal/winter.jpg']);
% 
% %% Spring (MAM)
% 
% for i = 1:3
%     for j = [3 4 5]
% sp_SA_diff_RF(:,:,i) = SA_interp_RARE(:,:,j)-SAi(:,:,j);
%     end
% end
% sp_SA_diff_RF = mean(sp_SA_diff_RF,3,'omitnan');
% 
% for i = 1:3
%     for j = [3 4 5]
% sp_SA_diff_GF(:,:,i) = SA_interp_GLORYS(:,:,j)-SAi(:,:,j);
%     end
% end
% sp_SA_diff_GF = mean(sp_SA_diff_GF,3,'omitnan');
% 
% %% Plotting Spring Differences
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% subplot(1,2,1)
% hold on
% contourf(XQ,YQ,sp_SA_diff_RF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,sp_SA_diff_RF,[-0.1,0,0.1,0.5,1,1.2],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.2 1.2]);
% cmocean('balance','pivot',0)
% 
% title('RARE - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'SA [g/kg]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% subplot(1,2,2)
% hold on
% contourf(XQ,YQ,sp_SA_diff_GF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,sp_SA_diff_GF,[-0.1,0,0.1,0.5,1,1.2],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.2 1.2]);
% cmocean('balance','pivot',0)
% 
% title('GLORYS - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'SA [g/kg]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% sgtitle('Spring Salinity Differences','FontSize',16,'FontWeight','bold')
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'Regridded Comparison/salinity_diff/seasonal/spring.jpg']);
% 
% %% Summer (JJA)
% 
% for i = 1:3
%     for j = [6 7 8]
% su_SA_diff_RF(:,:,i) = SA_interp_RARE(:,:,j)-SAi(:,:,j);
%     end
% end
% su_SA_diff_RF = mean(su_SA_diff_RF,3,'omitnan');
% 
% for i = 1:3
%     for j = [6 7 8]
% su_SA_diff_GF(:,:,i) = SA_interp_GLORYS(:,:,j)-SAi(:,:,j);
%     end
% end
% su_SA_diff_GF = mean(su_SA_diff_GF,3,'omitnan');
% 
% %% Plotting Summer Differences
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% subplot(1,2,1)
% hold on
% contourf(XQ,YQ,su_SA_diff_RF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,su_SA_diff_RF,[-0.1,0,0.1,0.5,1,1.2],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.2 1.2]);
% cmocean('balance','pivot',0)
% 
% title('RARE - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'SA [g/kg]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% subplot(1,2,2)
% hold on
% contourf(XQ,YQ,su_SA_diff_GF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,su_SA_diff_GF,[-0.1,0,0.1,0.5,1,1.2],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.2 1.2]);
% cmocean('balance','pivot',0)
% 
% title('GLORYS - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'SA [g/kg]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% sgtitle('Summer Salinity Differences','FontSize',16,'FontWeight','bold')
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'Regridded Comparison/salinity_diff/seasonal/summer.jpg']);
% 
% %% Autumn (SON)
% 
% for i = 1:3
%     for j = [9 10 11]
% au_SA_diff_RF(:,:,i) = SA_interp_RARE(:,:,j)-SAi(:,:,j);
%     end
% end
% au_SA_diff_RF = mean(au_SA_diff_RF,3,'omitnan');
% 
% for i = 1:3
%     for j = [9 10 11]
% au_SA_diff_GF(:,:,i) = SA_interp_GLORYS(:,:,j)-SAi(:,:,j);
%     end
% end
% au_SA_diff_GF = mean(au_SA_diff_GF,3,'omitnan');
% 
% %% Plotting Summer Differences
% 
% figure('color','w','position',[1300 50 1500 900],'papersize',[15 9]);% [21 15]
% subplot(1,2,1)
% hold on
% contourf(XQ,YQ,au_SA_diff_RF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,au_SA_diff_RF,[-0.1,0,0.1,0.5,1,1.2],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.2 1.2]);
% cmocean('balance','pivot',0)
% 
% title('RARE - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'SA [g/kg]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% subplot(1,2,2)
% hold on
% contourf(XQ,YQ,au_SA_diff_GF,10,'LineStyle','none')
% x = colorbar;
% fill(bathx,bathy,[0.5 0.5 0.5],'EdgeColor','none')
% [C,d] = contour(XQ,YQ,au_SA_diff_GF,[-0.1,0,0.1,0.5,1,1.2],'-k');
% clabel(C,d)
% 
% xlim([-8 -2])
% ylim([-1000 0])
% clim([-0.2 1.2]);
% cmocean('balance','pivot',0)
% 
% title('GLORYS - FSAOO','FontSize',20)
% xlabel('Longitude [ ^oE]','FontSize',20);
% ylabel('Depth [m]','FontSize',20);
% ylabel(x,'SA [g/kg]','FontSize',20);
% set(gcf,'color','w');
% hold off
% 
% sgtitle('Autumn Salinity Differences','FontSize',16,'FontWeight','bold')
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'Regridded Comparison/salinity_diff/seasonal/autumn.jpg']);