% FS Gridded Data vs RARE (monthly averages)
% Created by: Sara Vianco
% Masters Thesis
% 11 APR 2023

%% Preliminaries
close all
clear all

%%
cd('/Users/saravianco/Documents/Research/MATLAB/theo/');

info1 = ncinfo('Fram_Strait_gridded_monthly_mean_V_S_2003-2019-v1.0.nc');
info2 = ncinfo('Fram_Strait_time-series_2003-2019-v1.0.nc');
v = ncread("Fram_Strait_gridded_monthly_mean_V_S_2003-2019-v1.0.nc",'V_VEL');
    % v = v(:,:,1)';
lon = ncread("Fram_Strait_gridded_monthly_mean_V_S_2003-2019-v1.0.nc",'LONGITUDE')';
lat = ncread("Fram_Strait_gridded_monthly_mean_V_S_2003-2019-v1.0.nc",'LATITUDE');
     lat = mean(lat);
time = ncread("Fram_Strait_gridded_monthly_mean_V_S_2003-2019-v1.0.nc",'TIME');
pres = ncread("Fram_Strait_gridded_monthly_mean_V_S_2003-2019-v1.0.nc",'PRES');
    % pres = pres(:,:,1)';
temp = ncread("Fram_Strait_gridded_monthly_mean_V_S_2003-2019-v1.0.nc",'TEMP');
psal = ncread("Fram_Strait_gridded_monthly_mean_V_S_2003-2019-v1.0.nc",'PSAL');

base = datenum(1950,1,1,0,0,0);
date_num = datenum(time + base);

date_string = datestr(date_num,'mm/dd/yy');
lon_mat = repmat(-lon,56,1);

%% loading colormap cmocean
addpath("/Users/saravianco/Documents/Research/MATLAB/github_repo/");
addpath("/Users/saravianco/Documents/MATLAB/haversine/");

%% Converting pressure to depth (z)
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/html/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/library/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/thermodynamics_from_t/");
addpath("/Users/saravianco/Documents/Research/MATLAB/teos-10/gsw_matlab_v3_06_16/pdf/");

z = (gsw_z_from_p(pres,lat)); 

%% Intial video of the data
% writerObj = VideoWriter('EGS','MPEG-4');
% writerObj.FrameRate = 3;
% open(writerObj);
% 
% fig1 = figure(1);
%  for i = 1:length(v(:,1,1))
% 
% v_t = squeeze(v(i,:,:));
% z_t = squeeze(z(i,:,:));
% 
% surfc(lon_mat,z_t,v_t,'LineStyle','none');
% shading interp
% x = colorbar;
% % contour3(lon_mat,z_t,v_t,[-0.01,0],'-k')
% view([0,0,90])
% caxis([-0.25 0.11]);
% 
% title('East Greenland Shelf Flow 15SEP03-15AUG19','FontSize',16)
% xlabel('Longitude [^oW]','FontSize',14);
% ylabel('Depth [m]');
% ylabel(x,'V [m/s]');
% 
% F = getframe(fig1);
% writeVideo(writerObj, F);
% 
% end
% 
% close(writerObj);

%% A different video approach with contourf
% writerObj = VideoWriter('EGS_2','MPEG-4');
% writerObj.FrameRate = 3;
% open(writerObj);
% 
% fig2 = figure(2);
% hold on
% 
% for i = 1:length(v(:,1,1))
% 
% v_t = squeeze(v(i,:,:));
% z_t = squeeze(z(i,:,:));


% contourf(lon_mat,z_t,v_t);
% shading interp
% x = colorbar;
% hold on
% contour(lon_mat,z_t,v_t,[0,-0.1],'-k');
% clim([-0.25 0.11]);
% 
% title('East Greenland Shelf Flow SEP 2003-AUG 2019','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% F = getframe(fig2);
% writeVideo(writerObj, F);
% 
% end
% close(writerObj);

%% Separating V data into monthly means

formatIn = 'mm/dd/yy';
dv1 = datevec(date_string,formatIn);

% January
jani = find(dv1(:,2) == 1);
jan_v_total = v(jani,:,:);
jan_v_mean = squeeze(mean(jan_v_total,1,'omitnan'));

% February
febi = find(dv1(:,2) == 2);
feb_v_total = v(febi,:,:);
feb_v_mean = squeeze(mean(feb_v_total,1,'omitnan'));

% March
mari = find(dv1(:,2) == 3);
mar_v_total = v(mari,:,:);
mar_v_mean = squeeze(mean(mar_v_total,1,'omitnan'));

% April
apri = find(dv1(:,2) == 4);
apr_v_total = v(apri,:,:);
apr_v_mean = squeeze(mean(apr_v_total,1,'omitnan'));

% May
mayi = find(dv1(:,2) == 5);
may_v_total = v(mayi,:,:);
may_v_mean = squeeze(mean(may_v_total,1,'omitnan'));

% June
juni = find(dv1(:,2) == 6);
jun_v_total = v(juni,:,:);
jun_v_mean = squeeze(mean(jun_v_total,1,'omitnan'));

% July
juli = find(dv1(:,2) == 7);
jul_v_total = v(juli,:,:);
jul_v_mean = squeeze(mean(jul_v_total,1,'omitnan'));

% August
augi = find(dv1(:,2) == 8);
aug_v_total = v(augi,:,:);
aug_v_mean = squeeze(mean(aug_v_total,1,'omitnan'));

% September
sepi = find(dv1(:,2) == 9);
sep_v_total = v(sepi,:,:);
sep_v_mean = squeeze(mean(sep_v_total,1,'omitnan'));

% October
octi = find(dv1(:,2) == 10);
oct_v_total = v(octi,:,:);
oct_v_mean = squeeze(mean(oct_v_total,1,'omitnan'));

% November
novi = find(dv1(:,2) == 11);
nov_v_total = v(novi,:,:);
nov_v_mean = squeeze(mean(nov_v_total,1,'omitnan'));

% December
deci = find(dv1(:,2) == 12);
dec_v_total = v(deci,:,:);
dec_v_mean = squeeze(mean(dec_v_total,1,'omitnan'));

%% Individual Monthly plots
cd('/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_figures/');
 z_t = squeeze(z(1,:,:));
% 
% % January
% figure(1)
% hold on
% contourf(lon_mat,z_t,jan_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,h] = contour(lon_mat,z_t,jan_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: January','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(1),'01_theo.jpg');
% 
% % February
% figure(2)
% contourf(lon_mat,z_t,feb_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,feb_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: February','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(2),'02_theo.jpg');
% 
% % March
% figure(3)
% contourf(lon_mat,z_t,mar_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,mar_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: March','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(3),'03_theo.jpg');
% 
% % April
% figure(4)
% contourf(lon_mat,z_t,apr_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,apr_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: April','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(4),'04_theo.jpg');
% 
% % May
% figure(5)
% contourf(lon_mat,z_t,may_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,may_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: May','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(5),'05_theo.jpg');
% 
% % June
% figure(6)
% contourf(lon_mat,z_t,jun_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,jun_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: June','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(6),'06_theo.jpg');
% 
% % July
% figure(7)
% contourf(lon_mat,z_t,jul_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,jul_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: July','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(7),'07_theo.jpg');
% 
% % August
% figure(8)
% contourf(lon_mat,z_t,aug_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,aug_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: August','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(8),'08_theo.jpg');
% 
% % September
% figure(9)
% contourf(lon_mat,z_t,sep_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,sep_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: September','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(9),'09_theo.jpg');
% 
% % October
% figure(10)
% contourf(lon_mat,z_t,oct_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,oct_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% 
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: October','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(10),'10_theo.jpg');
% 
% % November
% figure(11)
% contourf(lon_mat,z_t,nov_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,nov_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% 
% title('EGS Mean Flow: November','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(11),'11_theo.jpg');
% 
% % December
% figure(12)
% contourf(lon_mat,z_t,dec_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,dec_v_mean,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% title('EGS Mean Flow: December','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(12),'12_theo.jpg');
% 
%% Total mean
% 
v_mean_total = squeeze(mean(v(:,:,:)));
% 
% figure(13)
% contourf(lon_mat,z_t,v_mean_total,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% [C,h] = contour(lon_mat,z_t,v_mean_total,[0.1,0,-0.05,-0.1,-0.15,-0.2],'-k');
% set(gca, 'XDir','reverse')
% clabel(C,h);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0) 
% 
% title('EGS Mean Flow (FSAOO)','FontSize',16)
% xlabel('Longitude [ ^oW]','FontSize',14);
% ylabel('Depth [m]','FontSize',14);
% ylabel(x,'V [m/s]','FontSize',14);
% 
% saveas(figure(13),'total_mean.jpg');

%% Monthly V Anomalies
% cd('/Users/saravianco/Documents/Research/MATLAB/theo/');
% 
% % January
% jan_v_anom = jan_v_mean-v_mean_total;
% 
% % February
% feb_v_anom = feb_v_mean-v_mean_total;
% 
% % March
% mar_v_anom = mar_v_mean-v_mean_total;
% 
% % April
% apr_v_anom = apr_v_mean-v_mean_total;
% 
% % May
% may_v_anom = may_v_mean-v_mean_total;
% 
% % June
% jun_v_anom = jun_v_mean-v_mean_total;
% 
% % July
% jul_v_anom = jul_v_mean-v_mean_total;
% 
% % August
% aug_v_anom = aug_v_mean-v_mean_total;
% 
% % September
% sep_v_anom = sep_v_mean-v_mean_total;
% 
% % October
% oct_v_anom = oct_v_mean-v_mean_total;
% 
% % November
% nov_v_anom = nov_v_mean-v_mean_total;
% 
% % December
% dec_v_anom = dec_v_mean-v_mean_total;

%% Imaging V Anomalies

% % January
% figure
% contourf(lon_mat,z_t,jan_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,jan_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous January EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/01a_theo.jpg');
% 
% % February
% figure
% contourf(lon_mat,z_t,feb_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,feb_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous February EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/02a_theo.jpg');
% 
% % March
% figure
% contourf(lon_mat,z_t,mar_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,mar_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous March EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/03a_theo.jpg');
% 
% % April
% figure
% contourf(lon_mat,z_t,apr_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,apr_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous April EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/04a_theo.jpg');
% 
% % May
% figure
% contourf(lon_mat,z_t,may_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,may_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous May EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/05a_theo.jpg');
% 
% % June
% figure
% contourf(lon_mat,z_t,jun_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,jun_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous June EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/06a_theo.jpg');
% 
% % July
% figure
% contourf(lon_mat,z_t,jul_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,jul_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous July EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/07a_theo.jpg');
% 
% % August
% figure
% contourf(lon_mat,z_t,aug_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,aug_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous August EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/08a_theo.jpg');
% 
% % September
% figure
% contourf(lon_mat,z_t,sep_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,sep_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous September EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/09a_theo.jpg');
% 
% % October
% figure
% contourf(lon_mat,z_t,oct_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,oct_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous October EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/10a_theo.jpg');
% 
% % November
% figure
% contourf(lon_mat,z_t,nov_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,nov_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous November EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/11a_theo.jpg');
% 
% % December
% figure
% contourf(lon_mat,z_t,dec_v_anom,8,'linestyle','none');
% shading interp
% x = colorbar;
% hold on
% set(gca, 'XDir','reverse')
% clim([-0.04 0.04]);
% [C,h] = contour(lon_mat,z_t,dec_v_anom,[0.04,0.02,0.01,0,-0.01,-0.02,-0.04],'-k');
% clabel(C,h);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oW]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,"V' [m/s]",'FontSize',14)
% title('Anomolous December EGS Flow','FontSize',16);
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_v_anom_figs/12a_theo.jpg');

%% Volumetric Transport

    addpath('/Volumes/TERABYTE/RARE/matching_theo');

    lon_1 = -6;
    lon_2 = -2;

    lon_1_idx = find(lon_mat(1,:)==-lon_1);
    lon_2_idx = find(lon_mat(1,:)==-lon_2);

    gridx_diff = lon_2_idx-lon_1_idx;

    lon_km = haversine([lat,lon_1],[lat,lon_2]);
    lon_m = lon_km*1000;

    dx = lon_m/gridx_diff;% horizontal distance of a gridcell

    z_pos = -z_t;

    idx_1000 = NaN(1,25);
    for i = 1:25
    idx_1000(1,i) = find(z_pos(:,i)<=1000,1,'last');

    end

    dz = NaN(size(z_pos));
    for i = 1:length(z_pos(1,:))
        for kk = 2:(length(z_pos(:,1))-1)
            dz(kk,i) = ((z_pos(kk+1,i)+z_pos(kk,i))/2)-((z_pos(kk,i)+z_pos(kk-1,i))/2);
            dz(1,i) = ((z_pos(2,i)+z_pos(1,i))/2)-(z_pos(1,i)/2);
        end
    end

for j = 1:25
    dz(56,j) = (z_pos(56,j)-z_pos(55,j))/2;

end

    max_idx = max(idx_1000);

    dz_egc = zeros(56,17);
    for i = 1:17
    dz_egc(1:idx_1000(i+8),i) = dz(1:idx_1000(i+8),i+8);
    end

    %% Calculating transport

    % January
    jan_egc = jan_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(jan_egc>0);
    jan_egc(row,col) = 0;
    jan_tx = (sum((dx*jan_egc.*dz_egc),'all'))/10^6;

    % February
    feb_egc = feb_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(feb_egc>0);
    feb_egc(row,col) = 0;
    feb_tx = (sum((dx*feb_egc.*dz_egc),'all'))/10^6;

    % March
    mar_egc = mar_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(mar_egc>0);
    mar_egc(row,col) = 0;
    mar_tx = (sum((dx*mar_egc.*dz_egc),'all'))/10^6;

    % April
    apr_egc = apr_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(apr_egc>0);
    apr_egc(row,col) = 0;
    apr_tx = (sum((dx*apr_egc.*dz_egc),'all'))/10^6;

    % May
    may_egc = may_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(may_egc>0);
    may_egc(row,col) = 0;
    may_tx = (sum((dx*may_egc.*dz_egc),'all'))/10^6;

    % June
    jun_egc = jun_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(jun_egc>0);
    jun_egc(row,col) = 0;
    jun_tx = (sum((dx*jun_egc.*dz_egc),'all'))/10^6;

    % July
    jul_egc = jul_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(jul_egc>0);
    jul_egc(row,col) = 0;
    jul_tx = (sum((dx*jul_egc.*dz_egc),'all'))/10^6;

    % August
    aug_egc = aug_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(aug_egc>0);
    aug_egc(row,col) = 0;
    aug_tx = (sum((dx*aug_egc.*dz_egc),'all'))/10^6;

    % September
    sep_egc = sep_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(sep_egc>0);
    sep_egc(row,col) = 0;
    sep_tx = (sum((dx*sep_egc.*dz_egc),'all'))/10^6;

    % October
    oct_egc = oct_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(oct_egc>0);
    oct_egc(row,col) = 0;
    oct_tx = (sum((dx*oct_egc.*dz_egc),'all'))/10^6;

    % November
    nov_egc = nov_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(nov_egc>0);
    nov_egc(row,col) = 0;
    nov_tx = (sum((dx*nov_egc.*dz_egc),'all'))/10^6;

    % December
    dec_egc = dec_v_mean(:,lon_1_idx:lon_2_idx);
    [row,col] = find(dec_egc>0);
    dec_egc(row,col) = 0;
    dec_tx = (sum((dx*dec_egc.*dz_egc),'all'))/10^6;

    %%
    monthly_tx = [jan_tx,feb_tx,mar_tx,apr_tx,may_tx,jun_tx,jul_tx,...
        aug_tx,sep_tx,oct_tx,nov_tx,dec_tx];
    month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

figure
plot(monthly_tx,'-k','linewi',2)
set(gca,'Ydir','reverse')
xlim([1 12])
xticks(1:12);
xticklabels(month)

ylabel('Transport [Sv]')
title('EGC Average Monthly Volumetric Transport (FSAOO)')
saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_tx/vol_tx_FSAOO.jpg');

%%
% 
% z_1000 = NaN(1,17);
% 
% for i = 1:17
% z_1000(1,i) = z_pos(idx_1000(i+8),i+8);
% end
% 
% diff_RARE = 1064.3-z_1000;
% diff_area_RARE = dx*diff_RARE;
% diff_atot_RARE = sum((dx*diff_RARE),'all');
% diff_tx_RARE = (-0.025*diff_atot_RARE)/10^6;
% 
% monthly_tx_adj = monthly_tx-diff_tx_RARE;
% 
% figure
% plot(monthly_tx_adj,'-k','linewi',2)
% set(gca,'Ydir','reverse')
% %ylim([-5.3 -2.7])
% xlim([1 12])
% xticks(1:12);
% xticklabels(month)
% 
% ylabel('Transport [Sv]')
% title('EGC Average Monthly Volumetric Transport (FSAOO)')
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_tx/vol_tx_FSAOO_adj.jpg');


%% Full Section

% lon_1 = lon(1);
% lon_2 = lon(end);
% 
% lon_1_idx = find(lon_mat(1,:)==-lon_1);
% lon_2_idx = find(lon_mat(1,:)==-lon_2);
% 
% gridx_diff = lon_2_idx-lon_1_idx;
% 
% lon_km = haversine([lat,lon_1],[lat,lon_2]);
% lon_m = lon_km*1000;
% 
% dx = lon_m/gridx_diff;% horizontal distance of a gridcell
% 
% z_pos = -z_t;
% 
% dz = NaN(size(z_pos));
% for i = 1:length(z_pos(1,:))
%     for j = 1:(length(z_pos(:,1))-1)
%         dz(j,i) = z_pos(j+1,i)-z_pos(j,i);
%     end
% end
% dz(isnan(dz)) = 0;
% 
% %% Calculating transport
% 
% % January
% jan_tx = (sum((dx*jan_v_mean.*dz),'all'))/10^6;
% 
% % February
% feb_tx = (sum((dx*feb_v_mean.*dz),'all'))/10^6;
% 
% % March
% mar_tx = (sum((dx*mar_v_mean.*dz),'all'))/10^6;
% 
% % April
% apr_tx = (sum((dx*apr_v_mean.*dz),'all'))/10^6;
% 
% % May
% may_tx = (sum((dx*may_v_mean.*dz),'all'))/10^6;
% 
% % June
% jun_tx = (sum((dx*jun_v_mean.*dz),'all'))/10^6;
% 
% % July
% jul_tx = (sum((dx*jul_v_mean.*dz),'all'))/10^6;
% 
% % August
% aug_tx = (sum((dx*aug_v_mean.*dz),'all'))/10^6;
% 
% % September
% sep_tx = (sum((dx*sep_v_mean.*dz),'all'))/10^6;
% 
% % October
% oct_tx = (sum((dx*oct_v_mean.*dz),'all'))/10^6;
% 
% % November
% nov_tx = (sum((dx*nov_v_mean.*dz),'all'))/10^6;
% 
% % December
% dec_tx = (sum((dx*dec_v_mean.*dz),'all'))/10^6;
% 
% % imaging
% monthly_tx = [jan_tx,feb_tx,mar_tx,apr_tx,may_tx,jun_tx,jul_tx,...
%     aug_tx,sep_tx,oct_tx,nov_tx,dec_tx];
% month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% 
% figure
% plot(monthly_tx,'-k','linewi',2)
% set(gca,'Ydir','reverse')
% ylim([-10 -3])
% xlim([1 12])
% xticks(1:12);
% xticklabels(month)
% 
% ylabel('Transport [Sv]')
% title('Full Section Average Monthly Volumetric Transport (FSAOO)')
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_tx/vol_tx_FSAOO_fullsct.jpg');

%% Density Sections

%% Separating Temperature variable into seperate months

% % January
% jani = find(dv1(:,2) == 1);
% jan_t_total = temp(jani,:,:);
% jan_t_mean = squeeze(mean(jan_t_total,1,'omitnan'));
% 
% % February
% febi = find(dv1(:,2) == 2);
% feb_t_total = temp(febi,:,:);
% feb_t_mean = squeeze(mean(feb_t_total,1,'omitnan'));
% 
% % March
% mari = find(dv1(:,2) == 3);
% mar_t_total = temp(mari,:,:);
% mar_t_mean = squeeze(mean(mar_t_total,1,'omitnan'));
% 
% % April
% apri = find(dv1(:,2) == 4);
% apr_t_total = temp(apri,:,:);
% apr_t_mean = squeeze(mean(apr_t_total,1,'omitnan'));
% 
% % May
% mayi = find(dv1(:,2) == 5);
% may_t_total = temp(mayi,:,:);
% may_t_mean = squeeze(mean(may_t_total,1,'omitnan'));
% 
% % June
% juni = find(dv1(:,2) == 6);
% jun_t_total = temp(juni,:,:);
% jun_t_mean = squeeze(mean(jun_t_total,1,'omitnan'));
% 
% % July
% juli = find(dv1(:,2) == 7);
% jul_t_total = temp(juli,:,:);
% jul_t_mean = squeeze(mean(jul_t_total,1,'omitnan'));
% 
% % August
% augi = find(dv1(:,2) == 8);
% aug_t_total = temp(augi,:,:);
% aug_t_mean = squeeze(mean(aug_t_total,1,'omitnan'));
% 
% % September
% sepi = find(dv1(:,2) == 9);
% sep_t_total = temp(sepi,:,:);
% sep_t_mean = squeeze(mean(sep_t_total,1,'omitnan'));
% 
% % October
% octi = find(dv1(:,2) == 10);
% oct_t_total = temp(octi,:,:);
% oct_t_mean = squeeze(mean(oct_t_total,1,'omitnan'));
% 
% % Notember
% novi = find(dv1(:,2) == 11);
% nov_t_total = temp(novi,:,:);
% nov_t_mean = squeeze(mean(nov_t_total,1,'omitnan'));
% 
% % December
% deci = find(dv1(:,2) == 12);
% dec_t_total = temp(deci,:,:);
% dec_t_mean = squeeze(mean(dec_t_total,1,'omitnan'));
% 
% %% Separating Prac. Salinity variable into seperate months
% 
% % January
% jani = find(dv1(:,2) == 1);
% jan_s_total = psal(jani,:,:);
% jan_s_mean = squeeze(mean(jan_s_total,1,'omitnan'));
% 
% % February
% febi = find(dv1(:,2) == 2);
% feb_s_total = psal(febi,:,:);
% feb_s_mean = squeeze(mean(feb_s_total,1,'omitnan'));
% 
% % March
% mari = find(dv1(:,2) == 3);
% mar_s_total = psal(mari,:,:);
% mar_s_mean = squeeze(mean(mar_s_total,1,'omitnan'));
% 
% % April
% apri = find(dv1(:,2) == 4);
% apr_s_total = psal(apri,:,:);
% apr_s_mean = squeeze(mean(apr_s_total,1,'omitnan'));
% 
% % May
% mayi = find(dv1(:,2) == 5);
% may_s_total = psal(mayi,:,:);
% may_s_mean = squeeze(mean(may_s_total,1,'omitnan'));
% 
% % June
% juni = find(dv1(:,2) == 6);
% jun_s_total = psal(juni,:,:);
% jun_s_mean = squeeze(mean(jun_s_total,1,'omitnan'));
% 
% % July
% juli = find(dv1(:,2) == 7);
% jul_s_total = psal(juli,:,:);
% jul_s_mean = squeeze(mean(jul_s_total,1,'omitnan'));
% 
% % August
% augi = find(dv1(:,2) == 8);
% aug_s_total = psal(augi,:,:);
% aug_s_mean = squeeze(mean(aug_s_total,1,'omitnan'));
% 
% % September
% sepi = find(dv1(:,2) == 9);
% sep_s_total = psal(sepi,:,:);
% sep_s_mean = squeeze(mean(sep_s_total,1,'omitnan'));
% 
% % October
% octi = find(dv1(:,2) == 10);
% oct_s_total = psal(octi,:,:);
% oct_s_mean = squeeze(mean(oct_s_total,1,'omitnan'));
% 
% % Notember
% novi = find(dv1(:,2) == 11);
% nov_s_total = psal(novi,:,:);
% nov_s_mean = squeeze(mean(nov_s_total,1,'omitnan'));
% 
% % December
% deci = find(dv1(:,2) == 12);
% dec_s_total = psal(deci,:,:);
% dec_s_mean = squeeze(mean(dec_s_total,1,'omitnan'));
% 
% 
% %% Concatenating into 3d arrays
% s_mean = cat(3,jan_s_mean, feb_s_mean, mar_s_mean, apr_s_mean, may_s_mean, jun_s_mean,... 
%     jul_s_mean, aug_s_mean, sep_s_mean, oct_s_mean, nov_s_mean, dec_s_mean);
% 
% t_mean = cat(3,jan_t_mean, feb_t_mean, mar_t_mean, apr_t_mean, may_t_mean, jun_t_mean,... 
%     jul_t_mean, aug_t_mean, sep_t_mean, oct_t_mean, nov_t_mean, dec_t_mean);
% 
%% Individual variables for density calculation
% lon_mat = repmat(lon,56,1);
% pres_set = squeeze(pres(1,:,:));
% 
% SA = NaN(size(s_mean));
% for i = 1:12
% [SA(:,:,i), in_ocean] = gsw_SA_from_SP(s_mean(:,:,i),pres_set,lon_mat,lat);
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
% rho_CT_exact(:,:,i) = gsw_rho_CT_exact(SA(:,:,i),CT(:,:,i),pres_set);
% end
% 
% rho_anom = rho_CT_exact-1000;

%% Adding Mooring locations
moorings_x = [-8 -6.5 -5.5 -5 -4 -3 -2];
moorings_y = zeros(1,7);

%% Imaging Density Contours

% % January
% figure
% hold on
% contourf(lon_mat,z_t,jan_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,1),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(jan_v_mean));
% imAlpha(isnan(jan_v_mean))=0;
% imagesc(jan_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean January EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/monthly_density/01_rho.jpg']);
% 
% % February
% figure
% hold on
% contourf(lon_mat,z_t,feb_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,2),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(feb_v_mean));
% imAlpha(isnan(feb_v_mean))=0;
% imagesc(feb_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean February EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/monthly_density/02_rho.jpg']);
% 
% % March
% figure
% hold on
% contourf(lon_mat,z_t,mar_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,3),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(mar_v_mean));
% imAlpha(isnan(mar_v_mean))=0;
% imagesc(mar_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean March EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/monthly_density/03_rho.jpg']);
% 
% % April
% figure
% hold on
% contourf(lon_mat,z_t,apr_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,4),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(apr_v_mean));
% imAlpha(isnan(apr_v_mean))=0;
% imagesc(apr_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean April EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/monthly_density/04_rho.jpg']);
% 
% % May
% figure
% hold on
% contourf(lon_mat,z_t,may_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,5),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(may_v_mean));
% imAlpha(isnan(may_v_mean))=0;
% imagesc(may_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean May EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/monthly_density/05_rho.jpg']);
% 
% % June
% figure
% hold on
% contourf(lon_mat,z_t,jun_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,6),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(jun_v_mean));
% imAlpha(isnan(jun_v_mean))=0;
% imagesc(jun_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean June EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%      'theo/monthly_density/06_rho.jpg']);
% 
% % July
% figure
% hold on
% contourf(lon_mat,z_t,jul_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,7),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(jul_v_mean));
% imAlpha(isnan(jul_v_mean))=0;
% imagesc(jul_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean July EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%      'theo/monthly_density/07_rho.jpg']);
% 
% % August
% figure
% hold on
% contourf(lon_mat,z_t,aug_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,8),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(aug_v_mean));
% imAlpha(isnan(aug_v_mean))=0;
% imagesc(aug_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean August EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%      'theo/monthly_density/08_rho.jpg']);
% 
% % September
% figure
% hold on
% contourf(lon_mat,z_t,sep_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,9),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(sep_v_mean));
% imAlpha(isnan(sep_v_mean))=0;
% imagesc(sep_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean September EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%      'theo/monthly_density/09_rho.jpg']);
% 
% % October
% figure
% hold on
% contourf(lon_mat,z_t,oct_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,10),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(oct_v_mean));
% imAlpha(isnan(oct_v_mean))=0;
% imagesc(oct_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean October EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%      'theo/monthly_density/10_rho.jpg']);
% 
% % November
% figure
% hold on
% contourf(lon_mat,z_t,nov_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,11),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(nov_v_mean));
% imAlpha(isnan(nov_v_mean))=0;
% imagesc(nov_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean November EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%      'theo/monthly_density/11_rho.jpg']);
% 
% % December
% figure
% hold on
% contourf(lon_mat,z_t,dec_v_mean,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(lon_mat,z_t,rho_anom(:,:,12),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% 
% imAlpha=ones(size(dec_v_mean));
% imAlpha(isnan(dec_v_mean))=0;
% imagesc(dec_v_mean,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean December EGS Velocity and Density (FSAOO)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%      'theo/monthly_density/12_rho.jpg']);

%% Volumetric Transport with Error bars

%January
jan_egc_mat = NaN(16,56,17);
jan_tx_mat = NaN(16,1);

for ii = 1:16
jan_egc_mat(ii,:,:) = jan_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(jan_egc_mat(:,:,:)>0);
jan_egc_mat(pos_idx) = 0;
jan_tx_mat(ii,1) = (sum((dx*(squeeze(jan_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% February
feb_egc_mat = NaN(16,56,17);
feb_tx_mat = NaN(16,1);

for ii = 1:16
feb_egc_mat(ii,:,:) = feb_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(feb_egc_mat(:,:,:)>0);
feb_egc_mat(pos_idx) = 0;
feb_tx_mat(ii,1) = (sum((dx*(squeeze(feb_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% March
mar_egc_mat = NaN(16,56,17);
mar_tx_mat = NaN(16,1);

for ii = 1:16
mar_egc_mat(ii,:,:) = mar_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(mar_egc_mat(:,:,:)>0);
mar_egc_mat(pos_idx) = 0;
mar_tx_mat(ii,1) = (sum((dx*(squeeze(mar_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% April
apr_egc_mat = NaN(16,56,17);
apr_tx_mat = NaN(16,1);

for ii = 1:16
apr_egc_mat(ii,:,:) = apr_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(apr_egc_mat(:,:,:)>0);
apr_egc_mat(pos_idx) = 0;
apr_tx_mat(ii,1) = (sum((dx*(squeeze(apr_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% May
may_egc_mat = NaN(16,56,17);
may_tx_mat = NaN(16,1);

for ii = 1:16
may_egc_mat(ii,:,:) = may_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(may_egc_mat(:,:,:)>0);
may_egc_mat(pos_idx) = 0;
may_tx_mat(ii,1) = (sum((dx*(squeeze(may_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% June
jun_egc_mat = NaN(16,56,17);
jun_tx_mat = NaN(16,1);

for ii = 1:16
jun_egc_mat(ii,:,:) = jun_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(jun_egc_mat(:,:,:)>0);
jun_egc_mat(pos_idx) = 0;
jun_tx_mat(ii,1) = (sum((dx*(squeeze(jun_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% July
jul_egc_mat = NaN(16,56,17);
jul_tx_mat = NaN(16,1);

for ii = 1:16
jul_egc_mat(ii,:,:) = jul_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(jul_egc_mat(:,:,:)>0);
jul_egc_mat(pos_idx) = 0;
jul_tx_mat(ii,1) = (sum((dx*(squeeze(jul_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% August
aug_egc_mat = NaN(16,56,17);
aug_tx_mat = NaN(16,1);

for ii = 1:16
aug_egc_mat(ii,:,:) = aug_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(aug_egc_mat(:,:,:)>0);
aug_egc_mat(pos_idx) = 0;
aug_tx_mat(ii,1) = (sum((dx*(squeeze(aug_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% September
sep_egc_mat = NaN(16,56,17);
sep_tx_mat = NaN(16,1);

for ii = 1:16
sep_egc_mat(ii,:,:) = sep_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(sep_egc_mat(:,:,:)>0);
sep_egc_mat(pos_idx) = 0;
sep_tx_mat(ii,1) = (sum((dx*(squeeze(sep_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% October
oct_egc_mat = NaN(16,56,17);
oct_tx_mat = NaN(16,1);

for ii = 1:16
oct_egc_mat(ii,:,:) = oct_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(oct_egc_mat(:,:,:)>0);
oct_egc_mat(pos_idx) = 0;
oct_tx_mat(ii,1) = (sum((dx*(squeeze(oct_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% November
nov_egc_mat = NaN(16,56,17);
nov_tx_mat = NaN(16,1);

for ii = 1:16
nov_egc_mat(ii,:,:) = nov_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(nov_egc_mat(:,:,:)>0);
nov_egc_mat(pos_idx) = 0;
nov_tx_mat(ii,1) = (sum((dx*(squeeze(nov_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end

% December
dec_egc_mat = NaN(16,56,17);
dec_tx_mat = NaN(16,1);

for ii = 1:16
dec_egc_mat(ii,:,:) = dec_v_total(ii,:,lon_1_idx:lon_2_idx);
pos_idx = find(dec_egc_mat(:,:,:)>0);
dec_egc_mat(pos_idx) = 0;
dec_tx_mat(ii,1) = (sum((dx*(squeeze(dec_egc_mat(ii,:,:))).*dz_egc),'all'))/10^6;
end


%% reconfiguring the full tx so that it matches

full_tx_fsaoo = NaN(17,12);
full_tx_fsaoo(2:end,1) = jan_tx_mat;
full_tx_fsaoo(2:end,2) = feb_tx_mat;
full_tx_fsaoo(2:end,3) = mar_tx_mat;
full_tx_fsaoo(2:end,4) = apr_tx_mat;
full_tx_fsaoo(2:end,5) = may_tx_mat;
full_tx_fsaoo(2:end,6) = jun_tx_mat;
full_tx_fsaoo(2:end,7) = jul_tx_mat;
full_tx_fsaoo(2:end,8) = aug_tx_mat;

full_tx_fsaoo(1:end-1,9) = sep_tx_mat;
full_tx_fsaoo(1:end-1,10) = oct_tx_mat;
full_tx_fsaoo(1:end-1,11) = nov_tx_mat;
full_tx_fsaoo(1:end-1,12) = dec_tx_mat;

%% Determining Errorbars/SD/95% Confidence bars
sd_vx_fsaoo = NaN(1,12);
mean_vx_fsaoo = NaN(1,12);

for jj = 1:12
sd_vx_fsaoo(1,jj) = std(full_tx_fsaoo(:,jj),'omitnan');
mean_vx_fsaoo(1,jj) = mean(full_tx_fsaoo(:,jj),'omitnan');
end

 %% Figures
% 
% % Volumetric Transport for each year, at each month. The black line
% % represents the average tx for that month
% c = cmocean('thermal',16);
% 
% figure
% hold on
% for ii = 1:16
% plot(full_tx(ii,:),'linewi',1,'Color',c(ii,:))
% end
% plot(mean_vx_fsaoo,'-k','linewi',3)
% xlim([1 12])
% ylim([-10 1])
% xticks(1:12);
% xticklabels(month)
% 
% legend('2004','2005','2006','2007','2008','2009','2010','2011','2012' ...
%     ,'2013','2014','2015','2016','2017','2018','2019','location','east outside')
% 
% ylabel('Transport [Sv]')
% title('EGC Monthly Volumetric Transport 2003-2019 (FSAOO)')
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_tx/FSAOO_EGC_vx_all.jpg');
% 
% %% Plotting the mean across all year for each month as well as error bars (1
% % SD)
% 
% figure
% hold on
% plot(mean_vx_fsaoo,'-k','linewi',2)
% shadedErrorBar([1:12],mean_vx_fsaoo,sd_vx_fsaoo,'lineprops','b')
% for ii = 1:16
% c = cmocean('thermal',16);
% scatter([1:12],full_tx(ii,:),50,c(ii,:),'o','filled')
% end
% xlim([1 12])
% ylim([-8 -1])
% xticks(1:12);
% xticklabels(month)
% legend('','','2004','2005','2006','2007','2008','2009','2010','2011','2012' ...
%     ,'2013','2014','2015','2016','2017','2018','2019','location','east outside')
% ylabel('Transport [Sv]')
% title('EGC Average Monthly Volumetric Transport 2003-2019 (FSAOO)')
% saveas(gcf,'/Users/saravianco/Documents/Research/MATLAB/theo/monthly_tx/FSAOO_EGC_vx_err.jpg');

%% Re-gridding the FSAOO Mooring Data

v_fsaoo = permute(v,[2,3,1]);
temp_fsaoo = permute(temp,[2,3,1]);
psal_fsaoo = permute(psal,[2,3,1]);
pres_fsaoo = permute(pres,[2,3,1]);

% With 1D interp

% interpolating for depth
xq = -(0:1:1000)';
% vq_test = interp1(z_t(:,1),v_test(:,1),xq);

v_interp = NaN(1001,25,192);
for i = 1:192
    for j = 1:25
v_interp(:,j,i) = interp1(z_t(:,j),v_fsaoo(:,j,i),xq);
    end
end

% interpolating for longitudinal resolution

yq = (-8:0.05:-2);

v_interp_fsaoo = NaN(1001,121,192);
for i = 1:192
    for j = 1:1001
v_interp_fsaoo(j,:,i) = interp1(lon,v_interp(j,:,i),yq);       
    end
end

%% Regridding temp/salinity

% temperature
t_interp = NaN(1001,25,192);
for i = 1:192
    for j = 1:25
t_interp(:,j,i) = interp1(z_t(:,j),temp_fsaoo(:,j,i),xq);
    end
end

% interpolating for longitudinal resolution

t_interp_fsaoo = NaN(1001,121,192);
for i = 1:192
    for j = 1:1001
t_interp_fsaoo(j,:,i) = interp1(lon,t_interp(j,:,i),yq);       
    end
end

% salinity
s_interp = NaN(1001,25,192);
for i = 1:192
    for j = 1:25
s_interp(:,j,i) = interp1(z_t(:,j),psal_fsaoo(:,j,i),xq);
    end
end

% interpolating for longitudinal resolution

s_interp_fsaoo = NaN(1001,121,192);
for i = 1:192
    for j = 1:1001
s_interp_fsaoo(j,:,i) = interp1(lon,s_interp(j,:,i),yq);       
    end
end

% pressure
p_interp = NaN(1001,25,192);
for i = 1:192
    for j = 1:25
p_interp(:,j,i) = interp1(z_t(:,j),pres_fsaoo(:,j,i),xq);
    end
end

% interpolating for longitudinal resolution

p_interp_fsaoo = NaN(1001,121,192);
for i = 1:192
    for j = 1:1001
p_interp_fsaoo(j,:,i) = interp1(lon,p_interp(j,:,i),yq);       
    end
end

% usable pressure for density calculations
pres_interp_2D = p_interp_fsaoo(:,:,1);

%% Separating into monthly variables
% * the i in the variable names represent 'interpolated'

% January
jan_v_mati = v_interp_fsaoo(:,:,jani);
jan_t_mati = t_interp_fsaoo(:,:,jani);
jan_s_mati = s_interp_fsaoo(:,:,jani);

jan_v_meani = mean(jan_v_mati,3,'omitnan');
jan_t_meani = mean(jan_t_mati,3,'omitnan');
jan_s_meani = mean(jan_s_mati,3,'omitnan');

% February
feb_v_mati = v_interp_fsaoo(:,:,febi);
feb_t_mati = t_interp_fsaoo(:,:,febi);
feb_s_mati = s_interp_fsaoo(:,:,febi);

feb_v_meani = mean(feb_v_mati,3,'omitnan');
feb_t_meani = mean(feb_t_mati,3,'omitnan');
feb_s_meani = mean(feb_s_mati,3,'omitnan');

% March
mar_v_mati = v_interp_fsaoo(:,:,mari);
mar_t_mati = t_interp_fsaoo(:,:,mari);
mar_s_mati = s_interp_fsaoo(:,:,mari);

mar_v_meani = mean(mar_v_mati,3,'omitnan');
mar_t_meani = mean(mar_t_mati,3,'omitnan');
mar_s_meani = mean(mar_s_mati,3,'omitnan');

% April
apr_v_mati = v_interp_fsaoo(:,:,apri);
apr_t_mati = t_interp_fsaoo(:,:,apri);
apr_s_mati = s_interp_fsaoo(:,:,apri);

apr_v_meani = mean(apr_v_mati,3,'omitnan');
apr_t_meani = mean(apr_t_mati,3,'omitnan');
apr_s_meani = mean(apr_s_mati,3,'omitnan');

% May
may_v_mati = v_interp_fsaoo(:,:,mayi);
may_t_mati = t_interp_fsaoo(:,:,mayi);
may_s_mati = s_interp_fsaoo(:,:,mayi);

may_v_meani = mean(may_v_mati,3,'omitnan');
may_t_meani = mean(may_t_mati,3,'omitnan');
may_s_meani = mean(may_s_mati,3,'omitnan');

% June
jun_v_mati = v_interp_fsaoo(:,:,juni);
jun_t_mati = t_interp_fsaoo(:,:,juni);
jun_s_mati = s_interp_fsaoo(:,:,juni);

jun_v_meani = mean(jun_v_mati,3,'omitnan');
jun_t_meani = mean(jun_t_mati,3,'omitnan');
jun_s_meani = mean(jun_s_mati,3,'omitnan');

% July
jul_v_mati = v_interp_fsaoo(:,:,juli);
jul_t_mati = t_interp_fsaoo(:,:,juli);
jul_s_mati = s_interp_fsaoo(:,:,juli);

jul_v_meani = mean(jul_v_mati,3,'omitnan');
jul_t_meani = mean(jul_t_mati,3,'omitnan');
jul_s_meani = mean(jul_s_mati,3,'omitnan');

% August
aug_v_mati = v_interp_fsaoo(:,:,augi);
aug_t_mati = t_interp_fsaoo(:,:,augi);
aug_s_mati = s_interp_fsaoo(:,:,augi);

aug_v_meani = mean(aug_v_mati,3,'omitnan');
aug_t_meani = mean(aug_t_mati,3,'omitnan');
aug_s_meani = mean(aug_s_mati,3,'omitnan');

% September
sep_v_mati = v_interp_fsaoo(:,:,sepi);
sep_t_mati = t_interp_fsaoo(:,:,sepi);
sep_s_mati = s_interp_fsaoo(:,:,sepi);

sep_v_meani = mean(sep_v_mati,3,'omitnan');
sep_t_meani = mean(sep_t_mati,3,'omitnan');
sep_s_meani = mean(sep_s_mati,3,'omitnan');

% October
oct_v_mati = v_interp_fsaoo(:,:,octi);
oct_t_mati = t_interp_fsaoo(:,:,octi);
oct_s_mati = s_interp_fsaoo(:,:,octi);

oct_v_meani = mean(oct_v_mati,3,'omitnan');
oct_t_meani = mean(oct_t_mati,3,'omitnan');
oct_s_meani = mean(oct_s_mati,3,'omitnan');

% November
nov_v_mati = v_interp_fsaoo(:,:,novi);
nov_t_mati = t_interp_fsaoo(:,:,novi);
nov_s_mati = s_interp_fsaoo(:,:,novi);

nov_v_meani = mean(nov_v_mati,3,'omitnan');
nov_t_meani = mean(nov_t_mati,3,'omitnan');
nov_s_meani = mean(nov_s_mati,3,'omitnan');

% December
dec_v_mati = v_interp_fsaoo(:,:,deci);
dec_t_mati = t_interp_fsaoo(:,:,deci);
dec_s_mati = s_interp_fsaoo(:,:,deci);

dec_v_meani = mean(dec_v_mati,3,'omitnan');
dec_t_meani = mean(dec_t_mati,3,'omitnan');
dec_s_meani = mean(dec_s_mati,3,'omitnan');

%% Density calculations
[LON,Z] = meshgrid(yq,xq);

% concatenating temp and sal into 3D arrays
s_meani = cat(3,jan_s_meani, feb_s_meani, mar_s_meani, apr_s_meani, may_s_meani, jun_s_meani,... 
    jul_s_meani, aug_s_meani, sep_s_meani, oct_s_meani, nov_s_meani, dec_s_meani);

t_meani = cat(3,jan_t_meani, feb_t_meani, mar_t_meani, apr_t_meani, may_t_meani, jun_t_meani,... 
    jul_t_meani, aug_t_meani, sep_t_meani, oct_t_meani, nov_t_meani, dec_t_meani);


SAi = NaN(size(s_meani));
for i = 1:12
[SAi(:,:,i), in_ocean] = gsw_SA_from_SP(s_meani(:,:,i),pres_interp_2D,LON,lat);
end

CTi = NaN(size(t_meani));
for i = 1:12
CTi(:,:,i) = gsw_CT_from_pt(SAi(:,:,i),t_meani(:,:,i));
end

rho_i = NaN(size(t_meani));
for i = 1:12
rho_i(:,:,i) = gsw_rho_CT_exact(SAi(:,:,i),CTi(:,:,i),pres_interp_2D);
end

rho_anom_i = rho_i-1000;

%% Concatenating Velocity Data

v_interp_FSAOO = cat(3,jan_v_meani,feb_v_meani,mar_v_meani,apr_v_meani,may_v_meani,...
    jun_v_meani,jul_v_meani,aug_v_meani,sep_v_meani,oct_v_meani,nov_v_meani,dec_v_meani);

%% Plotting Velocity and Density 

% % January
% figure
% hold on
% contourf(LON,Z,jan_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,1),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(jan_v_meani));
% imAlpha(isnan(jan_v_meani))=0;
% imagesc(jan_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean January EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/01.jpg']);
% 
% % February
% figure
% hold on
% contourf(LON,Z,feb_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,2),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(feb_v_meani));
% imAlpha(isnan(feb_v_meani))=0;
% imagesc(feb_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean February EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/02.jpg']);
% 
% % March
% figure
% hold on
% contourf(LON,Z,mar_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,3),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(mar_v_meani));
% imAlpha(isnan(mar_v_meani))=0;
% imagesc(mar_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean March EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/03.jpg']);
% 
% % April
% figure
% hold on
% contourf(LON,Z,apr_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,4),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(apr_v_meani));
% imAlpha(isnan(apr_v_meani))=0;
% imagesc(apr_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean April EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/04.jpg']);
% 
% % May
% figure
% hold on
% contourf(LON,Z,may_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,5),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(may_v_meani));
% imAlpha(isnan(may_v_meani))=0;
% imagesc(may_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean May EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/05.jpg']);
% 
% % June
% figure
% hold on
% contourf(LON,Z,jun_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,6),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(jun_v_meani));
% imAlpha(isnan(jun_v_meani))=0;
% imagesc(jun_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean June EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/06.jpg']);
% 
% % July
% figure
% hold on
% contourf(LON,Z,jul_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,7),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(jul_v_meani));
% imAlpha(isnan(jul_v_meani))=0;
% imagesc(jul_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean July EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/07.jpg']);
% 
% % August
% figure
% hold on
% contourf(LON,Z,aug_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,8),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(aug_v_meani));
% imAlpha(isnan(aug_v_meani))=0;
% imagesc(aug_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean August EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/08.jpg']);
% 
% % September
% figure
% hold on
% contourf(LON,Z,sep_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,9),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(sep_v_meani));
% imAlpha(isnan(sep_v_meani))=0;
% imagesc(sep_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean September EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/09.jpg']);
% 
% % October
% figure
% hold on
% contourf(LON,Z,oct_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,10),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(oct_v_meani));
% imAlpha(isnan(oct_v_meani))=0;
% imagesc(oct_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean October EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/10.jpg']);
% 
% % November
% figure
% hold on
% contourf(LON,Z,nov_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,11),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(nov_v_meani));
% imAlpha(isnan(nov_v_meani))=0;
% imagesc(nov_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean November EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/11.jpg']);
% 
% % December
% figure
% hold on
% contourf(LON,Z,dec_v_meani,8,'linestyle','none');
% shading interp
% x = colorbar;
% [C,d] = contour(LON,Z,rho_anom_i(:,:,12),[25:0.5:32.5],'-k');
% clabel(C,d);
% 
% plot(moorings_x,moorings_y,'.r','Marker','^','MarkerSize',8,'MarkerFaceColor','r');
% xline(moorings_x,'--r')
% % filling in the bathymetry with a grey color
% imAlpha=ones(size(dec_v_meani));
% imAlpha(isnan(dec_v_meani))=0;
% imagesc(dec_v_meani,'AlphaData',imAlpha);
% set(gca,'color',[0.5 0.5 0.5]);
% 
% xlim([-8 -2]);
% ylim([-1000 0]);
% clim([-0.25 0.11]);
% cmocean('balance','pivot',0)
% 
% xlabel('Longitude [^oE]','FontSize',14)
% ylabel('Depth [m]','FontSize',14)
% ylabel(x,'V [m/s]','FontSize',14)
% title('Mean December EGS Velocity and Density (FSAOO-Regridded)','FontSize',16);
% set(gcf,'color','w');
% 
% imwrite(getframe(gcf).cdata, ['/Users/saravianco/Documents/Research/MATLAB/' ...
%     'theo/regridded plots/v_rho/12.jpg'])