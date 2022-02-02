clear all;clc
 load('mycolor.mat');

%%%%%%%%%% Doppler veloity figure
% only for drawing figure
load('lat1.mat');
load('lon1.mat');
load('DopplerVelpres.mat') 
figure('Color',[1 1 1]);
h=imagesc(lon1,lat1,RadVel);set(gca,'YDir','normal');
colormap(rwb2);
% colorbar('YTickLabel',{'0','0.5','1','1.5','2'})
i=colorbar;% For presenting in pic
 set(h,'alphadata',~isnan(RadVel));
 set(gca,'XTickLabel',{'62.5W','62W','61.5W','61W','60.5W','60W','59.5W'}) 
 set(gca,'YTickLabel',{'38.4N','38.6N','38.8N','39N','39.2N','39.4N','39.6N','39.8N','40N','40.2N'}) 

 %%%%%%%% Doppler veloity 
% For calculating
load('lat1.mat');
load('lon1.mat');
load('DopplerVelscatt.mat') % 
h=imagesc(lon1,lat1,RadVel);set(gca,'YDir','normal');
% colormap(hsv);
colormap(rwb2);
colorbar('YTickLabel',{'0','0.5','1','1.5','2'})
% i=colorbar;set(get(i,'Title'),'string','m/s');

%%%%%%%%%%%% HYCOM magnitude of surface vector field
% 07/08 21:00 2021 Hycom 
 load('Hycom_LAT1.mat');
 load('Hycom_LON1.mat');
 load('Hycom_MAG1.mat');
figure('Color',[1 1 1]);
h=imagesc(LON1,LAT1,MAG1);set(gca,'YDir','normal');
% colormap(hsv);
colormap(rwb2);
 i=colorbar;%set(get(i,'Title'),'string','m/s');
set(h,'alphadata',~isnan(MAG1));
 set(gca,'XTickLabel',{'62.5W','62W','61.5W','61W','60.5W','60W','59.5W'}) 
 set(gca,'YTickLabel',{'38.4N','38.6N','38.8N','39N','39.2N','39.4N','39.6N','39.8N','40N','40.2N'}) 

 % 07/09 00:00 2021 Hycom
load('Hycom_LAT2.mat');
load('Hycom_LON2.mat');
load('Hycom_MAG2.mat');

figure('Color',[1 1 1]);
h=imagesc(LON2,LAT2,MAG2);set(gca,'YDir','normal');
% colormap(hsv);
colormap(rwb2);
 i=colorbar;%set(get(i,'Title'),'string','m/s');
set(h,'alphadata',~isnan(MAG2));
 set(gca,'XTickLabel',{'62.5W','62W','61.5W','61W','60.5W','60W','59.5W'}) 
 set(gca,'YTickLabel',{'38.4N','38.6N','38.8N','39N','39.2N','39.4N','39.6N','39.8N','40N','40.2N'}) 

 %%%%%%%%%%% Plot scatter for Doppler velocity and HYCOM magnitude
 % 07/08 21:00 2021 Hycom vs. Doppler velocity
  pmag1=(MAG1(:));
 pV1=abs(RadVel(:));
 pmag1(find(pV1==0))=[];
  pV1(find(pV1==0))=[];

samedegmag=corrcoef(pmag1,pV1);
 simmag1=samedegmag(2,1) 
 sz=5;
 figure('Color',[1 1 1]);
  c = linspace(1,100,length(pmag1));
 scatter(pmag1,pV1,sz,c);ylabel('Doppler V (m s^{-1})');
text(1.8,0.3,num2str(simmag1));text(1.5,0.3,'corrcoef=');
xlabel('Hycom mag (m s^{-1})');axis([0 2.2 0 2.2]);
hold on
plot([0,2.2],[0,2.2],'k--');
box on;
title('Hycom mag - Doppler Vel ');

% 07/09 00:00 2021 Hycom vs. Doppler velocity
 pmag2=(MAG2(:));
 pV2=abs(RadVel(:));
 pmag2(find(pV2==0))=[];
 pV2(find(pV2==0))=[];
samedegmag=corrcoef(pmag2,pV2);
 simmag2=samedegmag(2,1)
 sz=5;
 
 figure('Color',[1 1 1]);
  c = linspace(1,100,length(pmag2));
  scatter(pmag2,pV2,sz,c);%scatter(pmag2,pV2,sz,c);
  ylabel('Doppler V (m s^{-1})');
text(1.8,0.3,num2str(simmag2));text(1.5,0.3,'corrcoef=');
xlabel('Hycom mag (m s^{-1})');axis([0 2.2 0 2.2]);
hold on
plot([0,2.2],[0,2.2],'k--');
box on;
title('Hycom mag - Doppler Vel ');
 

