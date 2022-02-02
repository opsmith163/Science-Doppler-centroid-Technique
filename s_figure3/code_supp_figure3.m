%%%%%%%%%%%% Ekman current retrieval
 clear all; clc
 load('mycolor.mat');
 %%% SAR doppler 10/01  10:04
 filename=strcat('20211001T100446_20211001T100515_039924_04B97E_CDE5.nc');
  lonimg =ncread(filename,'lon');
 latimg =ncread(filename,'lat');
lonimg=(rot90(lonimg,1));
latimg=(rot90(latimg,1));
lon1=mean(lonimg)';
lat1=mean(latimg,2);

 RadVel =ncread(filename,'vv_001_owiRadVel');
EcmwfWindDirection =ncread(filename,'vv_001_owiEcmwfWindDirection');
EcmwfWindSpeed =ncread(filename,'vv_001_owiEcmwfWindSpeed');
WindDirection =ncread(filename,'vv_001_owiWindDirection');
WindSpeed =ncread(filename,'vv_001_owiWindSpeed');
EcmwfWindDirection=(rot90(EcmwfWindDirection,1));
EcmwfWindSpeed=(rot90(EcmwfWindSpeed,1));
WindDirection=(rot90(WindDirection,1));
WindSpeed=(rot90(WindSpeed,1));
RadVel=(rot90(RadVel,1));

 figure('Color',[1 1 1]);
h=imagesc(lon1,lat1,RadVel);set(gca,'YDir','normal');
colormap(rwb2);i=colorbar;
set(h,'alphadata',~isnan(RadVel));
 set(gca,'XTickLabel',{'63W','62.5W','62W','61.5W','61W'}); 
 set(gca,'YTickLabel',{'23.8N','24N','24.2N','24.4N','24.6N','24.8N','25N','25.2N','25.4N','25.6N','25.8N'}) 
% title('Doppler V (m s^{-1})');

pai=3.1416;
weu=EcmwfWindSpeed.*sind(360-EcmwfWindDirection);
wev=EcmwfWindSpeed.*cosd(360-EcmwfWindDirection);
wu=WindSpeed.*sind(360-WindDirection);
wv=WindSpeed.*cosd(360-WindDirection);

 figure('Color',[1 1 1]);
h=imagesc(lon1,lat1,wu);set(gca,'YDir','normal');
colormap(rwb2);i=colorbar;
set(h,'alphadata',~isnan(wu));
 set(gca,'XTickLabel',{'63W','62.5W','62W','61.5W','61W'}); 
 set(gca,'YTickLabel',{'23.8N','24N','24.2N','24.4N','24.6N','24.8N','25N','25.2N','25.4N','25.6N','25.8N'}) 
% title('Wind u (m s^{-1})');

 figure('Color',[1 1 1]);
h=imagesc(lon1,lat1,wv);set(gca,'YDir','normal');
colormap(rwb2);i=colorbar;
set(h,'alphadata',~isnan(wu));
 set(gca,'XTickLabel',{'63W','62.5W','62W','61.5W','61W'}); 
 set(gca,'YTickLabel',{'23.8N','24N','24.2N','24.4N','24.6N','24.8N','25N','25.2N','25.4N','25.6N','25.8N'}) 
% title('Wind v (m s^{-1})');



%%%%%%%% calculate wind stress
pair=1.25;Cd=2.6e-3;
westressx=pair*Cd*EcmwfWindSpeed.*weu;
westressy=pair*Cd*EcmwfWindSpeed.*wev;
wstressx=pair*Cd*WindSpeed.*wu;
wstressy=pair*Cd*WindSpeed.*wv;


%%%%%%%%%%% calculate Coriolis force parameter f
f=coriolisf(lat1);
plot(lat1,f);title('f - lat ');xlabel('lat');ylabel('f');


%%%%%%%%%% Ekman depth retrieval
pwater=1027;
[m,n]=size(RadVel);
for i=1:m
taomag(i,:)=sqrt(2)*sqrt(wstressx(i,:).^2+wstressy(i,:).^2).*pai/(pwater*f(i));
etaomag(i,:)=sqrt(2)*sqrt(westressx(i,:).^2+westressy(i,:).^2).*pai/(pwater*f(i));

end

RadVel=abs(RadVel);
 h1=fspecial('average',[3,3]);
 RadVel=imfilter(RadVel,h1);

DD=taomag./RadVel;
eDD=etaomag./RadVel;
DDshow=DD/(sqrt(2)*pai);
eDDshow=eDD/(sqrt(2)*pai);
eDD(find(eDD>230))=NaN;
DDshow(find(DDshow>100))=NaN;
eDDshow(find(eDDshow>50))=NaN;

 figure('Color',[1 1 1]);
h=imagesc(lon1,lat1,eDDshow);set(gca,'YDir','normal');
colormap(rwb2);i=colorbar;
set(h,'alphadata',~isnan(eDDshow));
 set(gca,'XTickLabel',{'63W','62.5W','62W','61.5W','61W'}); 
 set(gca,'YTickLabel',{'23.8N','24N','24.2N','24.4N','24.6N','24.8N','25N','25.2N','25.4N','25.6N','25.8N'}) 
% title('ECMWF Ekman depth (m)');


%%%%%%%%%% Calculate Ekman current at any depth
load('Hdepth.mat')
level=25;
depth=Hdepth(1:level);
[m,n]=size(RadVel);
for k=1:level
    for i=1:m
    for j=1:n
taozu(i,j,k)=sqrt(2)*exp(-pai*depth(k)/DD(i,j))*(wstressx(i,j)*sin(-pai*depth(k)/DD(i,j)+pai/6)+wstressy(i,j)*cos(-pai*depth(k)/DD(i,j)+pai/6)).*pai/(pwater*f(i)*DD(i,j));
taozv(i,j,k)=sqrt(2)*exp(-pai*depth(k)/DD(i,j))*(wstressy(i,j)*sin(-pai*depth(k)/DD(i,j)+pai/6)-wstressx(i,j)*cos(-pai*depth(k)/DD(i,j)+pai/6)).*pai/(pwater*f(i)*DD(i,j));
etaozu(i,j,k)=sqrt(2)*exp(-pai*depth(k)/eDD(i,j))*(westressx(i,j)*sin(-pai*depth(k)/eDD(i,j)+pai/6)+westressy(i,j)*cos(-pai*depth(k)/eDD(i,j)+pai/6)).*pai/(pwater*f(i)*eDD(i,j));
etaozv(i,j,k)=sqrt(2)*exp(-pai*depth(k)/eDD(i,j))*(westressy(i,j)*sin(-pai*depth(k)/eDD(i,j)+pai/6)-westressx(i,j)*cos(-pai*depth(k)/eDD(i,j)+pai/6)).*pai/(pwater*f(i)*eDD(i,j));
    end
end
end

layer=8;
mag=sqrt(taozu(:,:,layer).^2+taozv(:,:,layer).^2);
[msize,nsize]=size(RadVel);scale=4;

 figure('Color',[1 1 1]);
h=imagesc(lon1(1:scale:nsize),lat1(1:scale:msize),mag(1:scale:msize,1:scale:nsize));
i=colorbar; set(gca,'YDir','normal');
colormap(rwb2);
 set(gca,'XTickLabel',{'63W','62.5W','62W','61.5W','61W'}); 
 set(gca,'YTickLabel',{'23.8N','24N','24.2N','24.4N','24.6N','24.8N','25N','25.2N','25.4N','25.6N','25.8N'}) 
hold on
quiver(lon1(1:scale:nsize),lat1(1:scale:msize),taozu(1:scale:msize,1:scale:nsize,layer),taozv(1:scale:msize,1:scale:nsize,layer)) %���ƶ�άʸ����ͼ
% title('Ekman Sea Surface at 15m (m s^{-1})');



%% Validation

%%% Empirical model from Rio (2014)
% bf0=0.61;sitaf0=-30.75;bf15=0.25;sitaf15=-48.18;
bs0=0.7;sitas0=-32;bs15=0.3;sitas15=-57; % October parameter values

usekm0=bs0*(wstressx*cosd(sitas0)-wstressy*sind(sitas0));
vsekm0=bs0*(wstressx*sind(sitas0)+wstressy*cosd(sitas0));

usekm15=bs15*(wstressx*cosd(sitas15)-wstressy*sind(sitas15));
vsekm15=bs15*(wstressx*sind(sitas15)+wstressy*cosd(sitas15));

%%%%%%%% comparison between retrieval results and empirical model
pusekm0=usekm0(:);
pvsekm0=vsekm0(:);
pusekm15=usekm15(:);
pvsekm15=vsekm15(:);

ru=taozu(:,:,1);
rv=taozv(:,:,1);
pru=ru(:);
prv=rv(:);
ru15=taozu(:,:,8);
rv15=taozv(:,:,8);
pru15=ru15(:);
prv15=rv15(:);

pru(find(isnan(pru)))=0;
prv(find(isnan(prv)))=0;
pru15(find(isnan(pru15)))=0;
prv15(find(isnan(prv15)))=0;
pusekm0(find(isnan(pusekm0)))=0;
pvsekm0(find(isnan(pvsekm0)))=0;
pusekm15(find(isnan(pusekm15)))=0;
pvsekm15(find(isnan(pvsekm15)))=0;

pusekm0(find(pru==0))=[];pvsekm0(find(prv==0))=[];
pusekm15(find(pru15==0))=[];pvsekm15(find(prv15==0))=[];
pru(find(pru==0))=[];prv(find(prv==0))=[];
pru15(find(pru15==0))=[];prv15(find(prv15==0))=[];

samedeg=corrcoef(pusekm0,pru);
simu=samedeg(2,1)
samedeg=corrcoef(pvsekm0,prv);
simv=samedeg(2,1)
samedeg=corrcoef(pusekm15,pru15);
simu15=samedeg(2,1)
samedeg=corrcoef(pvsekm15,prv15);
simv15=samedeg(2,1)


 figure('Color',[1 1 1]);
scatter(pru15,pusekm15,'x','r');
text(0.05,-0.5,num2str(simu15));text(-0.3,-0.5,'u-comp. at 15m corrcoef=');%0446
hold on
scatter(prv15,pvsekm15,'o','k');
xlabel('Retrieved v (m s^{-1})');ylabel('Emprical v (m s^{-1})');%axis([0 2.2 0 2.2]);
text(0.05,-0.54,num2str(simv15));text(-0.3,-0.54,'v-comp. at 15m corrcoef=');%0446
box on
legend({'u component','v component'},2);legend('boxoff');
axis([-0.6 0.3 -0.6 0.3]);
plot([-0.6,0.3],[-0.6,0.3],'k--');
% title('Validation at 15 m ');









