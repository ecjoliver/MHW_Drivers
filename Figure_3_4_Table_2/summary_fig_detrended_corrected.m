% uses data generated in summaryFig_modesVSmhw_detrended
clear; close all
% files=dir('/Users/alexsengupta/CCRC/My Projects/OceanXtremes/Holbrook- Drivers/summary_fig/processed_MHW_SSTA.pc90.2degree.avhrr-only-v2.*.mat')
files=dir('/Users/alexsengupta/CCRC/My Projects/OceanXtremes/Holbrook- Drivers/corrected_drivers_code/processed_data/detrended_processed_mhw_severity*.mat')

   
for f=1:length(files)
    fn=files(f).name
    load(['/Users/alexsengupta/CCRC/My Projects/OceanXtremes/Holbrook- Drivers/corrected_drivers_code/processed_data/',fn])
    ind=findstr('.',fn);
    ind1=findstr('to',fn);
    mylon=str2num(fn(ind(2)+1:ind1(1)-1));
    mylat=str2num(fn(ind(3)+1:ind1(2)-1));
    
    ii=mylon/30*15+1;
    jj=(mylat+90)/20*10+1;
    matr_modePos_mhwInc_all(:,jj:jj+9,ii:ii+14)=matr_modePos_mhwInc;
    matr_modeNeg_mhwInc_all(:,jj:jj+9,ii:ii+14)=matr_modeNeg_mhwInc;
    matr_modePos_mhwDec_all(:,jj:jj+9,ii:ii+14)=matr_modePos_mhwDec;
    matr_modeNeg_mhwDec_all(:,jj:jj+9,ii:ii+14)=matr_modeNeg_mhwDec;
%     pcolor(squeeze(sum(matr_modeNeg_mhwDec_all)))
%     pause
end

%%%%% error in the saving of medianprop from
% summaryFig_modesVSmhw_detrended_corrected3 (corrected now, but need to
% rerun the code which takes days)
% code below to manually patch together
figure(1);clf;
pcolor(squeeze(sum(matr_modeNeg_mhwInc_all))+squeeze(sum(matr_modePos_mhwInc_all))+squeeze(sum(matr_modeNeg_mhwDec_all))+squeeze(sum(matr_modePos_mhwDec_all)));shading flat
caxis([0 .1])
matr_medianprop_all=zeros(9,90,180);
for f=[30 60 90 108]
    fn=files(f).name
    load(['/Users/alexsengupta/CCRC/My Projects/OceanXtremes/Holbrook- Drivers/corrected_drivers_code/processed_data/',fn])
    figure(2);clf;pcolor(squeeze(medianprop_all(1,:,:)));shading flat
    if f==30; matr_medianprop_all(:,:,1:105)=medianprop_all; end
    if f==60; matr_medianprop_all(:,:,1:150)=matr_medianprop_all(:,:,1:150)+medianprop_all; end
    if f==90; matr_medianprop_all(:,:,1:180)=matr_medianprop_all(:,:,1:180)+medianprop_all; end 
    if f==108; matr_medianprop_all(:,:,1:60)=matr_medianprop_all(:,:,1:60)+medianprop_all; end 
    figure(3);clf;pcolor(squeeze(matr_medianprop_all(1,:,:)));shading flat
    
end
close all
%%%%%%



%%
clear medianprop_all matr_modePos_mhwInc matr_modeNeg_mhwInc matr_modePos_mhwDec matr_modeNeg_mhwDec clear

% COLORMAPS SHOULD BE CENTERED AROUND 8?
% these maps show the regions where a+ve or -ve mode significantly enhances
% or suppresses the no. of MHW days
colormaps
lon=1:2:360;
lat=-89:2:89;
figure(100);clf
for m=1:9
    subplot(3,3,m)
    tmp=squeeze(matr_modePos_mhwInc_all(m,:,:))*100;
    tmp(tmp==0)=NaN;
    pcolor(lon,lat,tmp);shading flat
    title(index_name{m});caxis([3 13]);plotmap
end
colormap(centrewhite)
figure(101);clf
for m=1:9
    subplot(3,3,m)
    tmp=squeeze(matr_modeNeg_mhwInc_all(m,:,:))*100;
    tmp(tmp==0)=NaN;
    pcolor(lon,lat,tmp);shading flat
    title(index_name{m});caxis([3 13]);plotmap
end
colormap(centrewhite)
figure(102);clf
for m=1:9
    subplot(3,3,m)
    tmp=squeeze(matr_modePos_mhwDec_all(m,:,:))*100;
    tmp(tmp==0)=NaN;
    pcolor(lon,lat,tmp);shading flat
    title(index_name{m});caxis([3 13]);plotmap
end
colormap(centrewhite)
figure(103);clf
for m=1:9
    subplot(3,3,m)
    tmp=squeeze(matr_modeNeg_mhwDec_all(m,:,:))*100;
    tmp(tmp==0)=NaN;
    pcolor(lon,lat,tmp);shading flat
    title(index_name{m});caxis([3 13]);plotmap
end
colormap(centrewhite)

figure(104)
colorbar('h')
colormap(centrewhite);caxis([3 13])

%% New plots for revisions
medn=squeeze(matr_medianprop_all(1,:,:))*100; % same for all modes
medn(medn==0)=NaN;
figure(50);clf
pcolor(lon,lat,medn);shading flat ;plotmap;caxis([5 11])
pbaspect([1.8 1 1 ])
colormap(lbmap(21,'BlueRed'));colorbar
lat_lon_labels([],[],[]);
print -f50 -djpeg100 -r400 'images/median_MHW_days'

for m=1:9
inc_dec=squeeze(matr_modePos_mhwInc_all(m,:,:))*100 + squeeze(matr_modePos_mhwDec_all(m,:,:))*100;
inc_dec(inc_dec==0)=NaN;
figure(52);clf
pcolor(lon,lat,(inc_dec-medn)./medn*100);shading flat ;plotmap;caxis([-100 100])
% hold on
% contour(lon,lat,(inc_dec-medn)./medn*100,[-50 50],'k');colorbar
tmp1=(inc_dec-medn)./medn*100;
pbaspect([1.8 1 1 ])
colormap(lbmap(21,'BlueRed'));lat_lon_labels([],[],[]);colorbar;caxis([-80 80]);ylim([-80 80])
title(['Enhanncement and supression of MHW by positive phase of ',index_name{m}])

inc_dec=squeeze(matr_modeNeg_mhwInc_all(m,:,:))*100 + squeeze(matr_modeNeg_mhwDec_all(m,:,:))*100;
inc_dec(inc_dec==0)=NaN;
figure(53);clf
pcolor(lon,lat,(inc_dec-medn)./medn*100);shading flat ;plotmap;caxis([-100 100])
% hold on
% contour(lon,lat,(inc_dec-medn)./medn*100,[-50 50],'k');colorbar
tmp2=(inc_dec-medn)./medn*100;
pbaspect([1.8 1 1 ])
colormap(lbmap(21,'BlueRed'));lat_lon_labels([],[],[]);colorbar;caxis([-80 80]);ylim([-80 80])
title(['Enhanncement and supression of MHW by negative phase of ',index_name{m}])

print('-f52','-djpeg100','-r400',['images/MHW_days_enhance_suppress_POS_',index_name{m}])
print('-f53','-djpeg100','-r400',['images/MHW_days_enhance_suppress_NEG_',index_name{m}])
print('-f52','-depsc',['images/MHW_days_enhance_suppress_POS_',index_name{m}])
print('-f53','-depsc',['images/MHW_days_enhance_suppress_NEG_',index_name{m}])

figure(54);clf
pcolor(lon,lat,tmp1./(-tmp2));shading flat ;plotmap;
caxis([0.8 1.2]);colorbar;colormap(lbmap(21,'BlueRed'))


end

%%

print -f100 -depsc '../images/modePos_mhwInc'
print -f101 -depsc '../images/modeNeg_mhwInc'
print -f102 -depsc '../images/modePos_mhwDec'
print -f103 -depsc '../images/modeNeg_mhwDec'
print -f104 -depsc '../images/mod_mhw_cbar'

%%
mhwInc=cat(1,matr_modePos_mhwInc_all,matr_modeNeg_mhwInc_all);
mhwDec=cat(1,matr_modePos_mhwDec_all,matr_modeNeg_mhwDec_all);
index_name=[index_name index_name];
mhwInc(mhwInc==0)=NaN;
mhwDec(mhwDec==0)=NaN;

for i=1:length(lon)
    for j=1:length(lat)
        tmp=squeeze(mhwInc(:,j,i));
        mhwInc_max_no(j,i)=length(find(tmp==nanmax(tmp)));
        if mhwInc_max_no(j,i)>0
            mhwInc_max(j,i)=find(tmp==nanmax(tmp));
        else
            mhwInc_max(j,i)=NaN;
        end
        
        
        tmp=squeeze(mhwDec(:,j,i));
        mhwDec_min_no(j,i)=length(find(tmp==nanmin(tmp)));
        if mhwDec_min_no(j,i)>0
            mhwDec_min(j,i)=find(tmp==nanmin(tmp));
        else
            mhwDec_min(j,i)=NaN;
        end
    end
end

ind1=find(~isnan(mhwInc_max));
ind2=find(~isnan(mhwDec_min));
unique([mhwInc_max(ind1) ;mhwDec_min(ind2)])

%% custom colormap
cmap=[.9 0.1 0.1 %NAO
    0.0 .6 0.0
    0.2 0.2 .7
    .7 .7  0
    .7 0.1 .7
    0.0 .5 .7
    .3 .3 .3
    .7 .3 .15
    .3 .1 .3
    1 0.5 0.5 %NAO
    0.4 .8 0.5
    0.5 0.5 1
    1 1 0.3
    .75 0.45 .7
    0.3 .8 .8
    .8 .8 .8
    .9 .5 .3
    .5 .3 .5  ];

% cmap=[192 32 38
%     11 143 69
%     53 73 178
%     120 133 54
%     147 38 143
%     0 130 175
%     80 80 80
%     187 87 39
%     92 64 153
%     204 106 114
%     111 194 143
%     120 137 196
%     210 220 120
%     188 114 175
%     90 200 230
%     170 170 170
%     230 135 95
%     140 120 180]/255
    
    
    
 figure(4);clf
pcolor(1:19,1:2,[1:19;1:19]);colormap(cmap)
% print -f4 -depsc 'cmap'   
    
    
    
    
    
    
    
figure(2);clf
pcolor(lon,lat,mhwInc_max*NaN);shading flat
plotmap;colorbar
pbaspect([1.7 1 1]);ylim([-80 80])
% print -f2 -depsc 'blank'

figure(2);clf
pcolor(lon,lat,mhwInc_max);shading flat
plotmap;colorbar
pbaspect([1.7 1 1])
colormap(cmap);ylim([-80 80]);lat_lon_labels([],[],[]);

figure(3);clf
pcolor(lon,lat,mhwDec_min);shading flat
plotmap;colorbar
pbaspect([1.7 1 1])
colormap(cmap);ylim([-80 80]);lat_lon_labels([],[],[]);

print -f2 -depsc 'images/mhwInc_detrended_corrected'
print -f3 -depsc 'images/mhwDec_detrended_corrected'



% nao
% nino34
% pdo
% tpi
% anino
% sam
% modoki
% dmi
% npgo
%%
load('/Users/alexsengupta/CCRC/My Projects/OceanXtremes/CODE/climSST_allyears_20d_smoothed.mat')
figure(5);clf
tmp=squeeze(climSST(1,:,:));
tmp=tmp./tmp;
tmp(isnan(tmp))=0;
contourf(LON(1:4:end),LAT(1:4:end),tmp(1:4:end,1:4:end),[-99 .5 99]);
print -f5 -depsc 'blank'

xlon=[0 LON(1:4:end)' 360];
xlat=[-90 LAT(1:4:end)' 90];
tmp=tmp(1:4:end,1:4:end);
xtmp=ones(length(xlat),length(xlon));
xtmp(2:end-1,2:end-1)=tmp;
clf;contourf(xlon,xlat,xtmp,[-99 .5 99]);xlim([xlon(2) xlon(end-1)]);ylim([xlat(2) xlat(end-1)])
print -f5 -depsc 'blank'
% best to do this landmask in matlab2012

