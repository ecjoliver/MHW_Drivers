function varargout = MHWControlGUI_updated(varargin)
warning off
% move through MHW datset
% Can apply minimum area filter (Value=X represents X*X degrees minimum area)
% Can apply minimum intensity filter such that anomaly must exceed Value in oC

% NB to see SSTA, need to have minimum area filter checked

% MHWCONTROLGUI MATLAB code for MHWControlGUI.fig
%      MHWCONTROLGUI, by itself, creates a new MHWCONTROLGUI or raises the existing
%      singleton*.
%
%      H = MHWCONTROLGUI returns the handle to a new MHWCONTROLGUI or the handle to
%      the existing singleton*.
%
%      MHWCONTROLGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MHWCONTROLGUI.M with the given input arguments.
%
%      MHWCONTROLGUI('Property','Value',...) creates a new MHWCONTROLGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MHWControlGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MHWControlGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MHWControlGUI

% Last Modified by GUIDE v2.5 03-Aug-2017 13:16:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MHWControlGUI_updated_OpeningFcn, ...
    'gui_OutputFcn',  @MHWControlGUI_updated_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function MHWControlGUI_updated_plotdata(handles)

mydateindex=handles.myyear*10000 + handles.mymonth*100 + handles.myday;
indd=find(mydateindex==handles.dateindex);


set(handles.MYyear,'string',handles.myyear)
set(handles.MYmonth,'string',handles.mymonth)
set(handles.MYday,'string',handles.myday)
lon=handles.lon;
lat=handles.lat;

if isempty(indd)
    disp('No corresponding date available')
else
    
    if handles.pc90check.Value % if 90th pc on
        figure(1);clf
        MHW=zeros(length(handles.lat),length(handles.lon));
        handles.E=handles.E-30;
        if handles.W<=handles.E
            for i=handles.W:30:handles.E
                for j=[handles.S:20:handles.N]
                    fn=['/Users/z3045790/OneDrive/DATASETS/mhw_12_2018/mhw_severity.pc90.',num2str(i),'to',num2str(i+30),'.',num2str(j),'to',num2str(j+20),'.1981.2018.nc'];
                    mhw=squeeze(nc_varget(fn,'severity',[indd-1 0 0],[1 -1 -1]));
                    jj=(j+90)*4+1;
                    ii=i*4+1;
                    MHW(jj:jj+79,ii:ii+119)=mhw;
                end
            end
        else
            % cross 0oE
            startlons=[handles.W:30:359 0:30:handles.E];
            startlons(startlons>180)=startlons(startlons>180)-360;
            startlons=startlons+180;
            cc=0;
            for i=[handles.W:30:359 0:30:handles.E]
                cc=cc+1;
                for j=[handles.S:20:handles.N]
                    fn=['/Users/z3045790/OneDrive/DATASETS/mhw_12_2018/mhw_severity.pc90.',num2str(i),'to',num2str(i+30),'.',num2str(j),'to',num2str(j+20),'.1981.2018.nc'];
                    mhw=squeeze(nc_varget(fn,'severity',[indd-1 0 0],[1 -1 -1]));
                    jj=(j+90)*4+1;
                    ii=startlons(cc)*4+1;
                    MHW(jj:jj+79,ii:ii+119)=mhw;
                end
            end
            lon2=lon(721:end)-360;lon2=[lon2 ;lon(1:720)];
            lon=lon2;handles.lon=lon; clear lon2
        end
        MHW(MHW==0)=NaN;
        
        mask_area=MHW;mask_area(:)=1;
        mask_intensity=mask_area;
        
        if handles.AreaFilter.Value
            area=grid_area(lat,lon);
            % select minumum area for MHW
            min_area=(str2num(handles.areamin.String)*(2*pi*earthRadius/360))^2;
            
            mhw=MHW;
            mask=mhw; ind=find(mhw==0);mask(:)=1;mask(ind)=NaN;
            mhw(mhw==0)=NaN;
            mhw(isnan(mhw))=0;
            % if both area and intesity filter - apply the intesity filter
            % first
            if handles.IntensityFilter.Value
                ind=find(mhw<str2num(handles.intensitymin.String));
                mhw(ind)=0;
            end
            
            [B,L] = bwboundaries(mhw,'noholes');
            
            % if greater than minimum size add to new B & L and store area and cummulative intensity
            cc=1;
            clear B2 mhw_area cum_intense
            L2=L*0;
            for n=1:length(B)
                ind=find(L==n);
                tmp=sum(area(ind));
                if tmp>min_area
                    L2(ind)=cc;
                    B2{cc}=B{n};
                    mhw_area(cc)=tmp;
                    cum_intense(cc)=sum(mhw(ind).*area(ind));
                    cc=cc+1;
                end
            end
            L2(L2==0)=NaN;
            if exist('mhw_area','var')
                mhw_area=mhw_area/1000/1000; % km2
                cum_intense=cum_intense/1000/1000; %oC.km2
            end
            mask_area=L2./L2;
            
        end
        if handles.IntensityFilter.Value & ~handles.AreaFilter.Value
            ind=find(MHW>str2num(handles.intensitymin.String));
            mask_intensity(:)=NaN;
            mask_intensity(ind)=1;
        end
        
        pcolor(lon,lat,mask_area.*mask_intensity.*MHW);shading flat
        caxis([0 2])
        plotmap
        title(['y=',num2str(handles.myyear), ' m=',num2str(handles.mymonth), ' d=',num2str( handles.myday)])
        colorbar
        xW= handles.W; if handles.E<handles.W; xW=xW-360; end
        xlim([xW handles.E+30]);ylim([handles.S handles.N+20])
        
        if ~isempty(handles.X) & handles.AreaFilter.Value
            figure(1);plot(handles.X,handles.Y,'r+','MarkerSize',20,'LineWidth',3)
            %loop through heatwaves to find closest
            clear RNG
            if exist('mhw_area','var')
                for n=1:length(mhw_area)
                    [indj,indi]=find(L2==n);
                    [RNG(n), AZ] = distance(handles.Y,handles.X,mean(lat(indj)),mean(lon(indi)));
                end
                ind=find(RNG==min(RNG));% select closest
                clear boundary90
                if RNG(ind)<str2num(handles.min_distance.String) % make sure MHW is within certain radius
                    selectedMHWindex=find(L2==ind);
                    boundary90=B2{ind};
                    plot(lon(boundary90(:,2)),lat(boundary90(:,1)),'r') % plot boundary around closest MHW
                    tmparea=num2str(sum(area(selectedMHWindex))/1000/1000/1e6);
                    tmpint=num2str(max(MHW(selectedMHWindex)));
                    tmp90=num2str(prctile(MHW(selectedMHWindex),90));
                    set(handles.stat_area,'string',tmparea ); %km2
                    set(handles.stat_intensity,'string', tmpint); %oC
                    set(handles.pc90,'string', tmp90); %oC
                    %                 disp(['area: ',tmparea,' intensity: ',tmpint,' intensity90pc: ',tmp90,])
                else
                    set(handles.stat_area,'string', ''); %km2
                    set(handles.stat_intensity,'string', ''); %oC
                    set(handles.pc90,'string', ''); %oC
                end
            else
                set(handles.stat_area,'string', ''); %km2
                set(handles.stat_intensity,'string', ''); %oC
                set(handles.pc90,'string', ''); %oC
            end
        end
    end
    
    %%
    
    if handles.pc98check.Value % if 98th pc checked
        figure(2);clf
        MHW=zeros(length(handles.lat),length(handles.lon));
        if handles.W<=handles.E
            for i=handles.W:30:handles.E
                for j=[handles.S:20:handles.N]
%                     fn=['/Users/alexsengupta/Desktop/MHWpc98/MHW_SSTA.pc98.sst.avhrr-only-v2.1981-2014.i',num2str(i),'.j',num2str(j),'.nc'];
%                     fn=['/Volumes/HDD_MacRetina/DATASETS2/NOAA OI SST V2 High Resolution Dataset/REGIONAL_timeseries/heatwaves/MHW_SSTA.pc98.sst.avhrr-only-v2.1981-2014.i',num2str(i),'.j',num2str(j),'.nc'];
%                     mhw=squeeze(nc_varget(fn,'mhw',[indd-1 0 0],[1 -1 -1]));
%                     jj=j*4+1;
%                     ii=i*4+1;
%                     MHW(jj:jj+119,ii:ii+239)=mhw;
                    
                    fn=['/Users/z3045790/OneDrive/DATASETS/mhw_98pc/mhw_severity.pc98.',num2str(i),'to',num2str(i+30),'.',num2str(j),'to',num2str(j+20),'.1981.2019.nc'];
                    mhw=squeeze(nc_varget(fn,'severity',[indd-1 0 0],[1 -1 -1]));
                    jj=(j+90)*4+1;
                    ii=i*4+1;
                    MHW(jj:jj+79,ii:ii+119)=mhw;
                end
            end
        else
            % cross 0oE
            %             startlons=[handles.W:60:359 0:60:handles.E];
            %             startlons(startlons>180)=startlons(startlons>180)-360;
            %             startlons=startlons+180;
            cc=0;
            for i=[handles.W:30:359 0:30:handles.E]
                cc=cc+1;
                for j=[handles.S:20:handles.N]
% %                     fn=['/Users/alexsengupta/Desktop/MHWpc98/MHW_SSTA.pc98.sst.avhrr-only-v2.1981-2014.i',num2str(i),'.j',num2str(j),'.nc'];
%                     fn=['/Volumes/HDD_MacRetina/DATASETS2/NOAA OI SST V2 High Resolution Dataset/REGIONAL_timeseries/heatwaves/MHW_SSTA.pc98.sst.avhrr-only-v2.1981-2014.i',num2str(i),'.j',num2str(j),'.nc'];
%                     mhw=squeeze(nc_varget(fn,'mhw',[indd-1 0 0],[1 -1 -1]));
%                     jj=j*4+1;
%                     ii=startlons(cc)*4+1;
%                     MHW(jj:jj+119,ii:ii+239)=mhw;
                    
                    fn=['/Users/z3045790/OneDrive/DATASETS/mhw_98pc/mhw_severity.pc98.',num2str(i),'to',num2str(i+30),'.',num2str(j),'to',num2str(j+20),'.1981.2019.nc'];
                    mhw=squeeze(nc_varget(fn,'severity',[indd-1 0 0],[1 -1 -1]));
                    jj=(j+90)*4+1;
                    ii=startlons(cc)*4+1;
                    MHW(jj:jj+79,ii:ii+119)=mhw;
                end
            end
            %             lon2=lon(721:end)-360;lon2=[lon2 ;lon(1:720)];
            %             lon=lon2;handles.lon=lon; clear lon2
        end
        MHW(MHW==0)=NaN;
        
        mask_area=MHW;mask_area(:)=1;
        mask_intensity=mask_area;
        
        if handles.AreaFilter.Value
            area=grid_area(lat,lon);
            % select minumum area for MHW
            min_area=(str2num(handles.areamin.String)*(2*pi*earthRadius/360))^2;
            
            mhw=MHW;
            mask=mhw; ind=find(mhw==0);mask(:)=1;mask(ind)=NaN;
            mhw(mhw==0)=NaN;
            mhw(isnan(mhw))=0;
            % if both area and intesity filter - apply the intesity filter
            % first
            if handles.IntensityFilter.Value
                ind=find(mhw<str2num(handles.intensitymin.String));
                mhw(ind)=0;
            end
            
            [B,L] = bwboundaries(mhw,'noholes');
            
            % if greater than minimum size add to new B & L and store area and cummulative intensity
            cc=1;
            clear B2 mhw_area cum_intense
            L2=L*0;
            for n=1:length(B)
                ind=find(L==n);
                tmp=sum(area(ind));
                if tmp>min_area
                    L2(ind)=cc;
                    B2{cc}=B{n};
                    mhw_area(cc)=tmp;
                    cum_intense(cc)=sum(mhw(ind).*area(ind));
                    cc=cc+1;
                end
            end
            L2(L2==0)=NaN;
            if exist('mhw_area','var')
                mhw_area=mhw_area/1000/1000; % km2
                cum_intense=cum_intense/1000/1000; %oC.km2
            end
            mask_area=L2./L2;
            
        end
        if handles.IntensityFilter.Value & ~handles.AreaFilter.Value
            ind=find(MHW>str2num(handles.intensitymin.String));
            mask_intensity(:)=NaN;
            mask_intensity(ind)=1;
        end
        
        pcolor(lon,lat,mask_area.*mask_intensity.*MHW);shading flat
        caxis([0 2])
        plotmap
        title(['y=',num2str(handles.myyear), ' m=',num2str(handles.mymonth), ' d=',num2str( handles.myday)])
        colorbar
        %         xlim([handles.W handles.E+30]);ylim([handles.S handles.N+30])
        xW= handles.W; if handles.E<handles.W; xW=xW-360; end
        xlim([xW handles.E+30]);ylim([handles.S handles.N+20])
        
        if ~isempty(handles.X) & handles.AreaFilter.Value
            figure(2);plot(handles.X,handles.Y,'r+','MarkerSize',20,'LineWidth',3)
            %loop through heatwaves to find closest
            clear RNG
            if exist('mhw_area','var')
                for n=1:length(mhw_area)
                    [indj,indi]=find(L2==n);
                    [RNG(n), AZ] = distance(handles.Y,handles.X,mean(lat(indj)),mean(lon(indi)));
                end
                ind=find(RNG==min(RNG)); %select closest
                clear boundary98
                if RNG(ind)<str2num(handles.min_distance.String) % make sure MHW is within certain radius
                    selectedMHWindex=find(L2==ind);
                    boundary98=B2{ind};
                    
                    plot(lon(boundary98(:,2)),lat(boundary98(:,1)),'r')% plot boundary around closest MHW
                    tmparea=num2str(sum(area(selectedMHWindex))/1000/1000/1e6);
                    tmpint=num2str(max(MHW(selectedMHWindex)));
                    tmp90=num2str(prctile(MHW(selectedMHWindex),90));
                    set(handles.stat_area98,'string',tmparea ); %km2
                    set(handles.stat_intensity98,'string', tmpint); %oC
                    set(handles.pc9098,'string', tmp90); %oC
                    %                 disp(['area: ',tmparea,' intensity: ',tmpint,' intensity90pc: ',tmp90,])
                else
                    set(handles.stat_area98,'string', ''); %km2
                    set(handles.stat_intensity98,'string', ''); %oC
                    set(handles.pc9098,'string', ''); %oC
                end
            else
                set(handles.stat_area98,'string', ''); %km2
                set(handles.stat_intensity98,'string', ''); %oC
                set(handles.pc9098,'string', ''); %oC
            end
        end
    end
    
    
    %% SSTA
	if ~isempty(handles.X) & handles.ssta_check.Value  % if SSTA checked
        SSTA=zeros(length(handles.lat),length(handles.lon));
        if handles.W<=handles.E
            for i=handles.W:30:handles.E
                for j=[handles.S:20:handles.N]
%                     fn=['/Users/alexsengupta/Desktop/MHWpc98/MHW_SSTA.pc98.sst.avhrr-only-v2.1981-2014.i',num2str(i),'.j',num2str(j),'.nc'];
%                     fn=['/Volumes/HDD_MacRetina/DATASETS2/NOAA OI SST V2 High Resolution Dataset/REGIONAL_timeseries/heatwaves/MHW_SSTA.pc98.sst.avhrr-only-v2.1981-2014.i',num2str(i),'.j',num2str(j),'.nc'];
%                     mhw=squeeze(nc_varget(fn,'ssta',[indd-1 0 0],[1 -1 -1]));
%                     jj=j*4+1;
%                     ii=i*4+1;
%                     SSTA(jj:jj+119,ii:ii+239)=mhw;
                    
                    fn=['/Users/z3045790/OneDrive/DATASETS/mhw_98pc/mhw_severity.pc98.',num2str(i),'to',num2str(i+30),'.',num2str(j),'to',num2str(j+20),'.1981.2019.nc'];
                    mhw=squeeze(nc_varget(fn,'ssta',[indd-1 0 0],[1 -1 -1]));
                    jj=(j+90)*4+1;
                    ii=i*4+1;
                    SSTA(jj:jj+79,ii:ii+119)=mhw;
                end
            end
        else
            cc=0;
            for i=[handles.W:30:359 0:30:handles.E]
                cc=cc+1;
                for j=[handles.S:20:handles.N]
%                     fn=['/Users/alexsengupta/Desktop/MHWpc98/MHW_SSTA.pc98.sst.avhrr-only-v2.1981-2014.i',num2str(i),'.j',num2str(j),'.nc'];
%                     fn=['/Volumes/HDD_MacRetina/DATASETS2/NOAA OI SST V2 High Resolution Dataset/REGIONAL_timeseries/heatwaves/MHW_SSTA.pc98.sst.avhrr-only-v2.1981-2014.i',num2str(i),'.j',num2str(j),'.nc'];
%                     mhw=squeeze(nc_varget(fn,'ssta',[indd-1 0 0],[1 -1 -1]));
%                     jj=j*4+1;
%                     ii=startlons(cc)*4+1;
%                     SSTA(jj:jj+119,ii:ii+239)=mhw;
                    
                    fn=['/Users/z3045790/OneDrive/DATASETS/mhw_98pc/mhw_severity.pc98.',num2str(i),'to',num2str(i+30),'.',num2str(j),'to',num2str(j+20),'.1981.2019.nc'];
                    mhw=squeeze(nc_varget(fn,'ssta',[indd-1 0 0],[1 -1 -1]));
                    jj=(j+90)*4+1;
                    ii=startlons(cc)*4+1;
                    SSTA(jj:jj+79,ii:ii+119)=mhw;
                end
            end
        end
        SSTA(SSTA==0)=NaN;
        
        figure(3);clf
        if handles.AreaFilter.Value
            pcolor(lon,lat,SSTA);shading flat
            %             colormap(lbmap(21,'BlueRed'))
            colormap([0.3216    0.1882    0.1882
    0.2814    0.1717    0.2381
    0.2412    0.1551    0.2879
    0.2010    0.1385    0.3377
    0.1608    0.1219    0.3875
    0.1206    0.1053    0.4373
    0.0804    0.0887    0.4872
    0.0402    0.0721    0.5370
         0    0.0556    0.5868
         0    0.1230    0.6163
         0    0.1905    0.6458
         0    0.2579    0.6753
         0    0.3254    0.7049
         0    0.3929    0.7344
         0    0.4603    0.7639
         0    0.5278    0.7934
         0    0.5952    0.8229
         0    0.6627    0.8524
         0    0.7302    0.8819
         0    0.7976    0.9115
         0    0.8651    0.9410
         0    0.9325    0.9705
         0    1.0000    1.0000
    0.1111    1.0000    1.0000
    0.2222    1.0000    1.0000
    0.3333    1.0000    1.0000
    0.4444    1.0000    1.0000
    0.5556    1.0000    1.0000
    0.6667    1.0000    1.0000
    0.7778    1.0000    1.0000
    0.8889    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    0.8750
    1.0000    1.0000    0.7500
    1.0000    1.0000    0.6250
    1.0000    1.0000    0.5000
    1.0000    1.0000    0.3750
    1.0000    1.0000    0.2500
    1.0000    1.0000    0.1250
    1.0000    1.0000         0
    0.9909    0.9502         0
    0.9818    0.9003         0
    0.9728    0.8505         0
    0.9637    0.8007         0
    0.9546    0.7509         0
    0.9455    0.7010         0
    0.9364    0.6512         0
    0.9274    0.6014         0
    0.9183    0.5515         0
    0.9092    0.5017         0
    0.9001    0.4519         0
    0.8910    0.4021         0
    0.8819    0.3522         0
    0.8729    0.3024         0
    0.8638    0.2526         0
    0.8547    0.2027         0
    0.8077    0.1852    0.1118
    0.7606    0.1677    0.2235
    0.7136    0.1502    0.3353
    0.6666    0.1327    0.4471
    0.6195    0.1152    0.5588
    0.5725    0.0977    0.6706
    0.5255    0.0802    0.7824
    0.4784    0.0627    0.8941])
            caxis([-4 4])
            plotmap
            title(['y=',num2str(handles.myyear), ' m=',num2str(handles.mymonth), ' d=',num2str( handles.myday)])
            colorbar
            %             xlim([handles.W handles.E+30]);ylim([handles.S handles.N+30])
            xW= handles.W; if handles.E<handles.W; xW=xW-360; end
            xlim([xW handles.E+30]);ylim([handles.S handles.N+20])
            plot(handles.X,handles.Y,'r+','MarkerSize',20,'LineWidth',3)
            
            hold on;
            if exist('boundary90','var'); plot(lon(boundary90(:,2)),lat(boundary90(:,1)),'k') ; end
            if exist('boundary98','var')
                plot(lon(boundary98(:,2)),lat(boundary98(:,1)),'r')% plot boundary around closest MHW
                sstmask=SSTA*NaN; % create mask based on 98pc boundary for closest region
                sstmask(selectedMHWindex)=1;
                ssta=sstmask.*SSTA; % isolate SSTA in closest MHW region
                set(handles.stat_intensitySSTA,'string', num2str(nanmax(ssta(:)))); %oC
                disp(['y=',num2str(handles.myyear), ' m=',num2str(handles.mymonth), ' d=',num2str( handles.myday),'  area: ',tmparea,' intensity: ',tmpint,' intensity90pc: ',tmp90,' intensity SSTA: ',num2str(nanmax(ssta(:)))])
            else
                set(handles.stat_intensitySSTA,'string', '');
            end
        end
    end
    
    
    
    
    
end



% --- Executes just before MHWControlGUI_updated is made visible.
function MHWControlGUI_updated_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MHWControlGUI_updated (see VARARGIN)

fn='/Users/z3045790/OneDrive/DATASETS/mhw_12_2018/mhw_severity.pc90.0to30.-10to10.1981.2018.nc';
time=nc_varget(fn,'time');
dvec=datevec(time+datenum(1,1,1));
% fn='/Volumes/HDD_MacRetina/DATASETS2/NOAA OI SST V2 High Resolution Dataset/sst.day.mean.1981.v2.nc';
fn='/Users/z3045790/OneDrive/My Projects/OceanXtremes/noaa_oi_LonLat_info.nc'
handles.lon=nc_varget(fn,'lon');
handles.lat=nc_varget(fn,'lat');
handles.X=[];
handles.Y=[];

handles.dateindex=dvec(:,1)*10000 + dvec(:,2)*100+dvec(:,3);
MHW=zeros(length(handles.lat),length(handles.lon));

handles.myyear=1982;
handles.mymonth=1;
handles.myday=2;

% plotmap to chose regions
figure(3);clf
plotmap
xlim([0 360]);ylim([-90 90])
set(gca,'xtick',[0:30:360])
set(gca,'ytick',[0:20:180]-90)
grid on

W = input('Western limit? (0 30 60 90 120 150 180 210 240 270 300 330) [0]:');
if isempty(W); W = 0; end
E = input('Eastern limit? (30 60 90 120 150 180 210 240 270 300 330 360) [360]:');
if isempty(E); E = 360; end; %E=E-60; %if E<W; E==W; end
S = input('Southern limit? (-90 -70 -50 -30 -10 10 30 50 70) [-90]:');
if isempty(S); S = -90; end
N = input('Northern limit? (-70 -50 -30 -10 10 30 50 70 90) [90]:');
if isempty(N); N = 90; end; 
N=N-20; 
if N<S; N==S; end
handles.N=N;
handles.S=S;
handles.E=E;
handles.W=W;

MHWControlGUI_updated_plotdata(handles)
% Choose default command line output for MHWControlGUI_updated
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MHWControlGUI_updated wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MHWControlGUI_updated_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function MYyear_Callback(hObject, eventdata, handles)
% hObject    handle to MYyear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MYyear as text
%        str2double(get(hObject,'String')) returns contents of MYyear as a double


% --- Executes during object creation, after setting all properties.
function MYyear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MYyear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MYmonth_Callback(hObject, eventdata, handles)
% hObject    handle to MYmonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MYmonth as text
%        str2double(get(hObject,'String')) returns contents of MYmonth as a double


% --- Executes during object creation, after setting all properties.
function MYmonth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MYmonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MYday_Callback(hObject, eventdata, handles)
% hObject    handle to MYday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MYday as text
%        str2double(get(hObject,'String')) returns contents of MYday as a double


% --- Executes during object creation, after setting all properties.
function MYday_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MYday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GotoDate.
function GotoDate_Callback(hObject, eventdata, handles)

dn=datenum(str2num(get(handles.MYyear,'string')),str2num(get(handles.MYmonth,'string')),str2num(get(handles.MYday,'string')));

if dn<datenum(2016,4,20) & dn>datenum(1982,1,2)
    [handles.myyear,handles.mymonth,handles.myday]=datevec(dn);
    MHWControlGUI_updated_plotdata(handles)
else
    disp('no more days available in dataset')
end

guidata(hObject, handles);
% hObject    handle to GotoDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in FF_day.
function FF_day_Callback(hObject, eventdata, handles)
dn=datenum(handles.myyear,handles.mymonth,handles.myday)+1;
if dn<datenum(2016,4,20) % last availble date
    [handles.myyear,handles.mymonth,handles.myday]=datevec(dn);
    MHWControlGUI_updated_plotdata(handles)
else
    disp('no more days available in dataset')
end

guidata(hObject, handles);
% hObject    handle to FF_day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in FF_month.
function FF_month_Callback(hObject, eventdata, handles)
dn=datenum(handles.myyear,handles.mymonth,handles.myday)+30;
if dn<datenum(2014,5,1) % last availble date
    [handles.myyear,handles.mymonth,handles.myday]=datevec(dn);
    MHWControlGUI_updated_plotdata(handles)
else
    disp('no more days available in dataset')
end

guidata(hObject, handles);
% hObject    handle to FF_month (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in FF_year.
function FF_year_Callback(hObject, eventdata, handles)
dn=datenum(handles.myyear+1,handles.mymonth,handles.myday);
if dn<datenum(2014,5,1) % last availble date
    [handles.myyear,handles.mymonth,handles.myday]=datevec(dn);
    MHWControlGUI_updated_plotdata(handles)
else
    disp('no more days available in dataset')
end

guidata(hObject, handles);
% hObject    handle to FF_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in RW_day.
function RW_day_Callback(hObject, eventdata, handles)
dn=datenum(handles.myyear,handles.mymonth,handles.myday)-1;
if dn>datenum(1982,1,2) % last availble date
    [handles.myyear,handles.mymonth,handles.myday]=datevec(dn);
    MHWControlGUI_updated_plotdata(handles)
else
    disp('no more days available in dataset')
end

guidata(hObject, handles);
% hObject    handle to RW_day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in RW_year.
function RW_year_Callback(hObject, eventdata, handles)
dn=datenum(handles.myyear-1,handles.mymonth,handles.myday);
if dn>datenum(1982,1,2) % last availble date
    [handles.myyear,handles.mymonth,handles.myday]=datevec(dn);
    MHWControlGUI_updated_plotdata(handles)
else
    disp('no more days available in dataset')
end

guidata(hObject, handles);
% hObject    handle to RW_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in RW_month.
function RW_month_Callback(hObject, eventdata, handles)
dn=datenum(handles.myyear,handles.mymonth,handles.myday)-30;
if dn>datenum(1982,1,2) % last availble date
    [handles.myyear,handles.mymonth,handles.myday]=datevec(dn);
    MHWControlGUI_updated_plotdata(handles)
else
    disp('no more days available in dataset')
end

guidata(hObject, handles);
% hObject    handle to RW_month (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AreaFilter.
function AreaFilter_Callback(hObject, eventdata, handles)
MHWControlGUI_updated_plotdata(handles)
% hObject    handle to AreaFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AreaFilter


% --- Executes on button press in Exit.
function Exit_Callback(hObject, eventdata, handles)

close(handles.figure1)
close all
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function areamin_Callback(hObject, eventdata, handles)

% hObject    handle to areamin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of areamin as text
%        str2double(get(hObject,'String')) returns contents of areamin as a double


% --- Executes during object creation, after setting all properties.
function areamin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to areamin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function intensitymin_Callback(hObject, eventdata, handles)
% hObject    handle to intensitymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of intensitymin as text
%        str2double(get(hObject,'String')) returns contents of intensitymin as a double


% --- Executes during object creation, after setting all properties.
function intensitymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intensitymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in IntensityFilter.
function IntensityFilter_Callback(hObject, eventdata, handles)
MHWControlGUI_updated_plotdata(handles)
% hObject    handle to IntensityFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IntensityFilter


% --- Executes on button press in MHWlocation.
function MHWlocation_Callback(hObject, eventdata, handles)
figure(1)
[handles.X,handles.Y] = ginput(1)
MHWControlGUI_updated_plotdata(handles)
guidata(hObject, handles);
% hObject    handle to MHWlocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function stat_area_Callback(hObject, eventdata, handles)
% hObject    handle to stat_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stat_area as text
%        str2double(get(hObject,'String')) returns contents of stat_area as a double


% --- Executes during object creation, after setting all properties.
function stat_area_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stat_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stat_intensity_Callback(hObject, eventdata, handles)
% hObject    handle to stat_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stat_intensity as text
%        str2double(get(hObject,'String')) returns contents of stat_intensity as a double


% --- Executes during object creation, after setting all properties.
function stat_intensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stat_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pc90_Callback(hObject, eventdata, handles)
% hObject    handle to pc90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pc90 as text
%        str2double(get(hObject,'String')) returns contents of pc90 as a double


% --- Executes during object creation, after setting all properties.
function pc90_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pc90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pc98check.
function pc98check_Callback(hObject, eventdata, handles)
% hObject    handle to pc98check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MHWControlGUI_updated_plotdata(handles)
% Hint: get(hObject,'Value') returns toggle state of pc98check


% --- Executes on button press in ssta_check.
function ssta_check_Callback(hObject, eventdata, handles)
% hObject    handle to ssta_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MHWControlGUI_updated_plotdata(handles)
% Hint: get(hObject,'Value') returns toggle state of ssta_check


% --- Executes on button press in pc90check.
function pc90check_Callback(hObject, eventdata, handles)
% hObject    handle to pc90check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MHWControlGUI_updated_plotdata(handles)
% Hint: get(hObject,'Value') returns toggle state of pc90check



function stat_area98_Callback(hObject, eventdata, handles)
% hObject    handle to stat_area98 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stat_area98 as text
%        str2double(get(hObject,'String')) returns contents of stat_area98 as a double


% --- Executes during object creation, after setting all properties.
function stat_area98_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stat_area98 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stat_intensity98_Callback(hObject, eventdata, handles)
% hObject    handle to stat_intensity98 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stat_intensity98 as text
%        str2double(get(hObject,'String')) returns contents of stat_intensity98 as a double


% --- Executes during object creation, after setting all properties.
function stat_intensity98_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stat_intensity98 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pc9098_Callback(hObject, eventdata, handles)
% hObject    handle to pc9098 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pc9098 as text
%        str2double(get(hObject,'String')) returns contents of pc9098 as a double


% --- Executes during object creation, after setting all properties.
function pc9098_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pc9098 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stat_intensitySSTA_Callback(hObject, eventdata, handles)
% hObject    handle to stat_intensitySSTA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stat_intensitySSTA as text
%        str2double(get(hObject,'String')) returns contents of stat_intensitySSTA as a double


% --- Executes during object creation, after setting all properties.
function stat_intensitySSTA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stat_intensitySSTA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min_distance_Callback(hObject, eventdata, handles)
% hObject    handle to min_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_distance as text
%        str2double(get(hObject,'String')) returns contents of min_distance as a double


% --- Executes during object creation, after setting all properties.
function min_distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonEPS.
function buttonEPS_Callback(hObject, eventdata, handles)
% hObject    handle to buttonEPS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try;print(['images/ssta_',handles.MYday.String,'_',handles.MYmonth.String,'_',handles.MYyear.String,'.eps'],'-f3','-depsc');disp('Saved');end


% --- Executes on button press in FF_7day.
function FF_7day_Callback(hObject, eventdata, handles)
dn=datenum(handles.myyear,handles.mymonth,handles.myday)+7;
if dn<datenum(2016,4,20) % last availble date
    [handles.myyear,handles.mymonth,handles.myday]=datevec(dn);
    MHWControlGUI_updated_plotdata(handles)
else
    disp('no more days available in dataset')
end

guidata(hObject, handles);
% hObject    handle to FF_7day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
