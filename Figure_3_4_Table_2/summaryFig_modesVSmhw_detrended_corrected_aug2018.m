clear; close all
% code loads in 2degree MHW data and composites the number of MHW days for
% different climte modes , above and below a threshold (thresh)
% Also calculates significance based on a 500 iteration monte carlo at each
% location (using sythetic modes based on an AR(4) process)
% Claculates based on both a 5-95% and 1-99% significance level

%% set the threshold
thresh='zero'; 
% thresh='std05';
%%
                            
% load monthly indices
[time,index(1,:),index(2,:),index(3,:),index(4,:),index(5,:),index(6,:),index(7,:),index(8,:),index(9,:),index(10,:)]=load_modes(1982,2015);
index_name={'nao','amo','nino34','pdo','tpi','anino','sam','modoki','dmi','npgo'};
index(4,:)=smooth(index(4,:),120,'rloess'); % apply decadal smoothing to PDO
index(5,:)=smooth(index(5,:),120,'rloess'); % apply decadal smoothing to TPI

% normalize timeseries
for m=1:length(index_name)
    index(m,:)=normalise1D(index(m,:));
    %detrend index
    tmp=index(m,:);
    ind=find(~isnan(tmp));
    index(m,ind)=detrend(tmp(ind));
end
% remove AMO
index(2,:)=[];
index_name(2)=[];

figure(1);clf
for m=1:length(index_name)
    subplot(3,3,m)
    plot(time,index(m,:))
    xlim([1982 2015])
    title(index_name{m})
end


files=dir('/Users/alexsengupta/Desktop/MHW_new/2degree/mhw_severity.pc90.*.nc')


for f=1:length(files);
    % load SSTA and MHW data based on detrended SST in Eric's code
    fn=['/Users/alexsengupta/Desktop/MHW_new/2degree/',files(f).name];
    fno=['/Users/alexsengupta/CCRC/My Projects/OceanXtremes/Holbrook- Drivers/code2018/processed_data/detrended_processed_',files(f).name,'_',thresh,'.mat'];
    disp(f)
    disp(fno)
    if ~exist(fno,'file')
        ind=findstr('.',fn);
        mylon=fn(ind(2)+1:ind(3)-1);
        mylat=fn(ind(3)+1:ind(4)-1);
        ind=findstr('to',mylon);
        lon1=str2num(mylon(1:ind-1));
        ind=findstr('to',mylat);
        lat1=str2num(mylat(1:ind-1));
        ii=lon1/30*15+1;
        jj=(lat1+90)/20*10+1;
        
        mhwREG=nc_varget(fn,'severity');
        %         sstaREG=nc_varget(fn,'ssta');
        lat=nc_varget(fn,'lat');
        lon=nc_varget(fn,'lon');
        rtime=nc_varget(fn,'time');
        % "days since 01-01-01 00:00:00" ;
        
        dvec=datevec(datenum(1,1,0)+rtime);
        
        fulltime=dvec(:,1)+(dvec(:,2)-1)/12+dvec(:,3)/365;
        indx=find(fulltime>time(1) & fulltime<=time(end));
        fulltime=fulltime(indx);
        
        %%
        matr_modePos_mhwInc=zeros(length(index_name),length(lat),length(lon));
        matr_modeNeg_mhwInc=zeros(length(index_name),length(lat),length(lon));
        matr_modePos_mhwDec=zeros(length(index_name),length(lat),length(lon));
        matr_modeNeg_mhwDec=zeros(length(index_name),length(lat),length(lon));
        
        matr_modePos_mhwInc99=zeros(length(index_name),length(lat),length(lon));
        matr_modeNeg_mhwInc99=zeros(length(index_name),length(lat),length(lon));
        matr_modePos_mhwDec99=zeros(length(index_name),length(lat),length(lon));
        matr_modeNeg_mhwDec99=zeros(length(index_name),length(lat),length(lon));
        medianprop           =zeros(length(index_name),length(lat),length(lon));
        
        for i=1:length(lon)
            tic
            disp(num2str(i))
            for j=1:length(lat)
                
                mhw=squeeze(mhwREG(indx,j,i));
                %ssta=squeeze(sstaREG(:,j,i));
                if ~all(isnan(mhw))
                    
                    %                     ind=find(isnan(mhw));%ssta_mhw=ssta;ssta_mhw(ind)=NaN;
                    
                    for m=1:length(index_name)
                        index_interp=interp1(time,index(m,:),fulltime);
                        if strcmp(thresh,'zero'); threshold=0; end
                        if strcmp(thresh,'std05'); threshold=nanstd(index_interp)/2; end
                        NT=length(mhw);
                        
                        ind =find(index_interp>thresh);
                        tmp=mhw(ind);
                        total_days_pos=nansum(tmp./tmp);
                        
                        ind =find(index_interp<-thresh);
                        tmp=mhw(ind);
                        total_days_neg=nansum(tmp./tmp);
                        
                        % Monte Carlo
                        arcoeffs = arcov(index_interp(~isnan(index_interp)),4);
                        for n=1:500
                            y = filter(1,arcoeffs,randn(length(fulltime),1));y=normalise1D(y);
                            %         y = randn(length(fulltime),1);y=normalise1D(y);
                            ind =find(y>0);
                            tmp=mhw(ind);
                            MC_total_days_pos(n)=nansum(tmp./tmp);
                            MC_total_days_pos_prop(n)=nansum(tmp./tmp)/length(ind);
                        end
                        
                        % proportion of days with +/- mode in MHW
                        pc=prctile(MC_total_days_pos_prop,[5 95]); 
                        pc99=prctile(MC_total_days_pos_prop,[1 99]); 
                        
                        prop_modePos=total_days_pos/length(find(index_interp>thresh));
                        prop_modeNeg=total_days_neg/length(find(index_interp<-thresh));
                        
                        %% 5-95 percentiles
                        % positive mode leads to sig increase in mhw
                        if prop_modePos>pc(2)
                            matr_modePos_mhwInc(m,j,i)=prop_modePos;
                        end
                        % negative mode leads to sig increase in mhw
                        if prop_modeNeg>pc(2)
                            matr_modeNeg_mhwInc(m,j,i)=prop_modeNeg;
                        end
                        % positive mode leads to sig decrease in mhw
                        if prop_modePos<pc(1)
                            matr_modePos_mhwDec(m,j,i)=prop_modePos;
                        end
                        % negative mode leads to sig decrease in mhw
                        if prop_modeNeg<pc(1)
                            matr_modeNeg_mhwDec(m,j,i)=prop_modeNeg;
                        end
                        %% 1-99 percentiles
                        % positive mode leads to sig increase in mhw
                        if prop_modePos>pc99(2)
                            matr_modePos_mhwInc99(m,j,i)=prop_modePos;
                        end
                        % negative mode leads to sig increase in mhw
                        if prop_modeNeg>pc99(2)
                            matr_modeNeg_mhwInc99(m,j,i)=prop_modeNeg;
                        end
                        % positive mode leads to sig decrease in mhw
                        if prop_modePos<pc99(1)
                            matr_modePos_mhwDec99(m,j,i)=prop_modePos;
                        end
                        % negative mode leads to sig decrease in mhw
                        if prop_modeNeg<pc99(1)
                            matr_modeNeg_mhwDec99(m,j,i)=prop_modeNeg;
                        end
                        medianprop(m,j,i)=median(MC_total_days_pos_prop);
                    end
                end
                
            end
            toc
        end
        
%         matr_modePos_mhwInc_all(:,jj:jj+9,ii:ii+14)=matr_modePos_mhwInc;
%         matr_modeNeg_mhwInc_all(:,jj:jj+9,ii:ii+14)=matr_modeNeg_mhwInc;
%         matr_modePos_mhwDec_all(:,jj:jj+9,ii:ii+14)=matr_modePos_mhwDec;
%         matr_modeNeg_mhwDec_all(:,jj:jj+9,ii:ii+14)=matr_modeNeg_mhwDec;
%         
%         matr_modePos_mhwInc99_all(:,jj:jj+9,ii:ii+14)=matr_modePos_mhwInc99;
%         matr_modeNeg_mhwInc99_all(:,jj:jj+9,ii:ii+14)=matr_modeNeg_mhwInc99;
%         matr_modePos_mhwDec99_all(:,jj:jj+9,ii:ii+14)=matr_modePos_mhwDec99;
%         matr_modeNeg_mhwDec99_all(:,jj:jj+9,ii:ii+14)=matr_modeNeg_mhwDec99;
%         medianprop_all(:,jj:jj+9,ii:ii+14)=medianprop;

        save(fno,'matr_modePos_mhwInc','matr_modeNeg_mhwInc','matr_modePos_mhwDec','matr_modeNeg_mhwDec'...
            ,'matr_modePos_mhwInc99','matr_modeNeg_mhwInc99','matr_modePos_mhwDec99','matr_modeNeg_mhwDec99' ,'index_name','medianprop')
    end
    
end

