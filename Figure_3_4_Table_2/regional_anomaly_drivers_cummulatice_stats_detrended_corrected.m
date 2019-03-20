clear; close all

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

figure(1);clf
for m=1:length(index_name)
    subplot(5,2,m)
    plot(time,index(m,:))
    xlim([1982 2015])
    title(index_name{m})
end
% remove AMO
index(2,:)=[];
index_name(2)=[];


files=dir('/Users/alexsengupta/CCRC/My Projects/OceanXtremes/Holbrook- Drivers/corrected_drivers_code/processed_data/regional_ssta_severity*');
for reg=1:length(files)
    reg
    fn=files(reg).name;
    ind=findstr('.',fn);
    ind_=findstr('_',fn);
    region=fn(ind(2)+1:ind_(3)-1);
    REGION{reg}=region;

    fn=['/Users/alexsengupta/CCRC/My Projects/OceanXtremes/Holbrook- Drivers/corrected_drivers_code/processed_data/',fn];
    
    mhwtime=nc_varget(fn,'time');
    mhw=nc_varget(fn,'severity');
    ssta=nc_varget(fn,'ssta');
    dvec=datevec(datenum(1,1,0)+mhwtime);
    fulltime=dvec(:,1)+(dvec(:,2)-1)/12+dvec(:,3)/365;
    
    ind=find(isnan(mhw));ssta_mhw=ssta;ssta_mhw(ind)=NaN;
    
%     figure(8);clf
%     fill2c(fulltime,ssta,[.7 .7 .9],[.9 .7 .7]);hold on
%     fill2c(fulltime,ssta_mhw,'b',[.8 0 0])
%     xlim([1982 2015])
%     set(gca,'xtick',[1982:2:2015]);grid on;title(['region: ',region])
%     title(region)
%     pbaspect([4 1 1])
%     ylim([-max(ssta_mhw) max(ssta_mhw)])
%     set(gca,'fontsize',8)
    
%     ind=find(~isnan(ssta_mhw));
%     ind=diff(ind);
%     MHW_MEANDURATION(reg)=length(find(ind==1));
%     ind=find(ind>1);
%     MHW_COUNT(reg)=length(ind);
%     MHW_MEANDURATION(reg)=MHW_MEANDURATION(reg)/MHW_COUNT(reg);
    
    % Composite analysis. Look at look at MHW characteristics when indices are
    % +ve or -ve
    for m=1:length(index_name)
        index_interp=interp1(time,index(m,:),fulltime);

%         figure(10+m);clf
%         fill2c(fulltime,ssta,[.7 .7 .9],[.9 .7 .7]);hold on
%         fill2c(fulltime,ssta_mhw,'b',[.8 0 0])
%         xlim([1982 2015])
%         set(gca,'xtick',[1982:2015]);grid on;title(['region: ',region])
%         plot(fulltime,index_interp/nanstd(index_interp),'k')
%         title(index_name{m})
%         pbaspect([4 1 1])
        
        ind =find(index_interp>0);
        tmp=mhw(ind);
        total_days_pos=nansum(tmp./tmp);
        total_days_pos_prop=total_days_pos/length(ind);
%         tmp=ssta_mhw(ind);
%         total_icum_pos=nansum(tmp);
        
        %         %large events
        %         ind =find(normalise1D(index_interp)>1);
        %         tmp=mhw(ind);
        %         total_days_pos_large=nansum(tmp./tmp);
        
        ind =find(index_interp<0);
        tmp=mhw(ind);
        total_days_neg=nansum(tmp./tmp);
        total_days_neg_prop=total_days_neg/length(ind);
%         tmp=ssta_mhw(ind);
%         total_icum_neg=nansum(tmp);
        
        %         %large events
        %         ind =find(normalise1D(index_interp)<-1);
        %         tmp=mhw(ind);
        %         total_days_neg_large=nansum(tmp./tmp);
        
        % Monte Carlo
        arcoeffs = arcov(index_interp(~isnan(index_interp)),4);
        for n=1:1000
            y = filter(1,arcoeffs,randn(length(fulltime),1));y=normalise1D(y);
            %         y = randn(length(fulltime),1);y=normalise1D(y);
            ind =find(y>0);
            tmp=mhw(ind);
            MC_total_days_pos(n)=nansum(tmp./tmp);
            MC_total_days_pos_prop(n)=nansum(tmp./tmp)/length(ind);
%             tmp=ssta_mhw(ind);
%             MC_total_icum_pos(n)=nansum(tmp);
            
            %             %large
            %             ind =find(y>1);
            %             tmp=mhw(ind);
            %             MC_total_days_pos_large(n)=nansum(tmp./tmp);
            %             MC_total_days_pos_prop_large(n)=nansum(tmp./tmp)/length(ind);
            
            % No need as positive and negative are symmetrical
            %         ind =find(y<0);
            %         tmp=mhw(ind);
            %         MC_total_days_neg(n)=nansum(tmp./tmp);
            %         tmp=ssta_mhw(ind);
            %         MC_total_icum_neg(n)=nansum(tmp);
        end
        
        total_days_pos_prop_all(reg,m)=total_days_pos_prop;
        total_days_neg_prop_all(reg,m)=total_days_neg_prop;
        total_days_prop_MC_all(reg,m,:)=prctile(MC_total_days_pos_prop,[5 50 95]);
        clear *prop MC* *pos
    end
end
  

save 'processed_data/regional_proportions_plusMC.mat' index_name REGION total_days_pos_prop_all total_days_neg_prop_all total_days_prop_MC_all
%         
%         figure(50+reg);
%         subplot(2,1,1)
%         %         pc=prctile(MC_total_days_pos,[5 95]);plot([m m],pc,'color',[0.5 0.5 0.5]);
%         %         hold on
%         %         plot([m+0.15 m-0.15],[pc(1) pc(1)],'color',[0.5 0.5 0.5]);plot([m+0.15 m-0.15],[pc(2) pc(2)],'color',[0.5 0.5 0.5])
%         %         plot(m, total_days_pos,'ko','markerfacecolor','r');
%         %         plot(m, total_days_neg,'ko','markerfacecolor','k')
%         % proportion of days with +/- mode in MHW
%         
%         plot([m m],pc,'color',[0.9 0.5 0.5]);
%         hold on
%         plot([m+0.15 m-0.15],[pc(1) pc(1)],'color',[0.9 0.5 0.5]);
%         plot([m+0.15 m-0.15],[pc(2) pc(2)],'color',[0.9 0.5 0.5])
%         plot([m+0.05 m-0.05],median(MC_total_days_pos_prop)*[1 1],'color',[0.7 0.3 0.3])
%         plot(m, total_days_pos/length(find(index_interp>0)),'ko','markerfacecolor','r');
%         plot(m, total_days_neg/length(find(index_interp<0)),'ko','markerfacecolor','k')
%         
        %         %large
        %         pc=prctile(MC_total_days_pos_prop_large,[5 95]);
        %         plot([m+.2 m+.2],pc,'color',[0.9 0.5 0.5]);
        %         hold on
        %         plot(m+.2, total_days_pos_large/length(find(index_interp>1)),'xr');
        %         plot(m+.2, total_days_neg_large/length(find(index_interp<1)),'xk')
        
        
%         psig=1;nsig=1;
%         if total_days_pos/length(find(index_interp>0))>max(pc) | total_days_pos/length(find(index_interp>0))<min(pc)
%             psig=-1
%         end
%         if total_days_neg/length(find(index_interp<0))>max(pc) | total_days_neg/length(find(index_interp<0))<min(pc)
%             nsig=-1
%         end
%         table_days(reg,m,1)=total_days_pos/length(find(index_interp>0))*psig;
%         table_days(reg,m,2)=total_days_neg/length(find(index_interp<0))*nsig;
%         
%         subplot(2,1,2)
%         pc=prctile(MC_total_icum_pos,[5 95]);plot([m m],pc,'color',[0.5 0.5 0.5]);
%         hold on
%         plot([m+0.15 m-0.15],[pc(1) pc(1)],'color',[0.5 0.5 0.5]);plot([m+0.15 m-0.15],[pc(2) pc(2)],'color',[0.5 0.5 0.5])
%         plot(m, total_icum_pos,'o','markerfacecolor','r');
%         plot(m, total_icum_neg,'ko','markerfacecolor','k')
%         
%         psig=1;nsig=1;
%         if total_icum_pos>max(pc) | total_icum_pos<min(pc)
%             psig=-1;
%         end
%         if total_icum_neg>max(pc) | total_icum_neg<min(pc)
%             nsig=-1
%         end
%         table_icum(reg,m,1)=total_icum_pos*psig;
%         table_icum(reg,m,2)=total_icum_neg*nsig;
%         
%     end
%     figure(50+reg)
%     subplot(2,1,1)
%     set(gca,'xtick',1:9,'xticklabel',index_name)
%     rotateticklabel(gca,90);
%     xlim([0 10])
%     pbaspect([3 1 1])
%     title(['Days in heatwave: ',region])
%     
%     subplot(2,1,2)
%     set(gca,'xtick',1:9,'xticklabel',index_name)
%     rotateticklabel(gca,90);
%     xlim([0 10])
%     pbaspect([3 1 1])
%     title(['Cummulative intensity: ',region])
%     
%     print(['mkw_v_index_',region],'-depsc')
%     clear sst* mhw* total* MC* pc
% end
% 
% 
% for reg=1:length(files)
%     fn=files(reg).name;
%     ind=findstr('.',fn);
%     ind_=findstr('_',fn);
%     region=fn(ind(2)+1:ind_(2)-1);
%     disp([region,' ',num2str(table_days(reg,:,1))])
%     disp([region,' ',num2str(table_days(reg,:,2))])
% end
% 
% for reg=1:length(files)
%     fn=files(reg).name;
%     ind=findstr('.',fn);
%     ind_=findstr('_',fn);
%     region=fn(ind(2)+1:ind_(2)-1);
%     disp([region,' ',num2str(table_icum(reg,:,1))])
%     disp([region,' ',num2str(table_icum(reg,:,2))])
% end
