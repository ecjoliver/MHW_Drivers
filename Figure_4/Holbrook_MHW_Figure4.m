%Jessica Benthuysen
%29 Oct 2018
%Table 3 as a bar plot for those modes than enhance or reduce MHWs based on percentage of days.

%1 Benguela
%2 Leeuwin
%3 California
%4 Iberian/Canary
%5 Humboldt/Peru
%
%6 Gulf Stream
%7 Kuroshio
%8 Brazil-Malvinas Confluence
%9 Agulhas %
%10 Agulhas retroflection
%11 EAC
%12 EAC extension
%
%13 GBR
%14 Seychelles Is.
%15 Galapagos
%16 Bay of Bengal
%17 Caribbean
%
%18 Mediterranean
%19 Bering Sea
%20 NW Atlantic
%21 NE Pacific
%22 South Central Pacific 

%Cases: stores the values corresponding to the original 'Table 3', which has the percentage of days that the region experienced MHW conditions
%during a positive or negative phase of a given mode. 
%Positive values are where numbers are red (enhanced) and negative values are where numbers are blue (suppressed MHWs).
%The negative indicates that the percentage of days are suppressed but the values themselves are positive.

%Columns are listed below by positive and negative phases. Rows numbers correspond to the regions.
Cases = [-5.4	8.9	0	0	0	0	0	0	13.3	-0.5	0	0	0	10.3	0	0	0	0;   ... %1
         0	0	-2.0	13.5	-3.6	12.2	-2.9	12.0	0	0	9.6	-5.8	-2.3	14.9	0	0	11.2	-4.1; ...%2
	 0	0	16.4	-3.8	0	0	0	0	0	0	0	0	0	0	0	0	-5.2	15.3; ... %3
	 0	0	0	0	0	0	0	0	0	0	0	0	0	0	-4.4	0	0	0; ... %4
	 0	0	12.8	-4.7	0	0	0	0	0	0	0	0	0	0	11.8	0	0	0; ... %5
	 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ... %6
	 0	0	0	0	0	13.0	0	0	10.7	0	10.0	-5.8	0	0	0	0	11.3	-4.4; ... %7
	 0	0	0	0	0	0	0	0	-5.3	11.1	0	0	0	0	0	0	-5.7	10.4; ... %8
	 9.0	-5.5	0	0	0	0	9.5	0	0	0	0	0	0	0	0	0	0	0; ... %9
	 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ... %10
	 -5.2	11.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ... %11
	 0	0	0	0	0	0	0	0	0	-5.8	0	0	0	0	-6.1	0	0	0; ... %12
	 -4.5	10.3	0	0	0	0	0	0	0	0	0	-5.9	0	10.4	0	0	0	0; ... %13
	 0	0	11.5	-3.5	-4.4	10.8	0	0	0	9.7	0	0	0	0	12.3	-3.7	0	0; ... %14
	 0	0	18.0	-2.1	0	0	0	0	0	0	0	0	0	15.7	15.5	-5.7	0	0; ... %15
	 -4.6	9.8	0	0	-4.6	9.8	-3.8	10.1	0	0	9.0	-5.2	0	0	10.0	-4.9	10.1	0; ... %16
	 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ... %17
	 0	0	0	0	0	0	0	0	0	0	0	0	0	-4.9	10.0	-4.8	0	0; ... %18
	 0	0	0	0	14.4	-3.7	13.8	-5.1	-5.8	12.0	0	0	0	0	0	0	0	0; ... %19
	 -2.9	9.8	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ... %20
	 0	0	17.2	-5.1	17.5	-4.5	18.4	-4.8	-4.9	0	0	0	16.1	-4.5	0	0	-4.8	17.9; ... %21
	 0	0	14.0	-2.6	0	0	13.1	-4.0	0	0	-4.8	12.0	11.5	-4.1	0	0	0	0]; %22


yticklabels = {'Benguela'; 'Leeuwin'; 'California'; 'Iberian/Canary'; 'Humboldt/Peru'; ...
               'Gulf Stream'; 'Kuroshio'; 'Brazil-Malvinas'; 'Agulhas'; 'Agulhas retroflection'; 'EAC'; 'EAC extension'; ...
	       'Great Barrier Reef'; 'Seychelles'; 'Galapagos'; 'Bay of Bengal'; 'Caribbean'; ...
	       'Mediterranean'; 'Bering'; 'NW Atlantic'; 'NE Pacific'; 'S. Central Pacific'}; 


%18 Colors are set by modes (plus and then minus): NAO (+/-), Nino34, PDO, TPI/IPO, ATLN1, SAM, EMI, DMI, NPGO
Colors = [206 35 43;  229 119 118; 26 144 68; 108 184 119; 55 62 132; 112 115 168; ... %NAO NINO34 PDO
          169 174 58; 244 237 87; 145 59 132; 188 127 171; 32 116 160; 92 188 184; ... %TPI/IPO ANINO SAM
	  99 99 99;   203 203 201; 161 74 44; 230 131 90; 60 27 58; 117 71 110]/255;   %EMI DMI NPGO


%first column is the median number of days, overall proportion of days
CasesM = [7.2; 7.7; 10.2; 6.3; 8.7; 7.8; 7.9; 8.0; 7.2; 8.2; 8.5; 8.5; 7.4; 7.5; 10; 7.2; 7.3; 7.2; 9.2; 6.5; 11.2; 8.2]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configure bar positions and colours in the four regions.

barwidthM = [3; 6; 2; 1; 2; ...
             1; 4; 2; 2; 1; 1; 2; ...
	     3; 4; 3; 6; 1; ...
	     2; 3; 1; 6; 4]; 

barwidthP = [3; 6; 2; 0; 2; ...
             0; 4; 2; 2; 0; 1; 0; ...
	     2; 4; 3; 6; 0; ...
	     2; 3; 1; 6; 4]; 

%Include in the blank space in the below barwidth
barwidthN = [3; 6; 2; 1; 2; ... %EBCs
             1; 4; 2; 2; 1; 1; 2; ... %WBCs
	     3; 4; 3; 6; 1; ... %Tropics
	     2; 3; 1; 6; 4];  %MHLs

for IS = 1:4;
  if IS == 1; II = 1:5; elseif IS == 2; II = 6:12; 
  elseif IS == 3; II = 13:17; elseif IS == 4; II = 18:22; 
  end; 
  cntM = 0;   cntP = 0;   cntN = 0; 
  for iy = II
    xgridM{iy} = [1:barwidthM(iy)] + cntM; 
    cntM = cntM + barwidthM(iy) + 1; 
    xgridP{iy} = [1:barwidthP(iy)] + cntP; 
    cntP = cntP + barwidthM(iy) + 1; 
    xgridN{iy} = [1:barwidthN(iy)] + cntN; 
    cntN = cntN + barwidthM(iy) + 1; 
  end; 
end;

clear CASES CASES2 IXtotal IXtotal2 xticks

%Set-up correct spacing for the colored bars.
for IS = 1:4
  if IS == 1; II = 1:5; %EBCs
  elseif IS == 2; II = 6:12; %WBCs
  elseif IS == 3; II = 13:17; %Tropics
  elseif IS == 4; II = 18:22; %MHLs
  end; 
  for iy = II;  
    %if ismember(iy,[13]); addextra = 1; else; addextra = 0; end;  
    ii = find(Cases(iy,:) > 0);  %ii_enh{iy} = ii; 
    %if isempty(ii) == 0; ix = cnt:cnt + length(ii) - 1;
    %else; ix = cnt; %only one position, with median 
    %end;
    if isempty(ii) == 0; 
      CASES{iy} = Cases(iy,ii); 
      for ixx = 1:length(ii); COLORS{iy}(ixx,:) = Colors(ii(ixx),:); end;  
    end; 
    %Store values for the suppressed phase
    ii2 = find(Cases(iy,:) < 0);  %ii_supp{iy} = ii2; 
    %ix2 = cnt + 1:cnt + length(ii2); %add one so they are not plotted under the median
    if isempty(ii2) == 0
      CASES2{iy} = -Cases(iy,ii2); 
      for ixx = 1:length(ii2); COLORS2{iy}(ixx,:) = Colors(ii2(ixx),:); end;  
      if iy == 1; %Benguela, for 3rd bar; add in blank column; 
        CASES2{iy} = [-Cases(iy,ii2) 0]; COLORS2{iy} = [COLORS2{iy}; [1 1 1]]; 
      elseif iy == 5; %Humboldt/Peru, for 2nd bar- add in blank column; 
        CASES2{iy} = [-Cases(iy,ii2) 0]; COLORS2{iy} = [COLORS2{iy}; [1 1 1]];
      elseif iy == 7; %Kuroshio, for 1rst and 2nd bar- add in blank column; 
        CASES2{iy} = [0 0 -Cases(iy,ii2)]; COLORS2{iy} = [[1 1 1]; [1 1 1]; COLORS2{iy}]; 
      elseif iy == 9; %Agulhas, for 2nd bar- add in blank column; 
        CASES2{iy} = [-Cases(iy,ii2) 0]; COLORS2{iy} = [COLORS2{iy}; [1 1 1]]; 
      elseif iy == 13; %GBR, for 2nd bar- add in blank column; 
        CASES2{iy} = [-Cases(iy,ii2(1)) 0 -Cases(iy,ii2(2))]; COLORS2{iy} = [COLORS2{iy}(1,:); [1 1 1]; COLORS2{iy}(2,:)]; 
      elseif iy == 14; %Seychelles, for 3rd bar- add in blank column; 
        CASES2{iy} = [-Cases(iy,ii2(1:2)) 0 -Cases(iy,ii2(3))]; COLORS2{iy} = [COLORS2{iy}(1:2,:); [1 1 1]; COLORS2{iy}(3,:)]; 
      elseif iy == 15; %Galapagos, for 2nd bar- add in blank column; 
        CASES2{iy} = [-Cases(iy,ii2(1)) 0 -Cases(iy,ii2(2))]; COLORS2{iy} = [COLORS2{iy}(1,:); [1 1 1]; COLORS2{iy}(2,:)]; 
      elseif iy == 16; %Bay of Bengal, for 6th bar- add in blank column; 
        CASES2{iy} = [-Cases(iy,ii2) 0]; COLORS2{iy} = [COLORS2{iy}; [1 1 1]];
      elseif iy == 18; %Mediterranean, switch order of bars so they are consistent with the upper bar
        CASES2{iy} = [-Cases(iy,flip(ii2))]; COLORS2{iy} = [COLORS2{iy}(2,:); COLORS2{iy}(1,:)];
      elseif iy == 21; %NE Pacific, switch order of bars so they are consistent with the upper bars
        CASES2{iy} = [-Cases(iy,ii2([1:3 5 6 4]))]; COLORS2{iy} = [COLORS2{iy}([1:3 5 6 4],:)];
      end; 
    end; %if, isempty
    if iy == 6 | iy == 10 | iy == 17; %Gulf Stream | Agulghas retroflection | Caribbean , for 1rst bar- add in blank column; 
    %Originally empty, added here or else only the line will show up and then because of the order of 'hold on',
    %it will ensure that there are not blank spaces below the other bars along the y = 0 line. 
        CASES2{iy} = [0]; COLORS2{iy} = [1 1 1];
    end; 
    %Store values for up to the suppressed phase or the median.
    ii3 = 1:length(CASES2{iy}); %max([length(ii) length(ii2)]); %ii_total{iy} = 1:ii3; 
    for ixx = ii3
      if CASES2{iy}(ixx) == 0; CASESW{iy}(ixx) = CasesM(iy); 
      else; CASESW{iy}(ixx) = CASES2{iy}(ixx); 
      end; %if 
    end; %for, ixx
    %COLORSW{iy} = ones(length(ii3), 3); %white
  end;  %for, iy
end; %for, IS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure(1); clf; 
%set(f1, 'position', [2606 10 560 840]); %4 x 1 original
set(f1, 'position', [2606 100 1200 560]); 

for IS = 1:4; 
  if IS == 1; II = 1:5; %EBCs
    titlename = '(a)  Eastern Boundary Currents'; 
    xticklabels = {'Benguela'; 'Leeuwin'; 'California        '; 'Iberian/Canary'; '                Humboldt/Peru'}; 
  elseif IS == 2; II = 6:12; %WBCS
    titlename = '(b)  Western Boundary Currents'; 
    xticklabels = {'Gulf Stream'; 'Kuroshio'; 'Brazil-Malvinas'; 'Agulhas'; '  Agulhas r.f.'; '    EAC'; '    EAC extension'}; 
  elseif IS == 3; II = 13:17; %Tropics
    titlename = '(c)  Tropics'; 
    xticklabels = {'Great Barrier Reef'; 'Seychelles'; 'Galapagos'; 'Bay of Bengal'; 'Caribbean'}; 
  elseif IS == 4; II = 18:22; %MHLS
    titlename = '(d)  Middle and High Latitudes';
    xticklabels = {'Mediterranean'; 'Bering'; 'NW Atlantic'; 'NE Pacific'; 'S. Central Pacific retroflection'}; 
  end; %IS

  ss1 = subplot(2,2,IS); cla;
  pos1 = get(ss1, 'position'); 
  set(ss1, 'position', [pos1(1) pos1(2) pos1(3) pos1(3)/1.3]); %4 x 1 size: [.13 .7673 .7750 .1577]

  for iy = II;
    %Plot colored bar up to enhanced value.
    if isempty(CASES{iy}) == 0; 
      for ixx = 1:length(CASES{iy}); 
        b1 = bar(xgridP{iy}(ixx),CASES{iy}(ixx)); hold on; 
        b1.FaceColor = COLORS{iy}(ixx,:); b1.EdgeColor = 'none'; b1.BarWidth = 1;
      end; %for, ixx
    end; %if, isempty

    %Plot colored bar (for suppressed values) up to median line. 
    if isempty(CASES2{iy}) == 0; 
      for ixx = 1:length(CASES2{iy}); 
        b1 = bar(xgridN{iy}(ixx),CasesM(iy)); %bar to the median line 
        b1.FaceColor = COLORS2{iy}(ixx,:); b1.EdgeColor = 'none'; b1.BarWidth = 1;
      end; %for, ixx
    end; %if, isempty

    %Plot white bar from 0 up to suppressed value or up to median.
    if isempty(CASESW{iy}) == 0; 
      for ixx = 1:length(CASESW{iy}); 
        b1 = bar(xgridN{iy}(ixx),CASESW{iy}(ixx)); hold on; 
        b1.FaceColor = [1 1 1]; b1.EdgeColor = 'none'; b1.BarWidth = 1;
      end; %for, ixx
    end; %if, isempty
  
    %Plot median line
    pp1 = line([(xgridN{iy}(1) - 1/2) (xgridN{iy}(end) + 1/2)], [CasesM(iy) CasesM(iy)]); hold on;  
    pp1.Color = 'k'; pp1.LineWidth = 1;  
    
    xticks(iy) = mean(xgridN{iy}); 
  end; %for, iy
  %xlim([-.3 22.3]);

  set(gca, 'ytick', 0:5:20); set(gca, 'tickdir', 'out');  
  ylim([0 20]); ylabel('percentage of days');  
  tc = title(titlename, 'fontweight', 'normal');
  tcpos = get(tc, 'position'); set(gca, 'fontsize', 8); 
  set(tc, 'position', [tcpos(1) tcpos(2)+.6 0]); 
  set(gca, 'xtick', xticks(II), 'xticklabel', xticklabels); 
  set(gca, 'tickdir', 'out', 'ticklength', [.005 .005]);
end; %IS

print -f1 -depsc -r600 Figure4_Holbrook_original %before axes-formatting and legend are included 
print -f1 -dpng -r600 Figure4_Holbrook_original %before axes-formatting and legend are included 

%Note: Figure4_Holbrook_formatted includes the formatted text for the axes and legend





























%For some reason the colors are not completely/consistenly being output correctly...
%Need to check the colors manually.















yticklabels = {'Benguela'; 'Leeuwin'; 'California'; 'Iberian/Canary'; 'Humboldt/Peru'; ...
               'Gulf Stream'; 'Kuroshio'; 'Brazil-Malvinas'; 'Agulhas'; 'Agulhas retroflection'; 'EAC'; 'EAC extension'; ...
	       'Great Barrier Reef'; 'Seychelles'; 'Galapagos'; 'Bay of Bengal'; 'Caribbean'; ...
	       'Mediterranean'; 'Bering'; 'NW Atlantic'; 'NE Pacific'; 'S. Central Pacific'}; 



























return; 


