%% 
%
% This script will create Fig. 1 in Tao et al. 

% Author: Shengli Tao
% Email: sltao1990@gmail.com
% Date: First version in 08.2018. Formatted in 02.2021.


clear
clc


warning off
addpath('./util')

%% Define colors 

CT_grey=cbrewer('seq', 'Greys', 9);
CT_Dark=cbrewer('qual', 'Dark2', 9);
CT_Set2=cbrewer('qual', 'Set2', 9);
CT_PuOr=cbrewer('div', 'BrBG', 11);

CT_movemean_red=cbrewer('seq', 'Oranges', 9);
radar_line1=CT_movemean_red(6,:);
radar_line2=CT_movemean_red(8,:);


grey_final=[201 205 206]/255;

back_alpha=0.3; %%%% transparent background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Color set by wenwen %%%%%
green_wenwen=[155 213 216;95 197 202;31 159 164; 4 105 109]/255;
blue_wenwen=[186 220 249;121 189 243;3 118 206;0 85 153]/255; % 52 150 229;
warm_wenwen=[253 196 178;249 161 130; 240 114 69;226 65 8; 170 47 3]/255;
yellow_wenwen=[248 192 17]/255;

colors=[1 1 1; flipud(warm_wenwen);green_wenwen(1,:);green_wenwen(2,:)]; % color by wenwen


%% Figure
figure
set(gcf,'position',[46  46   950   600],...
                             'color','w','paperpositionmode','auto')
   

                         
%% Americas, Continental average  

axes1=axes('position',[.08 .62 .27-0.02 .3]);

load ./mat/Americas_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff %%% The final merged radar data-set. Provided as one single time series. Each row represents a month, each column a pixel.
load ./mat/Americas_tropic_belt_Final_time_series_LinearCDF_Noutlier_ERSSimpleDiff %%% The final merged radar data-set. Provided per sensor.

load ./mat/Americas_tropic_belt_HansenCover90_non_forest_count %%% For each pixel, we calculated its non-forest ratio using Hansen et al.' dataset.
load ./mat/Americas_tropical_belt_LUCC %%% For each pixel, we calculated its dominant LUCC using ESA CCI landcover images.
load ./mat/Americas_tropical_belt_fishnet %%% Pixel ID et al...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Averaging across pixels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Final_time_series_FullRadar_allpixel_average=nanmean(Final_time_series_FullRadar,2);
Final_time_series_FullRadar_allpixel_average_anomaly=Final_time_series_FullRadar_allpixel_average-nanmean(Final_time_series_FullRadar_allpixel_average);

% save('Americas_tropical_belt_Fig1_forBFAST','Final_time_series_FullRadar_allpixel_average_anomaly','-v7.3') %%%% same data for BFAST break point detection in R



ax = gca;
ax.YColor = 'k';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add background of CWD (monthly water deficit) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ./mat/Americas_tropical_belt_CHIRPS_19902018 %%% For each pixel, we calculated its monthly precipitation using CHIRPS data
load ./mat/Americas_tropical_belt_GLEAM_ET_19922018 %%% For each pixel, we calculated its monthly ET using GLEAM data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHIRPS_date_final=CHIRPS_date_final(25:348); % Start from 1992
CHIRPS_date_decimal=CHIRPS_date_decimal(25:348);
CHIRPS_data_final=CHIRPS_data_final(25:348,:);

CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);

%%%%%%%%%%%%%%%%%%%%%%%% WD of each month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Cummulative water deficit for each 0.25 degree cell.
CHIRPS_p_evap_allmonth=CHIRPS_data_final'-GLEAM_ET_data;       
p_evap=CHIRPS_p_evap_allmonth(:,1);       
p_evap(p_evap>0)=0;               

CHIRPS_WD_allmonth=[];
CHIRPS_WD_allmonth(:,1)=p_evap;
filecount=size(CHIRPS_p_evap_allmonth,2);
for i=2:length(CHIRPS_date_final)

        p_evap=CHIRPS_p_evap_allmonth(:,i);       

        Wn1=CHIRPS_WD_allmonth(:,i-1);  

        Wn=Wn1+p_evap;

        Wn(Wn>0)=0;  
        CHIRPS_WD_allmonth(:,i)=Wn;

end

        
%%%%%%%%%%%%%%%%%%%%%%%%% WD of every month from 1992.10 %%%%%%%%%%%%%%%%%%%%%%%%%%%  
temp_time=num2str(ERS2_date_final(1));
start_year=str2num(temp_time(1:4));
start_mon=str2num(temp_time(5:6));
temp_time=num2str(A_date_final(end));
end_year=str2num(temp_time(1:4));
end_mon=str2num(temp_time(5:6));

%%%%% create a full time series of WD from start date to end data at monthly time-step
WD_time_monthly=[];
WD_time_monthly_decimal=[];

for i=start_year:end_year
    if i==start_year
        For_loop_start_mon=start_mon;
        For_loop_end_mon=12;
    elseif i==end_year
        For_loop_start_mon=1;
        For_loop_end_mon=end_mon;
    else
        For_loop_start_mon=1;
        For_loop_end_mon=12;
    end


    for j=For_loop_start_mon:For_loop_end_mon

        temp_time=i*100+j;
        WD_time_monthly=[WD_time_monthly;temp_time];
        temp_time2=num2str(temp_time);
        temp_time_decimal=str2num(temp_time2(1:4))+str2num(temp_time2(5:6))/12;
        WD_time_monthly_decimal=[WD_time_monthly_decimal;temp_time_decimal];

    end

end


[c,d]=ismember(CHIRPS_date_final,WD_time_monthly);
WD_time_series_CHIRPS=nan(length(WD_time_monthly),size(A_data_final,2));
WD_time_series_CHIRPS(d(c),:)=CHIRPS_WD_allmonth(:,c)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WD_time_series_CHIRPS_drought_pixels=WD_time_series_CHIRPS(:,ind_1); 
%%%% named as "xxxxx_drought_pixels" but actually for all pixels. Lazy to change it....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WD_time_series_CHIRPS_drought_pixels=nanmean(WD_time_series_CHIRPS_drought_pixels,2);

%%%%%%%%%%%%%%%%%%%%%%%% bar figure as background, individual color %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    
%     min(WD_time_series_CHIRPS_drought_pixels)
    TRMM_bin=[-185:5:0];  

    CT_back=cbrewer('seq', 'Greys', length(TRMM_bin));
    CT_back=flipud(CT_back);

    colorCode=nan(length(WD_time_monthly),1);

    for irain = 1:length(TRMM_bin)-1
        ok = WD_time_series_CHIRPS_drought_pixels>=(TRMM_bin(irain)) & WD_time_series_CHIRPS_drought_pixels<=(TRMM_bin(irain+1));
        colorCode(ok) = irain;
    end
    

    for i=1:length(WD_time_monthly_decimal)
        bar(WD_time_monthly_decimal(i),10,1/12,'facecolor',CT_back(colorCode(i),:),'edgecolor','none','facealpha',back_alpha)
        hold on    
        bar(WD_time_monthly_decimal(i),-10,1/12,'facecolor',CT_back(colorCode(i),:),'edgecolor','none','facealpha',back_alpha)
        hold on
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Radar data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot1=plot(Final_time_monthly_decimal,Final_time_series_FullRadar_allpixel_average_anomaly,'-','color',radar_line1,'LineWidth',1); 

hold on

load ./mat/DBEST_Americas_Fig1_DBEST_Trend_Tropical_belt %%%% This is the DBEST generalized time-series obtained in R. Averaged values across all intact tropical rainforest pixels.
%%% DBEST removes the seasonal component of a time series

plot2=plot(Final_time_monthly_decimal,Trend,'-o','color',radar_line2,'LineWidth',2,'markersize',1); 


set(gca,'xtick',[1995+1/12:5:2015+1/12])
set(gca,'xticklabel',{'1995', '2000', '2005', '2010', '2015'})


set(gca,'ylim',[-0.3 0.3])
set(gca,'ytick',-0.1:0.1:0.1)



ylabel_pos=ylabel('Radar anomaly (dB)','FontSize',10); %,'FontWeight','bold 


box on

set(gca,'xlim',[1991 2020])
set(gca,'fontsize',10)

ylim1=ylim;
xlim1=xlim;
hold on


ax.XAxis.MinorTick = 'on'; 

h=gca;
h.XRuler.TickLength=[0.02 0.02];

%% Africa, Continental average.

CT_movemean_red=cbrewer('seq', 'BuGn', 9);
radar_line1=CT_movemean_red(6,:);
radar_line2=CT_movemean_red(8,:);


axes('position',[.36 .62 .27-0.02 .3]);

load ./mat/Africa_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff  %%%%% See above section (Americas, Continental average) for explanations. Same for codes below.
load ./mat/Africa_tropic_belt_Final_time_series_LinearCDF_Noutlier_ERSSimpleDiff

load ./mat/Africa_tropic_belt_HansenCover90_non_forest_count
load ./mat/Africa_tropical_belt_LUCC
load ./mat/Africa_tropical_belt_fishnet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Averaging across pixels     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Final_time_series_FullRadar_allpixel_average=nanmean(Final_time_series_FullRadar,2);
Final_time_series_FullRadar_allpixel_average_anomaly=Final_time_series_FullRadar_allpixel_average-nanmean(Final_time_series_FullRadar_allpixel_average);

% save('Africa_tropical_belt_Fig1_forBFAST','Final_time_series_FullRadar_allpixel_average_anomaly','-v7.3')

ax = gca;
ax.YColor = 'k';


%%%%%%%%%%%%%%%%%%%%%% Add background of CWD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ./mat/Africa_tropical_belt_CHIRPS_19902018
load ./mat/Africa_tropical_belt_GLEAM_ET_19922018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start from 1992
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHIRPS_date_final=CHIRPS_date_final(25:348);
CHIRPS_date_decimal=CHIRPS_date_decimal(25:348);
CHIRPS_data_final=CHIRPS_data_final(25:348,:);

CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);

%%%%%%%%%%%%%%%%%%%%%%%% WD of each month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Cummulative water deficit for each 0.25 degree cell.
CHIRPS_p_evap_allmonth=CHIRPS_data_final'-GLEAM_ET_data;       
p_evap=CHIRPS_p_evap_allmonth(:,1);       
p_evap(p_evap>0)=0;               

CHIRPS_WD_allmonth=[];
CHIRPS_WD_allmonth(:,1)=p_evap;
filecount=size(CHIRPS_p_evap_allmonth,2);
for i=2:length(CHIRPS_date_final)

        p_evap=CHIRPS_p_evap_allmonth(:,i);       

        Wn1=CHIRPS_WD_allmonth(:,i-1);

        Wn=Wn1+p_evap;

        Wn(Wn>0)=0;  
        CHIRPS_WD_allmonth(:,i)=Wn;

end
        
%%%%%%%%%%%%%%%%%%%%%%%%% WD of every month from 1992.10 %%%%%%%%%%%%%%%%%%%%%%%%%%%  
temp_time=num2str(ERS2_date_final(1));
start_year=str2num(temp_time(1:4));
start_mon=str2num(temp_time(5:6));
temp_time=num2str(A_date_final(end));
end_year=str2num(temp_time(1:4));
end_mon=str2num(temp_time(5:6));

%%%%% create a full time series of WD from start date to end data at monthly time-step
WD_time_monthly=[];
WD_time_monthly_decimal=[];

for i=start_year:end_year
    if i==start_year
        For_loop_start_mon=start_mon;
        For_loop_end_mon=12;
    elseif i==end_year
        For_loop_start_mon=1;
        For_loop_end_mon=end_mon;
    else
        For_loop_start_mon=1;
        For_loop_end_mon=12;
    end
    
    for j=For_loop_start_mon:For_loop_end_mon

        temp_time=i*100+j;
        WD_time_monthly=[WD_time_monthly;temp_time];
        temp_time2=num2str(temp_time);
        temp_time_decimal=str2num(temp_time2(1:4))+str2num(temp_time2(5:6))/12;
        WD_time_monthly_decimal=[WD_time_monthly_decimal;temp_time_decimal];

    end

end


[c,d]=ismember(CHIRPS_date_final,WD_time_monthly);
WD_time_series_CHIRPS=nan(length(WD_time_monthly),size(A_data_final,2));
WD_time_series_CHIRPS(d(c),:)=CHIRPS_WD_allmonth(:,c)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WD_time_series_CHIRPS_drought_pixels=WD_time_series_CHIRPS(:,ind_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WD_time_series_CHIRPS_drought_pixels=nanmean(WD_time_series_CHIRPS_drought_pixels,2);

%%%%%%%%%%%%%%%%%%%%%%%% bar figure as background, individual color %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    
%     min(WD_time_series_CHIRPS_drought_pixels)
    TRMM_bin=[-130:5:0];  
    TRMM_bin=[TRMM_bin 0];
    CT_back=cbrewer('seq', 'Greys', length(TRMM_bin));  
    CT_back=flipud(CT_back);

    colorCode=nan(length(WD_time_monthly),1);

    for irain = 1:length(TRMM_bin)-1
        ok = WD_time_series_CHIRPS_drought_pixels>=(TRMM_bin(irain)) & WD_time_series_CHIRPS_drought_pixels<=(TRMM_bin(irain+1));
        colorCode(ok) = irain;
    end
    

    for i=1:length(WD_time_monthly_decimal)
        bar(WD_time_monthly_decimal(i),10,1/12,'facecolor',CT_back(colorCode(i),:),'edgecolor','none','facealpha',back_alpha)
        hold on    
        bar(WD_time_monthly_decimal(i),-10,1/12,'facecolor',CT_back(colorCode(i),:),'edgecolor','none','facealpha',back_alpha)
        hold on
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Radar data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
plot(Final_time_monthly_decimal,Final_time_series_FullRadar_allpixel_average_anomaly,'-','color',radar_line1,'LineWidth',1) %CT_movemean_green(5,:)
hold on
load ./mat/DBEST_Africa_Fig1_DBEST_Trend_Tropical_belt
plot(Final_time_monthly_decimal,Trend,'-o','color',radar_line2,'LineWidth',2,'markersize',1)


set(gca,'ylim',[-0.33 0.33])
set(gca,'ytick',-0.1:0.1:0.1)
set(gca,'yticklabel',{})  

set(gca,'xlim',[1991 2020])

set(gca,'xtick',[1995+1/12:5:2015+1/12 ])
set(gca,'xticklabel',{'1995', '2000', '2005', '2010', '2015'})

set(gca,'FontSize',10)
box on

ylim1=ylim;
xlim1=xlim;
hold on


ax.XAxis.MinorTick = 'on';

h=gca;
h.XRuler.TickLength=[0.02 0.02];

%% Asia, Continental average 


CT_movemean_red=cbrewer('seq', 'PuBu', 9);
radar_line1=CT_movemean_red(6,:);
radar_line2=CT_movemean_red(8,:);

% axes('position',[.07 .03 .4 .32])
axes('position',[.64 .62 .27-0.02 .3]);

load ./mat/Asia_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff %%%%% See above section (Americas, Continental average) for explanations. Same for codes below.
load ./mat/Asia_tropic_belt_Final_time_series_LinearCDF_Noutlier_ERSSimpleDiff

load ./mat/Asia_tropic_belt_HansenCover90_non_forest_count
load ./mat/Asia_tropical_belt_LUCC
load ./mat/Asia_tropical_belt_fishnet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Averaging across pixels   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Final_time_series_FullRadar_allpixel_average=nanmean(Final_time_series_FullRadar,2);
Final_time_series_FullRadar_allpixel_average_anomaly=Final_time_series_FullRadar_allpixel_average-nanmean(Final_time_series_FullRadar_allpixel_average);

% save('Asia_tropical_belt_Fig1_forBFAST','Final_time_series_FullRadar_allpixel_average_anomaly','-v7.3')

ax = gca;
ax.YColor = 'k';


%%%%%%%%%%%%%%%%%%%%% Add background of CWD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ./mat/Asia_tropical_belt_CHIRPS_19902018
load ./mat/Asia_tropical_belt_GLEAM_ET_19922018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start from 1992
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHIRPS_date_final=CHIRPS_date_final(25:348);
CHIRPS_date_decimal=CHIRPS_date_decimal(25:348);
CHIRPS_data_final=CHIRPS_data_final(25:348,:);

CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);

%%%%%%%%%%%%%%%%%%%%%%%% WD of each month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Cummulative water deficit for each 0.25 degree cell.
CHIRPS_p_evap_allmonth=CHIRPS_data_final'-GLEAM_ET_data;       
p_evap=CHIRPS_p_evap_allmonth(:,1);       
p_evap(p_evap>0)=0;               

CHIRPS_WD_allmonth=[];
CHIRPS_WD_allmonth(:,1)=p_evap;
filecount=size(CHIRPS_p_evap_allmonth,2);
for i=2:length(CHIRPS_date_final)

        p_evap=CHIRPS_p_evap_allmonth(:,i);       

        Wn1=CHIRPS_WD_allmonth(:,i-1); 

        Wn=Wn1+p_evap;

        Wn(Wn>0)=0;  
        CHIRPS_WD_allmonth(:,i)=Wn;

end
    
%%%%%%%%%%%%%%%%%%%%%%%%% WD of every month from 1992.10 %%%%%%%%%%%%%%%%%%%%%%%%%%%  
temp_time=num2str(ERS2_date_final(1));
start_year=str2num(temp_time(1:4));
start_mon=str2num(temp_time(5:6));
temp_time=num2str(A_date_final(end));
end_year=str2num(temp_time(1:4));
end_mon=str2num(temp_time(5:6));

%%%%% create a full time series of WD from start date to end data at monthly time-step
WD_time_monthly=[];
WD_time_monthly_decimal=[];

for i=start_year:end_year
    if i==start_year
        For_loop_start_mon=start_mon;
        For_loop_end_mon=12;
    elseif i==end_year
        For_loop_start_mon=1;
        For_loop_end_mon=end_mon;
    else
        For_loop_start_mon=1;
        For_loop_end_mon=12;
    end
    
    for j=For_loop_start_mon:For_loop_end_mon

        temp_time=i*100+j;
        WD_time_monthly=[WD_time_monthly;temp_time];
        temp_time2=num2str(temp_time);
        temp_time_decimal=str2num(temp_time2(1:4))+str2num(temp_time2(5:6))/12;
        WD_time_monthly_decimal=[WD_time_monthly_decimal;temp_time_decimal];

    end

end


[c,d]=ismember(CHIRPS_date_final,WD_time_monthly);
WD_time_series_CHIRPS=nan(length(WD_time_monthly),size(A_data_final,2));
WD_time_series_CHIRPS(d(c),:)=CHIRPS_WD_allmonth(:,c)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WD_time_series_CHIRPS_drought_pixels=WD_time_series_CHIRPS(:,ind_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


WD_time_series_CHIRPS_drought_pixels=nanmean(WD_time_series_CHIRPS_drought_pixels,2);

%%%%%%%%%%%%%%%%%%%%%%%% bar figure as background, individual color %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    
%     min(WD_time_series_CHIRPS_drought_pixels)
    TRMM_bin=[-63:2:0];  

    
    CT_back=cbrewer('seq', 'Greys', length(TRMM_bin));
    CT_back=flipud(CT_back);

    colorCode=nan(length(WD_time_monthly),1);

    for irain = 1:length(TRMM_bin)-1
        ok = WD_time_series_CHIRPS_drought_pixels>=(TRMM_bin(irain)) & WD_time_series_CHIRPS_drought_pixels<=(TRMM_bin(irain+1));
        colorCode(ok) = irain;
    end


    for i=1:length(WD_time_monthly_decimal)
        bar(WD_time_monthly_decimal(i),10,1/12,'facecolor',CT_back(colorCode(i),:),'edgecolor','none','facealpha',back_alpha)
        hold on    
        bar(WD_time_monthly_decimal(i),-10,1/12,'facecolor',CT_back(colorCode(i),:),'edgecolor','none','facealpha',back_alpha)
        hold on
    end


end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Radar data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(Final_time_monthly_decimal,Final_time_series_FullRadar_allpixel_average_anomaly,'-','color',radar_line1,'LineWidth',1)  
hold on
load ./mat/DBEST_Asia_Fig1_DBEST_Trend_Tropical_belt
plot(Final_time_monthly_decimal,Trend,'-o','color',radar_line2,'LineWidth',2,'markersize',1)


set(gca,'xtick',[1995+1/12:5:2015+1/12])
set(gca,'xticklabel',{'1995', '2000', '2005', '2010', '2015'})


set(gca,'FontSize',10)

box on

set(gca,'xlim',[1991 2020])
set(gca,'ylim',[-0.32 0.32])
set(gca,'ytick',-0.1:0.1:0.1)
set(gca,'yticklabel',{})  

ylim1=ylim;
xlim1=xlim;
hold on



ax.XAxis.MinorTick = 'on'; 
h=gca;
h.XRuler.TickLength=[0.02 0.02];

%% Legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Legend on the right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('position',[.91 .62 .010 .3]);
yyaxis right
ax = gca;
ax.YColor = 'k';

B = imagesc(meshgrid(1:length(TRMM_bin),0:5)');
alpha(B,.3)
set(gca,'xtick',[])
set(gca,'ytick',[])

colormap(CT_back)

yyaxis left
ax = gca;
ax.YColor = 'k';
set(gca,'xtick',[])
set(gca,'ytick',[])




%% Spatial pattern of Radar trend, pixel level

pos2 = [0.05 0.3 0.87 0.3];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

%% Americas, Pixel trend

load ./mat/Americas_tropical_belt_fishnet  %%%%% See above section (Americas, Continental average) for explanations. Same for codes below.
load ./mat/Americas_tropic_belt_HansenCover90_non_forest_count
load ./mat/Americas_tropical_belt_LUCC

xll=min(Am_cell_centers(Am_cell_area>0.062,1));
yll=min(Am_cell_centers(Am_cell_area>0.062,2));

xright=max(Am_cell_centers(Am_cell_area>0.062,1));
yup=max(Am_cell_centers(Am_cell_area>0.062,2));

nrows=int32((yup-yll)/0.25+1);
ncols=int32((xright-xll)/0.25+1);
cellsize=0.25;

R_map=[1/0.25 yup xll];


%%%% This is the DBEST generalized trend calculated in R. Trend for each intact tropical rainforest pixel is provided.
load ./mat/DBEST_Americas_Trend_Tropical_belt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal

all_scores=all_scores(ind_1);
Am_cell_centers=Am_cell_centers(ind_1,:);

% 100*sum(all_scores<0)/length(all_scores)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% catogorize the trends %%%%%%%%%%%%%%%%%%%%%%
P=nan(nrows,ncols);
[row, col] = setpostn(P, R_map, Am_cell_centers(:,2), Am_cell_centers(:,1));
row=double(nrows)-double(row);

P(sub2ind(size(P),row,col))=all_scores;

P_new=nan(size(P));
P_new(isnan(P))=0;
P_new(P<-1.3*0.001)=1;
P_new(P>=-1.3*0.001 & P<-1*0.001)=2;
P_new(P>=-1*0.001 & P<-0.6*0.001)=3;
P_new(P>=-0.6*0.001 & P<-0.3*0.001)=4;
P_new(P>=-0.3*0.001 & P<0)=5;
P_new(P>=0 & P<0.3*0.001)=6;
P_new(P>=0.3*0.001)=7;


% geotiffwrite('pixel_trend_radar.tif',flipud(P_new),R_map)
% geotiffwrite('pixel_trend_radar_original.tif',flipud(P),R_map)
 
axesm('MapProjection','mercator','maplatlimit',[-19 15],'maplonlimit',[-100 172],...  
'ParallelLabel','off','PlabelMeridian','west','MeridianLabel','off','MLabelParallel','south',...
'FontSize',6,'FontWeight','normal','PLineLocation',20,'MLineLocation',20);

setm(gca,'Frame','on')
setm(gca,'FEdgeColor',[0 0 0])
setm(gca,'FLineWidth',0.5)

latitudes = yll:cellsize:yup;  
latitudes=flipud(latitudes');  %%%%%%%%%%% remember to flipud the latitudes %%%%%%%%%%%%%%%%%%%%%%
longitudes = xll:cellsize:xright;
[latGrid, lonGrid] = meshgrat(latitudes,longitudes);
geoshow(latGrid,lonGrid,double(P_new),'DisplayType','texturemap');  %%%% texturemap

caxis([0 8]);

colormap(gca,colors);

%% Legend

chandle = colorbar('Location','SouthOutside','FontSize',6,'FontWeight','normal'); % This line places the colorbar
set(get(chandle,'ylabel'),'String','Trend in radar signal','FontSize',22,'FontWeight','normal'); % Set the colorbar's label

set(chandle, 'ylim', [1 8])

lab = [{};{};{};{};{};{};{}];
set(chandle,'YTick',[1.5:1:7.5],'yticklabel',lab,'fontsize',8);
set(chandle, 'TickLength', [0 0]);


x=get(chandle,'Position');
x(1)=0.56;
x(2)=0.393;
x(3)=0.17;
x(4)=.016;
set(chandle,'Position',x)

load('./mat/world_country_new','world_country_new')
S = geoshape(world_country_new);
hold on
geoshow(S, 'DefaultFaceColor', [1 1 1],'DefaultEdgeColor', 'k','Linewidth',0.5,'FaceAlpha',0)
hold on
tightmap

%% Africa, Pixel trend

load ./mat/Africa_tropical_belt_fishnet %%%%% See above section (Americas, Continental average) for explanations. Same for codes below.
load ./mat/Africa_tropic_belt_HansenCover90_non_forest_count
load ./mat/Africa_tropical_belt_LUCC

xll=min(Af_cell_centers(Af_cell_area>0.062,1));
yll=min(Af_cell_centers(Af_cell_area>0.062,2));

xright=max(Af_cell_centers(Af_cell_area>0.062,1));
yup=max(Af_cell_centers(Af_cell_area>0.062,2));

nrows=int32((yup-yll)/0.25+1);
ncols=int32((xright-xll)/0.25+1);
cellsize=0.25;

R_map=[1/0.25 yup xll];


load ./mat/DBEST_Africa_Trend_Tropical_belt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal

all_scores=all_scores(ind_1);
Af_cell_centers=Af_cell_centers(ind_1,:);

% 100*sum(all_scores<0)/length(all_scores)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%% catogorize the trends %%%%%%%%%%%%%%%%%%%%%%
P=nan(nrows,ncols);
[row, col] = setpostn(P, R_map, Af_cell_centers(:,2), Af_cell_centers(:,1));
row=double(nrows)-double(row);

P(sub2ind(size(P),row,col))=all_scores;


P_new=nan(size(P));
P_new(isnan(P))=0;

P_new(P<-1.3*0.001)=1;
P_new(P>=-1.3*0.001 & P<-1*0.001)=2;
P_new(P>=-1*0.001 & P<-0.6*0.001)=3;
P_new(P>=-0.6*0.001 & P<-0.3*0.001)=4;
P_new(P>=-0.3*0.001 & P<0)=5;
P_new(P>=0 & P<0.3*0.001)=6;
P_new(P>=0.3*0.001)=7;


latitudes = yll:cellsize:yup;  
latitudes=flipud(latitudes');  
longitudes = xll:cellsize:xright;
[latGrid, lonGrid] = meshgrat(latitudes,longitudes);
geoshow(latGrid,lonGrid,double(P_new),'DisplayType','texturemap');

caxis([0 8]);  

colormap(gca,colors);


load('./mat/world_country_new','world_country_new')
S = geoshape(world_country_new);
hold on
geoshow(S, 'DefaultFaceColor', [1 1 1],'DefaultEdgeColor', 'k','Linewidth',0.5,'FaceAlpha',0)
hold on
tightmap

%% Asia, Pixel trend

load ./mat/Asia_tropical_belt_fishnet  %%%%% See above section (Americas, Continental average) for explanations. Same for codes below.
load ./mat/Asia_tropic_belt_HansenCover90_non_forest_count
load ./mat/Asia_tropical_belt_LUCC


xll=min(Asia_cell_centers(Asia_cell_area>0.062,1));
yll=min(Asia_cell_centers(Asia_cell_area>0.062,2));

xright=max(Asia_cell_centers(Asia_cell_area>0.062,1));
yup=max(Asia_cell_centers(Asia_cell_area>0.062,2));

nrows=int32((yup-yll)/0.25+1);
ncols=int32((xright-xll)/0.25+1);
cellsize=0.25;

R_map=[1/0.25 yup xll];


load ./mat/DBEST_Asia_Trend_Tropical_belt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal

all_scores=all_scores(ind_1);
Asia_cell_centers=Asia_cell_centers(ind_1,:);

% 100*sum(all_scores<0)/length(all_scores)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% catogorize the trends %%%%%%%%%%%%%%%%%%%%%%
P=nan(nrows,ncols);
[row, col] = setpostn(P, R_map, Asia_cell_centers(:,2), Asia_cell_centers(:,1));
row=double(nrows)-double(row);

P(sub2ind(size(P),row,col))=all_scores;


P_new=nan(size(P));
P_new(isnan(P))=0;

P_new(P<-1.3*0.001)=1;
P_new(P>=-1.3*0.001 & P<-1*0.001)=2;
P_new(P>=-1*0.001 & P<-0.6*0.001)=3;
P_new(P>=-0.6*0.001 & P<-0.3*0.001)=4;
P_new(P>=-0.3*0.001 & P<0)=5;
P_new(P>=0 & P<0.3*0.001)=6;
P_new(P>=0.3*0.001)=7;


latitudes = yll:cellsize:yup;  
latitudes=flipud(latitudes');  
longitudes = xll:cellsize:xright;
[latGrid, lonGrid] = meshgrat(latitudes,longitudes);
geoshow(latGrid,lonGrid,double(P_new),'DisplayType','texturemap');

caxis([0 8]); %same as the max value

colormap(gca,colors);


load('./mat/world_country_new','world_country_new')
S = geoshape(world_country_new);
hold on
geoshow(S, 'DefaultFaceColor', [1 1 1],'DefaultEdgeColor', 'k','Linewidth',0.5,'FaceAlpha',0)
hold off 
tightmap
box off
axis off




%% Print the figure
% print(gcf,'-dtiff','-r300',strcat('./Figure1_submit_2020_rgbcolor.tif'))
