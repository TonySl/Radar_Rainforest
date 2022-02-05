%% 
%
% This script will create Fig. 1 in Tao et al. 
% All input files can be downloaded at: https://doi.org/10.6084/m9.figshare.14061428.v3

% Author: Shengli Tao
% Email: sltao@pku.edu.cn
% Date: First version in 08.2018. Formatted in 01.2022.


clear
clc


warning off
addpath('./util')
addpath('./mat')


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
blue_wenwen=[186 220 249;121 189 243;3 118 206;0 85 153]/255;  
warm_wenwen=[253 196 178;249 161 130; 240 114 69;226 65 8; 170 47 3]/255;
yellow_wenwen=[248 192 17]/255;

colors=[1 1 1; flipud(warm_wenwen);green_wenwen(1,:);green_wenwen(2,:)];  

%%
figure
set(gcf,'position',[46  46   950   600],...
                             'color','w','paperpositionmode','auto')
   

%% %%%%%% Americas continental average  
axes1=axes('position',[.08 .45 .27-0.02 .3]);

load Americas_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff %%% The final merged radar data-set. Provided as one single time series. Each row represents a month, each column a pixel.
load Americas_tropic_belt_Final_time_series_LinearCDF_Noutlier_ERSSimpleDiff %%% The final merged radar data-set. Provided per sensor.
load Americas_tropical_belt_LUCC %%% For each radar pixel, we calculated its dominant LUCC using ESA CCI landcover images.
load Americas_tropical_belt_fishnet %%% Radar pixel ID, x, y, et al...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Hansen ITR mask
% load Americas_tropic_belt_HansenCover90_non_forest_count %%% For each pixel, we calculated its non-forest ratio using Hansen et al.' dataset.
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<5 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Vancutsem ITR mask
load Am_intact_vancutsem_count %%%% For each Radar pixel, we calculated its number of undisturbed 30-m cells within it
size2=(0.25/0.000269495)^2; %%%% size of the Vancutsem pixels (~30-m)
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
%%%% index of ITR pixel
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
% disp(sum(ind_1)) %%% How many ITR pixels?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allcell_intact_ratio_selected=allcell_intact_ratio(ind_1);
% disp(sum(allcell_intact_ratio_selected>97))


Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1);  %%%% Subset to intact tropical rainforest


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Averaging across pixels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Final_time_series_FullRadar_allpixel_average=nanmean(Final_time_series_FullRadar,2);
Final_time_series_FullRadar_allpixel_average_anomaly=Final_time_series_FullRadar_allpixel_average-nanmean(Final_time_series_FullRadar_allpixel_average);


ax = gca;
ax.YColor = 'k';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add CWD (monthly water deficit) as background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load D:\Research\Drought\Tropical_belt\mat\Americas_tropical_belt_CHIRPS_19902018 %%% For each pixel, we calculated its monthly precipitation using CHIRPS data
load D:\Research\Drought\Tropical_belt\mat\Americas_tropical_belt_GLEAM_ET_19922018 %%% For each pixel, we calculated its monthly ET using GLEAM data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CHIRPS_date_final=CHIRPS_date_final(25:348); %%% start from 1992
CHIRPS_date_decimal=CHIRPS_date_decimal(25:348);
CHIRPS_data_final=CHIRPS_data_final(25:348,:);

CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);

%%%%%%%%%%%%%%%%%%%%%%%% WD of each month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Water deficit for each 0.25 degree cell.
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
WD_time_series_CHIRPS_drought_pixels=WD_time_series_CHIRPS(:,ind_1);  
%%%% named as "xxxxx_drought_pixels" but actually for all pixels. Lazy to change it....

% histogram(nanmin(WD_time_series_CHIRPS_drought_pixels,2))
min_WD_pixel=nanmin(WD_time_series_CHIRPS_drought_pixels,[],1);
% disp(-3*std(min_WD_pixel,'omitnan'))
WD_time_series_CHIRPS_drought_pixels(:,min_WD_pixel<-3*std(min_WD_pixel,'omitnan'))=NaN; %%%% outlier removal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WD_time_series_CHIRPS_drought_pixels=nanmean(WD_time_series_CHIRPS_drought_pixels,2);
% WD_time_series_CHIRPS_drought_pixels_anomaly=WD_time_series_CHIRPS_drought_pixels-nanmean(WD_time_series_CHIRPS_drought_pixels);

% wd_std=(WD_time_series_CHIRPS_drought_pixels-nanmean(WD_time_series_CHIRPS_drought_pixels))/nanstd(WD_time_series_CHIRPS_drought_pixels);
% WD_time_series_CHIRPS_drought_pixels(wd_std>-1.0)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Radar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis left
ax = gca;
ax.YColor = 'k';

plot1=plot(Final_time_monthly_decimal,Final_time_series_FullRadar_allpixel_average_anomaly,'-','color',radar_line1,'LineWidth',1); %CT_movemean_yellow(5,:)red_light_final

hold on
load DBEST_Americas_continent_Vancutsemratio95 %%%% DBEST trend calculated from R codes.  DBEST removes the seasonal component of a time series
plot2=plot(Final_time_monthly_decimal,Trend,'-o','color',radar_line2,'LineWidth',2,'markersize',1); %%%%red_dark_final


set(gca,'xtick',1995+1/12:5:2015+1/12)
% set(gca,'xticklabel','')
set(gca,'xticklabel',{'1995', '2000', '2005', '2010', '2015'})


set(gca,'ylim',[-0.32 0.32])
set(gca,'ytick',-0.1:0.1:0.1)



ylabel_pos=ylabel('Radar anomaly (dB)','FontSize',10);
set(ylabel_pos, 'Units', 'Normalized', 'Position', [-0.15, 0.42, 0]);

box on
set(gca,'xlim',[1991 2020])
set(gca,'fontsize',10)

ylim1=ylim;

ax.XAxis.MinorTick = 'on'; 


% grid on
h=gca;
h.XRuler.TickLength=[0.02 0.02];

ylim([ylim1(1),ylim1(2)+0.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis right
ax = gca;
ax.YColor = 'k';

bar_h=bar(Final_time_monthly_decimal,WD_time_series_CHIRPS_drought_pixels);
bar_h.FaceColor = CT_grey(4,:);
min_wd=min(WD_time_series_CHIRPS_drought_pixels);
set(gca,'ylim',[min_wd*(2.6) 0])
set(gca,'ytick',[round(min_wd) round(min_wd/2) 0])


%% %%%%%% Africa continental average 
%%%%% See above section (Americas continental average) for code explanations.

CT_movemean_red=cbrewer('seq', 'BuGn', 9);
radar_line1=CT_movemean_red(6,:);
radar_line2=CT_movemean_red(8,:);


axes('position',[.38 .45 .27-0.02 .3]);

load Africa_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
load Africa_tropic_belt_Final_time_series_LinearCDF_Noutlier_ERSSimpleDiff
load Africa_tropical_belt_LUCC
load Africa_tropical_belt_fishnet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Hansen ITR mask
% load Africa_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<5 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Vancutsem ITR mask
load Af_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
% disp(sum(ind_1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allcell_intact_ratio_selected=allcell_intact_ratio(ind_1);
% disp(sum(allcell_intact_ratio_selected>97))


Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Final_time_series_FullRadar_allpixel_average=nanmean(Final_time_series_FullRadar,2);
Final_time_series_FullRadar_allpixel_average_anomaly=Final_time_series_FullRadar_allpixel_average-nanmean(Final_time_series_FullRadar_allpixel_average);


ax = gca;
ax.YColor = 'k';


%%%%%%%%%%%%%%%%%%%%%% Add background of CWD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load D:\Research\Drought\Tropical_belt\mat\Africa_tropical_belt_CHIRPS_19902018
load D:\Research\Drought\Tropical_belt\mat\Africa_tropical_belt_GLEAM_ET_19922018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start from 1992
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHIRPS_date_final=CHIRPS_date_final(25:348);
CHIRPS_date_decimal=CHIRPS_date_decimal(25:348);
CHIRPS_data_final=CHIRPS_data_final(25:348,:);

CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);

%%%%%%%%%%%%%%%%%%%%%%%% WD of each month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Water deficit for each 0.25 degree cell.
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
WD_time_series_CHIRPS_drought_pixels=WD_time_series_CHIRPS(:,ind_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% histogram(nanmin(WD_time_series_CHIRPS_drought_pixels,2))
min_WD_pixel=nanmin(WD_time_series_CHIRPS_drought_pixels,[],1);
% disp(-3*std(min_WD_pixel,'omitnan'))
WD_time_series_CHIRPS_drought_pixels(:,min_WD_pixel<-3*std(min_WD_pixel,'omitnan'))=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WD_time_series_CHIRPS_drought_pixels=nanmean(WD_time_series_CHIRPS_drought_pixels,2);
% WD_time_series_CHIRPS_drought_pixels_anomaly=WD_time_series_CHIRPS_drought_pixels-nanmean(WD_time_series_CHIRPS_drought_pixels);

% wd_std=(WD_time_series_CHIRPS_drought_pixels-nanmean(WD_time_series_CHIRPS_drought_pixels))/nanstd(WD_time_series_CHIRPS_drought_pixels);
% WD_time_series_CHIRPS_drought_pixels(wd_std>-1.0)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Radar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis left
ax = gca;
ax.YColor = 'k';

  
plot(Final_time_monthly_decimal,Final_time_series_FullRadar_allpixel_average_anomaly,'-','color',radar_line1,'LineWidth',1) %CT_movemean_green(5,:)
hold on
load DBEST_Africa_continent_Vancutsemratio95
plot(Final_time_monthly_decimal,Trend,'-o','color',radar_line2,'LineWidth',2,'markersize',1)


set(gca,'ylim',[-0.32 0.32])
set(gca,'ytick',-0.1:0.1:0.1)
set(gca,'yticklabel',{}) 

set(gca,'xlim',[1991 2020])

set(gca,'xtick',[1995+1/12:5:2015+1/12 ])

set(gca,'xticklabel',{'1995', '2000', '2005', '2010', '2015'})

set(gca,'FontSize',10)


ylim1=ylim;

hold on

ax.XAxis.MinorTick = 'on'; 

% grid on
h=gca;
h.XRuler.TickLength=[0.02 0.02];

ylim([ylim1(1),ylim1(2)+0.1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis right
ax = gca;
ax.YColor = 'k';

bar_h=bar(Final_time_monthly_decimal,WD_time_series_CHIRPS_drought_pixels);
bar_h.FaceColor = CT_grey(4,:);
min_wd=min(WD_time_series_CHIRPS_drought_pixels);
set(gca,'ylim',[min_wd*(2.6) 0])
set(gca,'ytick',[round(min_wd) round(min_wd/2) 0])

%% %%%%%% Asia continental average 
%%%%% See above section (Americas continental average) for code explanations.

CT_movemean_red=cbrewer('seq', 'PuBu', 9);
radar_line1=CT_movemean_red(6,:);
radar_line2=CT_movemean_red(8,:);

axes('position',[.68 .45 .27-0.02 .3]);

load Asia_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
load Asia_tropic_belt_Final_time_series_LinearCDF_Noutlier_ERSSimpleDiff

load Asia_tropical_belt_LUCC
load Asia_tropical_belt_fishnet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hansen ITR mask
% load Asia_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<5 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal

%%%% Vancutsem ITR mask
load Asia_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
% disp(sum(ind_1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allcell_intact_ratio_selected=allcell_intact_ratio(ind_1);
% disp(sum(allcell_intact_ratio_selected>97))


Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1);  %%% subset to ITR rainforest pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Final_time_series_FullRadar_allpixel_average=nanmean(Final_time_series_FullRadar,2);
Final_time_series_FullRadar_allpixel_average_anomaly=Final_time_series_FullRadar_allpixel_average-nanmean(Final_time_series_FullRadar_allpixel_average);


ax = gca;
ax.YColor = 'k';


%%%%%%%%%%%%%%%%%%%%% Add background of CWD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load D:\Research\Drought\Tropical_belt\mat\Asia_tropical_belt_CHIRPS_19902018
load D:\Research\Drought\Tropical_belt\mat\Asia_tropical_belt_GLEAM_ET_19922018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start from 1992
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHIRPS_date_final=CHIRPS_date_final(25:348);
CHIRPS_date_decimal=CHIRPS_date_decimal(25:348);
CHIRPS_data_final=CHIRPS_data_final(25:348,:);

CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);

%%%%%%%%%%%%%%%%%%%%%%%% WD of each month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Water deficit for each 0.25 degree cell.
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
WD_time_series_CHIRPS_drought_pixels=WD_time_series_CHIRPS(:,ind_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% histogram(nanmin(WD_time_series_CHIRPS_drought_pixels,2))
min_WD_pixel=nanmin(WD_time_series_CHIRPS_drought_pixels,[],1);
% disp(-3*std(min_WD_pixel,'omitnan'))
WD_time_series_CHIRPS_drought_pixels(:,min_WD_pixel<-3*std(min_WD_pixel,'omitnan'))=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WD_time_series_CHIRPS_drought_pixels=nanmean(WD_time_series_CHIRPS_drought_pixels,2);
% WD_time_series_CHIRPS_drought_pixels_anomaly=WD_time_series_CHIRPS_drought_pixels-nanmean(WD_time_series_CHIRPS_drought_pixels);

% wd_std=(WD_time_series_CHIRPS_drought_pixels-nanmean(WD_time_series_CHIRPS_drought_pixels))/nanstd(WD_time_series_CHIRPS_drought_pixels);
% WD_time_series_CHIRPS_drought_pixels(wd_std>-1.0)=NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Radar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis left
ax = gca;
ax.YColor = 'k';

plot(Final_time_monthly_decimal,Final_time_series_FullRadar_allpixel_average_anomaly,'-','color',radar_line1,'LineWidth',1) %CT_movemean_blue
hold on
load DBEST_Asia_continent_Vancutsemratio95
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



ax.XAxis.MinorTick = 'on';  

% grid on
h=gca;
h.XRuler.TickLength=[0.02 0.02];

ylim([ylim1(1),ylim1(2)+0.1])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis right
ax = gca;
ax.YColor = 'k';

bar_h=bar(Final_time_monthly_decimal,WD_time_series_CHIRPS_drought_pixels);
bar_h.FaceColor = CT_grey(4,:);
min_wd=min(WD_time_series_CHIRPS_drought_pixels);
set(gca,'ylim',[min_wd*(2.6) 0])
set(gca,'ytick',[round(min_wd) round(min_wd/2) 0])
ylabel_pos=ylabel('{\itD} (mm)','FontSize',10');
set(ylabel_pos, 'Units', 'Normalized', 'Position', [1.11, 0.76, 0]);


%%
% print(gcf,'-dtiff','-r300',strcat('./Fig1.tif')) 
