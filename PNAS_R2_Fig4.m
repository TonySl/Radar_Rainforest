%% 
% This script will create Fig. 4 in Tao et al. 


% Author: Shengli Tao, assisstant professor in Peking University
% Email: sltao@pku.edu.cn
% Date: First version in 08.2018. Formatted in May.2022.


clc
clear

addpath('./util')
addpath('./mat')
warning off


legacy_length=12*2;   %%% two-year drought legecy effect.
mcwdstd_thre=-1.0; %%% z-score threshold
detrend_or_not=0; %%% whether detrend the signal or not (if the trend itself is caused by droughts, then probably no need to detrend). 
ET100mm_or_not=0; %%% GLEAM monthly ET or 100mm ET detrend the signal or not. 



figure
set(gcf,'position',[46  46   950   550],...
                             'color','w','paperpositionmode','auto')
                         
scatter_size=30;
scatter_alpha=0.9;
line_width=1;


%%  %%%%% color set by wenwen %%%%%
green_wenwen=[155 213 216;95 197 202;31 159 164; 4 105 109]/255;
blue_wenwen=[186 220 249;121 189 243;3 118 206;0 85 153]/255; % 52 150 229;
warm_wenwen=[253 196 178;249 161 130; 240 114 69;226 65 8; 170 47 3]/255;
yellow_wenwen=[248 192 17]/255;   


%% Americas 
    
load Americas_tropical_belt_fishnet %%% Pixel ID et al...
load Americas_tropical_belt_LUCC %%% For each pixel, we calculated its dominant LUCC using ESA CCI landcover images.
load Americas_tropical_belt_GLEAM_ET_19922018  %%% For each pixel, we calculated its monthly ET using GLEAM data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Hansen ITR mask 
% load Americas_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<5 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal


%%%%% Vancutsem ITR mask
load Am_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
disp(sum(ind_1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load Americas_tropical_belt_CHIRPS_19902018 %%% For each pixel, we calculated its monthly precipitation using CHIRPS data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start from 1992
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHIRPS_date_final=CHIRPS_date_final(25:348);
CHIRPS_date_decimal=CHIRPS_date_decimal(25:348);
CHIRPS_data_final=CHIRPS_data_final(25:348,:);

% CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
% CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
% CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);

%%%%%%%%%%%%%%%%%%%%%%%% WD of each month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Water deficit for each 0.25 degree cell.
if ET100mm_or_not
    CHIRPS_p_evap_allmonth=CHIRPS_data_final'-100;    
else
    CHIRPS_p_evap_allmonth=CHIRPS_data_final'-GLEAM_ET_data;    
end

p_evap=CHIRPS_p_evap_allmonth(:,1);       
p_evap(p_evap>0)=0;               

CHIRPS_WD_allmonth=[];
CHIRPS_WD_allmonth(:,1)=p_evap;
% filecount=size(CHIRPS_p_evap_allmonth,2);

for i=2:length(CHIRPS_date_final)

        p_evap=CHIRPS_p_evap_allmonth(:,i);       
        Wn1=CHIRPS_WD_allmonth(:,i-1); 
        Wn=Wn1+p_evap;
        Wn(Wn>0)=0;  
        CHIRPS_WD_allmonth(:,i)=Wn;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% maximum WD for each year %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);


CHIRPS_MWD_year=unique(CHIRPS_hydro_year_index);
CHIRPS_MWD_allpixel=nan(size(CHIRPS_WD_allmonth,1),length(CHIRPS_MWD_year));
CHIRPS_MWD_allpixel_minmonth=nan(size(CHIRPS_WD_allmonth,1),length(CHIRPS_MWD_year));

for i=1:length(CHIRPS_MWD_year)
%     disp(TRMM_MWD_year(i));
    one_year_allpixel_WD=CHIRPS_WD_allmonth(:,CHIRPS_hydro_year_index==CHIRPS_MWD_year(i));
    [one_year_allpixel_MWD,min_month]=nanmin(one_year_allpixel_WD,[],2);    % min WD of 12 months
    CHIRPS_MWD_allpixel(:,i)= one_year_allpixel_MWD;
    CHIRPS_MWD_allpixel_minmonth(:,i)= min_month;
end



drought_years=1992:2018;  %%%% potential drought years

zscore_mcwd_all=nan(size(CHIRPS_WD_allmonth,1),length(CHIRPS_MWD_year));
for i_droughtyear=1:length(drought_years)
    
    drought_year=drought_years(i_droughtyear);   
    CHIRPS_MWD_droughtyear=CHIRPS_MWD_allpixel(:,CHIRPS_MWD_year==drought_year);     
    CHIRPS_MWD_background=CHIRPS_MWD_allpixel(:,CHIRPS_MWD_year~=drought_year);
    Background_month_std=nanstd(CHIRPS_MWD_background,0,2);
    Background_month_mean=nanmean(CHIRPS_MWD_background,2);
    zscore_mcwd_all(:,i_droughtyear)=(CHIRPS_MWD_droughtyear-Background_month_mean)./Background_month_std;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Drought resistance and resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

load Americas_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1); %%%%%%%%%%%% subset to ITR


zscore_mcwd_all=zscore_mcwd_all(ind_1,:); %%%%%%%%%%%% subset to ITR
CHIRPS_MWD_allpixel_minmonth=CHIRPS_MWD_allpixel_minmonth(ind_1,:); %%%%%%%%%%%% subset to ITR



all_pixel_all_drought_severity=nan(size(Final_time_series_FullRadar,2),length(1992:2018));                          
all_pixel_all_drought_radar_resistance=nan(size(Final_time_series_FullRadar,2),length(1992:2018));
all_pixel_all_drought_radar_resilience=nan(size(Final_time_series_FullRadar,2),length(1992:2018));


for i=1:1:size(zscore_mcwd_all,1)

    ecoregion_fullRadar=Final_time_series_FullRadar(:,i);

    if detrend_or_not

        %%%%%%%%%%%%%%% Detrend won't take NaN %%%%%%%%%%%%%%
        ecoregion_fullRadar_detrend=detrend_nan(ecoregion_fullRadar);
        ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar_detrend,1,100);

    else

        ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar,1,100);

    end



    each_drought_positin = CHIRPS_MWD_allpixel_minmonth(i,:);


    this_pixel_all_drought_severity=[];
    this_pixel_all_drought_resistance=[];
    this_pixel_all_drought_resilience=[];
    
    this_pixel_all_drought_time=[];


    localmax_01ind_radar=islocalmax(ecoregion_fullRadar);

    for iiidrought=1:length(1992:2018)


        drought_center_time_1=(1991+iiidrought)*100+each_drought_positin(iiidrought);
        drought_center_time=(1991+iiidrought)+each_drought_positin(iiidrought)/12.0;

        drought_index_in_rain=find(CHIRPS_date_final==drought_center_time_1);


        if floor(drought_center_time)>=2019 %%%% our radar dataset ends in 2018.xx
            continue
        end


        one_drought_severity=zscore_mcwd_all(i,iiidrought);  

        if one_drought_severity>mcwdstd_thre %%%%% a threshold on z-score of mcwd
            continue
        end

        one_drought_severity=abs(one_drought_severity);


        this_pixel_all_drought_severity=[this_pixel_all_drought_severity;one_drought_severity];
        this_pixel_all_drought_time=[this_pixel_all_drought_time;floor(drought_center_time)];
        
        findmin_length=3; %%% Minimum radar value with +- 3 months of central drought time

        if min(drought_index_in_rain)-findmin_length<=0
            start_index=1;
        else
            start_index=min(drought_index_in_rain)-findmin_length;
        end

        if max(drought_index_in_rain)+findmin_length>length(CHIRPS_date_final)
            end_index=length(CHIRPS_date_final);
        else
            end_index=max(drought_index_in_rain)+findmin_length;
        end

        drought_time_in_rain=CHIRPS_date_final(start_index:end_index);
        drought_time_in_rain_decimal=CHIRPS_date_decimal(start_index:end_index);

        if sum(ismember(Final_time_monthly,drought_time_in_rain))==0

            Lowest_radar_value_in_drought=NaN;
            resistance_one_drought=NaN;
            resilience_one_drought=NaN;
            recover_time=NaN;

        else
            [Lowest_radar_value_in_drought,indmin]=nanmin(ecoregion_fullRadar_1_100_normalize(ismember(Final_time_monthly,drought_time_in_rain)));           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% resistance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lowest_ind_in_radar_time_series=Final_time_monthly==drought_time_in_rain(indmin);

            if find(lowest_ind_in_radar_time_series==1)<=legacy_length
                lowest_Before6mon_ind_in_rain_time_series=1:find(lowest_ind_in_radar_time_series==1); 
            else
                lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1)-legacy_length:find(lowest_ind_in_radar_time_series==1); 
            end

            temp_ind2=zeros(length(localmax_01ind_radar),1);
            temp_ind2(lowest_Before6mon_ind_in_rain_time_series)=1;
            final_Before6mon_ind_01=localmax_01ind_radar & temp_ind2;


            [before_drought_max,ind22]=(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));

            local_max_ind_before_drought=find(final_Before6mon_ind_01==1);
            time_to_drought=find(lowest_ind_in_radar_time_series==1)-local_max_ind_before_drought(ind22);

            if isempty(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)))
                recover_time=NaN;
            else
                temp2233=ecoregion_fullRadar_1_100_normalize>nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01));
                ind_recovertime=find(temp2233==1)-find(lowest_ind_in_radar_time_series==1);
                recover_time=min(ind_recovertime(ind_recovertime>0));

                if isempty(recover_time)
                    recover_time=NaN;
                end

            end

            resistance_one_drought_temp=((Lowest_radar_value_in_drought)-(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01))))/(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));
            if isempty(resistance_one_drought_temp)
                resistance_one_drought=NaN;
            else
                resistance_one_drought=resistance_one_drought_temp;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(Final_time_monthly)-find(lowest_ind_in_radar_time_series==1)<=legacy_length
                lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):length(Final_time_monthly); 
            else
                lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):find(lowest_ind_in_radar_time_series==1)+legacy_length; 
            end     
            temp_ind2=zeros(length(localmax_01ind_radar),1);
            temp_ind2(lowest_Before6mon_ind_in_rain_time_series)=1;
            final_after6mon_ind_01=localmax_01ind_radar & temp_ind2;


            [after_drought_max,ind22]=(nanmax(ecoregion_fullRadar_1_100_normalize(final_after6mon_ind_01)));
            local_max_ind_after_drought=find(final_after6mon_ind_01==1);
            time_from_drought=local_max_ind_after_drought(ind22)-find(lowest_ind_in_radar_time_series==1);


            resilience_one_drought_temp=((nanmax(ecoregion_fullRadar_1_100_normalize(final_after6mon_ind_01)))-(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01))))/(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));

            if isempty(resilience_one_drought_temp)
                resilience_one_drought=NaN;
            else
                resilience_one_drought=resilience_one_drought_temp;
            end


        end

        this_pixel_all_drought_resistance=[this_pixel_all_drought_resistance;resistance_one_drought];
        this_pixel_all_drought_resilience=[this_pixel_all_drought_resilience;resilience_one_drought];


    end


    if ~isempty(this_pixel_all_drought_time)


    [Groups,uniqueyear]=findgroups(this_pixel_all_drought_time); 
    all_pixel_all_drought_severity(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_severity, Groups)';
    all_pixel_all_drought_radar_resistance(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_resistance, Groups)'; 
    all_pixel_all_drought_radar_resilience(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_resilience, Groups)'; 
    
    end

end   


%%%%%% 1992 1998 2018 should be removed because of insufficient data %%%%%%
all_pixel_all_drought_radar_resistance(:,[1 7 27])=NaN;
all_pixel_all_drought_radar_resilience(:,[1 7 27])=NaN;
all_pixel_all_drought_severity(:,[1 7 27])=NaN;


ecoregion_drought_severity=nanmedian(all_pixel_all_drought_severity);
ecoregion_drought_resilience=nanmedian(all_pixel_all_drought_radar_resilience);
ecoregion_drought_resistance=nanmedian(all_pixel_all_drought_radar_resistance);


%% Figure America

CT_movmean_red=cbrewer('seq', 'Oranges', 9);
radar_line2=CT_movmean_red(8,:);
line_colors=radar_line2;


disp('America......')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1=ecoregion_drought_resistance;
datain=[(1:length(data1))' data1'];
[taub tau h sig_resistance Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
disp('Americas Resistance')
disp([taub sig_resistance])
taub_resistance_allecoregions=taub;



data1=ecoregion_drought_resilience;
datain=[(1:length(data1))' data1'];
[taub tau h sig_resilience Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
disp('Americas Resilience')
disp([taub sig_resilience])
taub_resilience_allecoregions=taub;
sig_resilience=roundn(sig_resilience,-2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resistance %%%%%%%%%%%%%%%%%
pos2 = [0.05+0.08 0.32+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

plot_time=1:length(1992:2018);
plot_temp345=ecoregion_drought_resistance;

idx3=isnan(plot_temp345); 
s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,line_colors,'filled');  %%% shade_colors  line_colors
s1.MarkerFaceAlpha = scatter_alpha;
set(gca,'xtick',(1992:3:2019)-1992+1)
set(gca,'xticklabel',{})

box on

set(gca,'xlim',[0 28])
set(gca,'ytick',[-1:0.2:-0.4])
set(gca,'ylim',[-1.3 -0.2])  


set(gca,'fontsize',9)
ylabel('Resistance','fontsize',10)

xlim1=xlim;
if sig_resistance>0.05
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,-0.36,['\tau = ' num2str(roundn(taub_resistance_allecoregions,-2)) ', {\itP} > 0.05'],'fontweight','bold','color',line_colors)  
else
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,-0.36,['\tau = ' num2str(roundn(taub_resistance_allecoregions,-2)) ', {\itP} < 0.05'],'fontweight','bold','color',line_colors)  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resilience %%%%%%%%%%%%%%%%%

pos2 = [0.05+0.08 0.11+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 


plot_time=1:length(1992:2018);
plot_temp345=ecoregion_drought_resilience;

idx3=isnan(plot_temp345); 

s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,line_colors,'filled');
s1.MarkerFaceAlpha = scatter_alpha;


set(gca,'xtick',(1992:3:2019)-1992+1)
set(gca,'xticklabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016','2019'},'fontsize',9)
set(gca,'XTickLabelRotation',90)


box on


set(gca,'fontsize',9)
ylabel('Resilience','fontsize',10); 
set(gca,'xlim',[0 28])
set(gca,'ylim',[-0.19 0.22])
set(gca,'ytick',-0.1:0.1:0.1)

xlim1=xlim;
if sig_resilience>0.05 
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,0.162,['\tau = '  num2str(roundn(taub_resilience_allecoregions,-2))  ', {\itP} > 0.05'],'fontweight','bold','color',line_colors)  
else
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,0.162,['\tau = ' num2str(roundn(taub_resilience_allecoregions,-2)) ', {\itP} < 0.05'],'fontweight','bold','color',line_colors)  
end





%% Africa ##################################################################################

load Africa_tropical_belt_fishnet
load Africa_tropic_belt_HansenCover90_non_forest_count
load Africa_tropical_belt_LUCC
load Africa_tropical_belt_GLEAM_ET_19922018 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Af_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load  Africa_tropical_belt_CHIRPS_19902018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start from 1992
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHIRPS_date_final=CHIRPS_date_final(25:348);
CHIRPS_date_decimal=CHIRPS_date_decimal(25:348);
CHIRPS_data_final=CHIRPS_data_final(25:348,:);

% CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
% CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
% CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);

%%%%%%%%%%%%%%%%%%%%%%%% WD of each month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Water deficit for each 0.25 degree cell.
if ET100mm_or_not
    CHIRPS_p_evap_allmonth=CHIRPS_data_final'-100;       
else    
    CHIRPS_p_evap_allmonth=CHIRPS_data_final'-GLEAM_ET_data;    
end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% maximum WD for each year %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);


CHIRPS_MWD_year=unique(CHIRPS_hydro_year_index);
CHIRPS_MWD_allpixel=nan(size(CHIRPS_WD_allmonth,1),length(CHIRPS_MWD_year));
CHIRPS_MWD_allpixel_minmonth=nan(size(CHIRPS_WD_allmonth,1),length(CHIRPS_MWD_year));

for i=1:length(CHIRPS_MWD_year)

    one_year_allpixel_WD=CHIRPS_WD_allmonth(:,CHIRPS_hydro_year_index==CHIRPS_MWD_year(i));
    [one_year_allpixel_MWD,min_month]=nanmin(one_year_allpixel_WD,[],2);    % min WD of 12 months
    CHIRPS_MWD_allpixel(:,i)= one_year_allpixel_MWD;
    CHIRPS_MWD_allpixel_minmonth(:,i)= min_month;
end



drought_years=1992:2018; %%%%% potential drought years

zscore_mcwd_all=nan(size(CHIRPS_WD_allmonth,1),length(CHIRPS_MWD_year));
for i_droughtyear=1:length(drought_years)
    
    drought_year=drought_years(i_droughtyear);    
    CHIRPS_MWD_droughtyear=CHIRPS_MWD_allpixel(:,CHIRPS_MWD_year==drought_year);      
    CHIRPS_MWD_background=CHIRPS_MWD_allpixel(:,CHIRPS_MWD_year~=drought_year);
    Background_month_std=nanstd(CHIRPS_MWD_background,0,2);
    Background_month_mean=nanmean(CHIRPS_MWD_background,2);
    zscore_mcwd_all(:,i_droughtyear)=(CHIRPS_MWD_droughtyear-Background_month_mean)./Background_month_std;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pixel level resistance and resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load  Africa_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff

Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1);

zscore_mcwd_all=zscore_mcwd_all(ind_1,:); %%%%%%%%%%%% subset to ITR
CHIRPS_MWD_allpixel_minmonth=CHIRPS_MWD_allpixel_minmonth(ind_1,:); %%%%%%%%%%%% subset to ITR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_pixel_all_drought_severity=nan(size(Final_time_series_FullRadar,2),length(1992:2018));                          
all_pixel_all_drought_radar_resistance=nan(size(Final_time_series_FullRadar,2),length(1992:2018));
all_pixel_all_drought_radar_resilience=nan(size(Final_time_series_FullRadar,2),length(1992:2018));



for i=1:1:size(zscore_mcwd_all,1)


    ecoregion_fullRadar=Final_time_series_FullRadar(:,i);

    if detrend_or_not

        %%%%%%%%%%%%%%% Detrend won't take NaN %%%%%%%%%%%%%%
        ecoregion_fullRadar_detrend=detrend_nan(ecoregion_fullRadar);
        ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar_detrend,1,100);

    else

        ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar,1,100);

    end
    

    each_drought_positin = CHIRPS_MWD_allpixel_minmonth(i,:);


    this_pixel_all_drought_severity=[];
    this_pixel_all_drought_resistance=[];
    this_pixel_all_drought_resilience=[];
    
    this_pixel_all_drought_time=[];

    localmax_01ind_radar=islocalmax(ecoregion_fullRadar);

    for iiidrought=1:length(1992:2018)


        drought_center_time_1=(1991+iiidrought)*100+each_drought_positin(iiidrought);
        drought_center_time=(1991+iiidrought)+each_drought_positin(iiidrought)/12.0;

        drought_index_in_rain=find(CHIRPS_date_final==drought_center_time_1);


        if floor(drought_center_time)>=2019
            continue
        end



        one_drought_severity=zscore_mcwd_all(i,iiidrought);  

        if one_drought_severity>mcwdstd_thre
            continue
        end

        one_drought_severity=abs(one_drought_severity);


        this_pixel_all_drought_severity=[this_pixel_all_drought_severity;one_drought_severity];
        this_pixel_all_drought_time=[this_pixel_all_drought_time;floor(drought_center_time)];
        
        findmin_length=3; %%% Minimum radar value with +- 3 months of central drought time

        if min(drought_index_in_rain)-findmin_length<=0
            start_index=1;
        else
            start_index=min(drought_index_in_rain)-findmin_length;
        end

        if max(drought_index_in_rain)+findmin_length>length(CHIRPS_date_final)
            end_index=length(CHIRPS_date_final);
        else
            end_index=max(drought_index_in_rain)+findmin_length;
        end

        drought_time_in_rain=CHIRPS_date_final(start_index:end_index);
        drought_time_in_rain_decimal=CHIRPS_date_decimal(start_index:end_index);

        if sum(ismember(Final_time_monthly,drought_time_in_rain))==0

            Lowest_radar_value_in_drought=NaN;
            resistance_one_drought=NaN;
            resilience_one_drought=NaN;
            recover_time=NaN;

        else
            [Lowest_radar_value_in_drought,indmin]=nanmin(ecoregion_fullRadar_1_100_normalize(ismember(Final_time_monthly,drought_time_in_rain)));           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% resistance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lowest_ind_in_radar_time_series=Final_time_monthly==drought_time_in_rain(indmin);

            if find(lowest_ind_in_radar_time_series==1)<=legacy_length
                lowest_Before6mon_ind_in_rain_time_series=1:find(lowest_ind_in_radar_time_series==1); 
            else
                lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1)-legacy_length:find(lowest_ind_in_radar_time_series==1); 
            end

            temp_ind2=zeros(length(localmax_01ind_radar),1);
            temp_ind2(lowest_Before6mon_ind_in_rain_time_series)=1;
            final_Before6mon_ind_01=localmax_01ind_radar & temp_ind2;


            [before_drought_max,ind22]=(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));

            local_max_ind_before_drought=find(final_Before6mon_ind_01==1);
            time_to_drought=find(lowest_ind_in_radar_time_series==1)-local_max_ind_before_drought(ind22);

            if isempty(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)))
                recover_time=NaN;
            else
                temp2233=ecoregion_fullRadar_1_100_normalize>nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01));
                ind_recovertime=find(temp2233==1)-find(lowest_ind_in_radar_time_series==1);
                recover_time=min(ind_recovertime(ind_recovertime>0));

                if isempty(recover_time)
                    recover_time=NaN;
                end

            end

            resistance_one_drought_temp=((Lowest_radar_value_in_drought)-(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01))))/(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));
           if isempty(resistance_one_drought_temp)
                resistance_one_drought=NaN;
            else
                resistance_one_drought=resistance_one_drought_temp;
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(Final_time_monthly)-find(lowest_ind_in_radar_time_series==1)<=legacy_length
                lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):length(Final_time_monthly); 
            else
                lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):find(lowest_ind_in_radar_time_series==1)+legacy_length; 
            end     
            temp_ind2=zeros(length(localmax_01ind_radar),1);
            temp_ind2(lowest_Before6mon_ind_in_rain_time_series)=1;
            final_after6mon_ind_01=localmax_01ind_radar & temp_ind2;


            [after_drought_max,ind22]=(nanmax(ecoregion_fullRadar_1_100_normalize(final_after6mon_ind_01)));
            local_max_ind_after_drought=find(final_after6mon_ind_01==1);
            time_from_drought=local_max_ind_after_drought(ind22)-find(lowest_ind_in_radar_time_series==1);


            resilience_one_drought_temp=((nanmax(ecoregion_fullRadar_1_100_normalize(final_after6mon_ind_01)))-(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01))))/(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));

            if isempty(resilience_one_drought_temp)
                resilience_one_drought=NaN;
            else
                resilience_one_drought=resilience_one_drought_temp;
            end


        end

        this_pixel_all_drought_resistance=[this_pixel_all_drought_resistance;resistance_one_drought];
        this_pixel_all_drought_resilience=[this_pixel_all_drought_resilience;resilience_one_drought];


    end

    if ~isempty(this_pixel_all_drought_time) 

    [Groups,uniqueyear]=findgroups(this_pixel_all_drought_time); 
    all_pixel_all_drought_severity(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_severity, Groups)';
    all_pixel_all_drought_radar_resistance(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_resistance, Groups)'; 
    all_pixel_all_drought_radar_resilience(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_resilience, Groups)'; 
    
    end


end



%%%%%%% 1992 2018 should be removed because of insufficient data %%%%%%

all_pixel_all_drought_radar_resistance(:,[1  27])=NaN;
all_pixel_all_drought_radar_resilience(:,[1  27])=NaN;
all_pixel_all_drought_severity(:,[1 27])=NaN;


ecoregion_drought_severity=nanmedian(all_pixel_all_drought_severity);
ecoregion_drought_resilience=nanmedian(all_pixel_all_drought_radar_resilience);
ecoregion_drought_resistance=nanmedian(all_pixel_all_drought_radar_resistance);


%% Figure Africa

CT_movmean_red=cbrewer('seq', 'BuGn', 9);
radar_line2=CT_movmean_red(8,:);
line_colors=radar_line2;


disp('Africa......')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1=ecoregion_drought_resistance;
datain=[(1:length(data1))' data1'];
[taub tau h sig_resistance Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
disp('Africa Resistance')
disp([taub sig_resistance])
taub_resistance_allecoregions=taub;


data1=ecoregion_drought_resilience;
datain=[(1:length(data1))' data1'];
[taub tau h sig_resilience Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
disp('Africa Resilience')
disp([taub sig_resilience])
taub_resilience_allecoregions=taub;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resistance %%%%%%%%%%%%%%%%%


pos2 = [0.32+0.08 0.32+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 


plot_time=1:length(1992:2018);
plot_temp345=ecoregion_drought_resistance;

idx3=isnan(plot_temp345); 

s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,line_colors,'filled');
s1.MarkerFaceAlpha = scatter_alpha;
hold on


set(gca,'xtick',(1992:3:2019)-1992+1)
set(gca,'xticklabel',{})
set(gca,'XTickLabelRotation',90)

set(gca,'xlim',[0 28])
box on

set(gca,'fontsize',9)
set(gca,'ytick',[-1:0.2:-0.4])
set(gca,'ylim',[-1.3 -0.2]) 
set(gca,'yticklabel',{})



xlim1=xlim;
if sig_resistance>0.05
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,-0.36,['\tau = ' num2str(roundn(taub_resistance_allecoregions,-2)) ', {\itP} > 0.05'],'fontweight','bold','color',line_colors)  
else
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,-0.36,['\tau = ' num2str(roundn(taub_resistance_allecoregions,-2)) ', {\itP} < 0.05'],'fontweight','bold','color',line_colors)  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resilience %%%%%%%%%%%%%%%%%

pos2 = [0.32+0.08 0.11+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

plot_time=1:length(1992:2018);
plot_temp345=ecoregion_drought_resilience;

idx3=isnan(plot_temp345); 

s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,line_colors,'filled');
s1.MarkerFaceAlpha = scatter_alpha;

set(gca,'xtick',(1992:3:2019)-1992+1)
set(gca,'xticklabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016','2019'},'fontsize',9)
set(gca,'XTickLabelRotation',90)


box on

set(gca,'xlim',[0 28])
set(gca,'fontsize',9)
set(gca,'ylim',[-0.19 0.22])
set(gca,'ytick',-0.1:0.1:0.1)
set(gca,'yticklabel',{})



xlim1=xlim;
if sig_resilience>0.05
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,0.162,['\tau = '  num2str(roundn(taub_resilience_allecoregions,-2))  ', {\itP} > 0.05'],'fontweight','bold','color',line_colors)  
else
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,0.162,['\tau = ' num2str(roundn(taub_resilience_allecoregions,-2)) ', {\itP} < 0.05'],'fontweight','bold','color',line_colors)  
end





%% Asia ##################################################################################


load  Asia_tropical_belt_fishnet
load  Asia_tropic_belt_HansenCover90_non_forest_count
load  Asia_tropical_belt_LUCC
load  Asia_tropical_belt_GLEAM_ET_19922018


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Asia_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(sum(ind_1))

load  Asia_tropical_belt_CHIRPS_19902018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WD calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start from 1992
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHIRPS_date_final=CHIRPS_date_final(25:348);
CHIRPS_date_decimal=CHIRPS_date_decimal(25:348);
CHIRPS_data_final=CHIRPS_data_final(25:348,:);

% CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
% CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
% CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);

%%%%%%%%%%%%%%%%%%%%%%%% WD of each month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Water deficit for each 0.25 degree cell.
if ET100mm_or_not
    CHIRPS_p_evap_allmonth=CHIRPS_data_final'-100;        
else    
    CHIRPS_p_evap_allmonth=CHIRPS_data_final'-GLEAM_ET_data;    
end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% maximum WD for each year %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHIRPS_hydro_mon_index=repmat((1:12),1,length(1992:2018))';
CHIRPS_hydro_year_index=repmat((1992:2018),12,1);
CHIRPS_hydro_year_index=reshape(CHIRPS_hydro_year_index,12*length(1992:2018),1);


CHIRPS_MWD_year=unique(CHIRPS_hydro_year_index);
CHIRPS_MWD_allpixel=nan(size(CHIRPS_WD_allmonth,1),length(CHIRPS_MWD_year));
CHIRPS_MWD_allpixel_minmonth=nan(size(CHIRPS_WD_allmonth,1),length(CHIRPS_MWD_year));

for i=1:length(CHIRPS_MWD_year)
    one_year_allpixel_WD=CHIRPS_WD_allmonth(:,CHIRPS_hydro_year_index==CHIRPS_MWD_year(i));
    [one_year_allpixel_MWD,min_month]=nanmin(one_year_allpixel_WD,[],2);    % min WD of 12 months
    CHIRPS_MWD_allpixel(:,i)= one_year_allpixel_MWD;
    CHIRPS_MWD_allpixel_minmonth(:,i)= min_month;
end



drought_years=1992:2018; %%%%% potential drought years

zscore_mcwd_all=nan(size(CHIRPS_WD_allmonth,1),length(CHIRPS_MWD_year));

for i_droughtyear=1:length(drought_years)
    
    drought_year=drought_years(i_droughtyear);   
    CHIRPS_MWD_droughtyear=CHIRPS_MWD_allpixel(:,CHIRPS_MWD_year==drought_year);
    CHIRPS_MWD_background=CHIRPS_MWD_allpixel(:,CHIRPS_MWD_year~=drought_year);
    Background_month_std=nanstd(CHIRPS_MWD_background,0,2);
    Background_month_mean=nanmean(CHIRPS_MWD_background,2);
    zscore_mcwd_all(:,i_droughtyear)=(CHIRPS_MWD_droughtyear-Background_month_mean)./Background_month_std;


end




load  Asia_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1);


zscore_mcwd_all=zscore_mcwd_all(ind_1,:); %%%%%%%%%%%% subset to ITR
CHIRPS_MWD_allpixel_minmonth=CHIRPS_MWD_allpixel_minmonth(ind_1,:); %%%%%%%%%%%% subset to ITR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_pixel_all_drought_severity=nan(size(Final_time_series_FullRadar,2),length(1992:2018));                          
all_pixel_all_drought_radar_resistance=nan(size(Final_time_series_FullRadar,2),length(1992:2018));
all_pixel_all_drought_radar_resilience=nan(size(Final_time_series_FullRadar,2),length(1992:2018));




for i=1:1:size(zscore_mcwd_all,1)


    ecoregion_fullRadar=Final_time_series_FullRadar(:,i);

    if detrend_or_not

        %%%%%%%%%%%%%%% Detrend won't take NaN %%%%%%%%%%%%%%
        ecoregion_fullRadar_detrend=detrend_nan(ecoregion_fullRadar);
        ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar_detrend,1,100);

    else

        ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar,1,100);

    end

 
    each_drought_positin = CHIRPS_MWD_allpixel_minmonth(i,:);


    this_pixel_all_drought_severity=[];
    this_pixel_all_drought_resistance=[];
    this_pixel_all_drought_resilience=[];

    this_pixel_all_drought_time=[];

    localmax_01ind_radar=islocalmax(ecoregion_fullRadar);

    for iiidrought=1:length(1992:2018)


        drought_center_time_1=(1991+iiidrought)*100+each_drought_positin(iiidrought);
        drought_center_time=(1991+iiidrought)+each_drought_positin(iiidrought)/12.0;

        drought_index_in_rain=find(CHIRPS_date_final==drought_center_time_1);


        if floor(drought_center_time)>=2019
            continue
        end


        one_drought_severity=zscore_mcwd_all(i,iiidrought);  

        if one_drought_severity>mcwdstd_thre
            continue
        end

        one_drought_severity=abs(one_drought_severity);


        this_pixel_all_drought_severity=[this_pixel_all_drought_severity;one_drought_severity];
        this_pixel_all_drought_time=[this_pixel_all_drought_time;floor(drought_center_time)];
        
        
        findmin_length=3; %%% Minimum radar value with +- 3 months of central drought time


        if min(drought_index_in_rain)-findmin_length<=0
            start_index=1;
        else
            start_index=min(drought_index_in_rain)-findmin_length;
        end

        if max(drought_index_in_rain)+findmin_length>length(CHIRPS_date_final)
            end_index=length(CHIRPS_date_final);
        else
            end_index=max(drought_index_in_rain)+findmin_length;
        end

        drought_time_in_rain=CHIRPS_date_final(start_index:end_index);
        drought_time_in_rain_decimal=CHIRPS_date_decimal(start_index:end_index);

        if sum(ismember(Final_time_monthly,drought_time_in_rain))==0

            Lowest_radar_value_in_drought=NaN;
            resistance_one_drought=NaN;
            resilience_one_drought=NaN;
            recover_time=NaN;

        else
            [Lowest_radar_value_in_drought,indmin]=nanmin(ecoregion_fullRadar_1_100_normalize(ismember(Final_time_monthly,drought_time_in_rain)));           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% resistance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lowest_ind_in_radar_time_series=Final_time_monthly==drought_time_in_rain(indmin);

            if find(lowest_ind_in_radar_time_series==1)<=legacy_length
                lowest_Before6mon_ind_in_rain_time_series=1:find(lowest_ind_in_radar_time_series==1); 
            else
                lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1)-legacy_length:find(lowest_ind_in_radar_time_series==1); 
            end

            temp_ind2=zeros(length(localmax_01ind_radar),1);
            temp_ind2(lowest_Before6mon_ind_in_rain_time_series)=1;
            final_Before6mon_ind_01=localmax_01ind_radar & temp_ind2;


            [before_drought_max,ind22]=(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));

            local_max_ind_before_drought=find(final_Before6mon_ind_01==1);
            time_to_drought=find(lowest_ind_in_radar_time_series==1)-local_max_ind_before_drought(ind22);

            if isempty(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)))
                recover_time=NaN;
            else
                temp2233=ecoregion_fullRadar_1_100_normalize>nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01));
                ind_recovertime=find(temp2233==1)-find(lowest_ind_in_radar_time_series==1);
                recover_time=min(ind_recovertime(ind_recovertime>0));

                if isempty(recover_time)
                    recover_time=NaN;
                end

            end

            resistance_one_drought_temp=((Lowest_radar_value_in_drought)-(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01))))/(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));
            if isempty(resistance_one_drought_temp)
                resistance_one_drought=NaN;
            else
                resistance_one_drought=resistance_one_drought_temp;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(Final_time_monthly)-find(lowest_ind_in_radar_time_series==1)<=legacy_length
                lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):length(Final_time_monthly); 
            else
                lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):find(lowest_ind_in_radar_time_series==1)+legacy_length; 
            end     
            temp_ind2=zeros(length(localmax_01ind_radar),1);
            temp_ind2(lowest_Before6mon_ind_in_rain_time_series)=1;
            final_after6mon_ind_01=localmax_01ind_radar & temp_ind2;


            [after_drought_max,ind22]=(nanmax(ecoregion_fullRadar_1_100_normalize(final_after6mon_ind_01)));
            local_max_ind_after_drought=find(final_after6mon_ind_01==1);
            time_from_drought=local_max_ind_after_drought(ind22)-find(lowest_ind_in_radar_time_series==1);


            resilience_one_drought_temp=((nanmax(ecoregion_fullRadar_1_100_normalize(final_after6mon_ind_01)))-(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01))))/(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));

            if isempty(resilience_one_drought_temp)
                resilience_one_drought=NaN;
            else
                resilience_one_drought=resilience_one_drought_temp;
            end


        end

        this_pixel_all_drought_resistance=[this_pixel_all_drought_resistance;resistance_one_drought];
        this_pixel_all_drought_resilience=[this_pixel_all_drought_resilience;resilience_one_drought];

    end


    if ~isempty(this_pixel_all_drought_time) 


    [Groups,uniqueyear]=findgroups(this_pixel_all_drought_time); 
    all_pixel_all_drought_severity(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_severity, Groups)';
    all_pixel_all_drought_radar_resistance(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_resistance, Groups)'; 
    all_pixel_all_drought_radar_resilience(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_resilience, Groups)'; 
  
    end


end

%%%%%% 1992 2018 should be removed because of insufficient data %%%%%%

all_pixel_all_drought_radar_resistance(:,[1  27])=NaN;
all_pixel_all_drought_radar_resilience(:,[1  27])=NaN;
all_pixel_all_drought_severity(:,[1  27])=NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% color each bar by drought severity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ecoregion_drought_severity=nanmedian(all_pixel_all_drought_severity);
ecoregion_drought_resilience=nanmedian(all_pixel_all_drought_radar_resilience);
ecoregion_drought_resistance=nanmedian(all_pixel_all_drought_radar_resistance);


%% Figure Asia

 
CT_movmean_red=cbrewer('seq', 'PuBu', 9);
radar_line2=CT_movmean_red(8,:);
line_colors=radar_line2;


disp('Asia......')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1=ecoregion_drought_resistance;
datain=[(1:length(data1))' data1'];
[taub tau h sig_resistance Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
disp('Asia Resistance')
disp([taub sig_resistance])
taub_resistance_allecoregions=taub;



data1=ecoregion_drought_resilience;
datain=[(1:length(data1))' data1'];
[taub tau h sig_resilience Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
disp('Asia Resilience')
disp([taub sig_resilience])
taub_resilience_allecoregions=taub;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resistance %%%%%%%%%%%%%%%%%

pos2 = [0.59+0.08 0.32+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 


plot_time=1:length(1992:2018);
plot_temp345=ecoregion_drought_resistance;

idx3=isnan(plot_temp345); 

s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,line_colors,'filled');
s1.MarkerFaceAlpha = scatter_alpha;


set(gca,'xtick',(1992:3:2019)-1992+1)
set(gca,'xticklabel',{})
set(gca,'XTickLabelRotation',90)

box on
set(gca,'fontsize',9)


set(gca,'ytick',[-1:0.2:-0.4])
set(gca,'ylim',[-1.3 -0.2]) 
set(gca,'yticklabel',{})
set(gca,'xlim',[0 28])


xlim1=xlim;
if sig_resistance>0.05
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,-0.36,['\tau = ' num2str(roundn(taub_resistance_allecoregions,-2)) ', {\itP} > 0.05'],'fontweight','bold','color',line_colors)  
else
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,-0.36,['\tau = ' num2str(roundn(taub_resistance_allecoregions,-2)) ', {\itP} < 0.05'],'fontweight','bold','color',line_colors)  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resilience %%%%%%%%%%%%%%%%%

pos2 = [0.59+0.08 0.11+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 


plot_time=1:length(1992:2018);
plot_temp345=ecoregion_drought_resilience;

idx3=isnan(plot_temp345); 
s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,line_colors,'filled');
s1.MarkerFaceAlpha = scatter_alpha;


box on
set(gca,'fontsize',9)


set(gca,'ylim',[-0.19 0.22])
set(gca,'ytick',-0.1:0.1:0.1)

set(gca,'yticklabel',{})
set(gca,'xlim',[0 28])

set(gca,'xtick',(1992:3:2019)-1992+1)
set(gca,'xticklabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016','2019'},'fontsize',9)
set(gca,'XTickLabelRotation',90)

xlim1=xlim;
if sig_resilience>0.05
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,0.162,['\tau = '  num2str(roundn(taub_resilience_allecoregions,-2))  ', {\itP} > 0.05'],'fontweight','bold','color',line_colors)  
else
    text(xlim1(1)+(xlim1(2)-xlim1(1))*0.42,0.162,['\tau = ' num2str(roundn(taub_resilience_allecoregions,-2)) ', {\itP} < 0.05'],'fontweight','bold','color',line_colors)  
end



%%
tightfig
% print(gcf,'-dtiff','-r300',strcat('./PNAS_R2_Fig4.tif'))

