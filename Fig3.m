
%% 
%
% This script will create Fig. 3 in Tao et al. 
% All input (and output) files can be downloaded at: https://doi.org/10.6084/m9.figshare.14061428.v3

% Author: Shengli Tao
% Email: sltao@pku.edu.cn
% Date: First version in 08.2018. Formatted in 01.2022.

clc
clear

addpath('./util')
addpath('./mat')
warning off

legacy_length=12*2;   %%% two-year drought legecy effect.

wd_std_threhold=-1.0; %%%%% standarized cwd threshold to define drought

detrend_or_not=0; %%% whether detrend the signal or not. 

regional_or_continental=0; %%% 0 regional, 1 continental



figure
set(gcf,'position',[46  46   950   550],...
                             'color','w','paperpositionmode','auto')



%% Americas. Pixel and regional levels drought resistance and resilience

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

load Americas_tropical_belt_Ecoregion %%% For each pixel, we labeled its belonging ecoregion.
Olson_ID=cell_dominantEcoregion;

load Americas_tropical_belt_CHIRPS_19902018 %%% For each pixel, we calculated its monthly precipitation using CHIRPS data

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Pixel and regional levels drought resistance and resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

load Americas_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1); %%%%%%%%%%%% subset to ITR
Olson_ID=Olson_ID(ind_1); %%%%%%%%%%%% subset to ITR
CHIRPS_WD_allmonth=CHIRPS_WD_allmonth(ind_1,:); %%%%%%%%%%%% subset to ITR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Regions were defined by aggregating ecoregions %%%%%%%%%%%%%%%%%%%%%%%%% 

ecoregion_id_americas{1,1}=[60142;  60163;  60107; 60132; 60158; 60143]; %%%%  Northwest Amazonia
ecoregion_id_americas{2,1}=[60125; 60182; 60124;60169;60173]; %%%% East Amazonia (Guiana Shield)
ecoregion_id_americas{3,1}=[60135; 60168; 60180; 60170; 60157; 60140; 60180; 60170; 60138; 60212]; %%%% South Amazonia
ecoregion_id_americas{4,1}=[60166;  60133; 60157]; %%%% Southwest Amazonia

if regional_or_continental
    
    ecoregion_id_americas{1,1}=[ecoregion_id_americas{1,1};ecoregion_id_americas{2,1};ecoregion_id_americas{3,1};ecoregion_id_americas{4,1}];
    ecoregion_id_americas{2,1}=[];
    ecoregion_id_americas{3,1}=[];
    ecoregion_id_americas{4,1}=[];
    
end

ecoregion_drought_severity=[];
ecoregion_drought_resistance=[];
ecoregion_drought_resilience=[];
ecoregion_drought_recovertime=[];



for i_ecoregion=1:length(ecoregion_id_americas)

    all_pixel_all_drought_severity=nan(size(Final_time_series_FullRadar,2),length(1992:2018));                          
    all_pixel_all_drought_radar_resistance=nan(size(Final_time_series_FullRadar,2),length(1992:2018));
    all_pixel_all_drought_radar_resilience=nan(size(Final_time_series_FullRadar,2),length(1992:2018));
   

    all_slope_resilience=nan(size(Final_time_series_FullRadar,2),1);
    all_slope_resistaance=nan(size(Final_time_series_FullRadar,2),1);

    valid_pixel_count_all=0;
    
    for i=1:1:size(Final_time_series_FullRadar,2)


        if ~ismember(Olson_ID(i),ecoregion_id_americas{i_ecoregion,1})
            continue
        end


        valid_pixel_count_all=valid_pixel_count_all+1;


        ecoregion_fullRadar=Final_time_series_FullRadar(:,i);
        
        
        
        if detrend_or_not
            
            %%%%%%%%%%%%%%% Detrend won't take NaN %%%%%%%%%%%%%%
            ecoregion_fullRadar_detrend=detrend_nan(ecoregion_fullRadar);
            ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar_detrend,1,100);
        
        else
            
            ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar,1,100);
        
        end

        wd_one=CHIRPS_WD_allmonth(i,:);
        wd_std=(wd_one-nanmean(wd_one))/nanstd(wd_one);
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get each drought  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        measurements = regionprops(logical(wd_std<wd_std_threhold), 'Area');
        theLengths = [measurements.Area];
        
        position_wd_std=find(wd_std<wd_std_threhold);
        each_drought_positin = mat2cell(position_wd_std,1,theLengths);
        

        this_pixel_all_drought_severity=[];
        this_pixel_all_drought_radar_decrease=[];
        this_pixel_all_drought_time=[];

        this_pixel_all_drought_resistance=[];
        this_pixel_all_drought_resilience=[];

        localmax_01ind_radar=islocalmax(ecoregion_fullRadar);
        
        
        for iiidrought=1:length(each_drought_positin)

            drought_index_in_rain=each_drought_positin{iiidrought};

            drought_center_time=mean(CHIRPS_date_decimal(drought_index_in_rain));

            if floor(drought_center_time)>=2019  %%%% our radar dataset ends in 2018.xx
                continue
            end
                         
    
            one_drought_severity=abs(nansum(wd_std(drought_index_in_rain))); %nanmin nansum nanmean


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

                resistance_one_drought_temp=((Lowest_radar_value_in_drought)-(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01))))/(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));
                
                if isempty(resistance_one_drought_temp)
    %                  
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

        if ~isempty(this_pixel_all_drought_time) %%%% if drought happened to this pixel, some pixels didn't suffer from droughts at all
    

        [Groups,uniqueyear]=findgroups(this_pixel_all_drought_time); %%%%% Deal with the very rare situation that two droughts occured in one year.
        %%% Take the mean values
        all_pixel_all_drought_severity(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_severity, Groups)';
        all_pixel_all_drought_radar_resistance(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_resistance, Groups)'; 
        all_pixel_all_drought_radar_resilience(i,ismember(1992:2018,uniqueyear))=splitapply(@mean, this_pixel_all_drought_resilience, Groups)'; 
        all_pixel_all_drought_year(i,ismember(1992:2018,uniqueyear))=uniqueyear;
        
        end

    end
    
    %%%%%% 1992 and 2018 are the first and last year, respectively. Resistance and resilience can't be calculated.    
    %%%%%% 1998 has too few monthly data available in tropical Americas due to ERS image quality %%%%%% 
      
    all_pixel_all_drought_radar_resistance(:,[1 7 27])=NaN;
    all_pixel_all_drought_radar_resilience(:,[1 7 27])=NaN;
    all_pixel_all_drought_severity(:,[1 7 27])=NaN;
    

   
    drought_severity=nanmedian(all_pixel_all_drought_severity);   
    


    ecoregion_drought_severity=[ecoregion_drought_severity;drought_severity];    
    ecoregion_drought_resilience=[ecoregion_drought_resilience;nanmedian(all_pixel_all_drought_radar_resilience)];
    ecoregion_drought_resistance=[ecoregion_drought_resistance;nanmedian(all_pixel_all_drought_radar_resistance)];
    
end

%% Figure Americas


scatter_alpha=0.7;
scatter_size=30;


reds = cbrewer('seq','OrRd',10);
colors=reds([3 5 8 10],:); 

pts_colors=colors;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% severity, resistance %%%%%%%%%%%%%%%%%

pos2 = [0.05+0.08 0.69+0.0865 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

for i=[1 4 2 3] %%% four regions
    

        drought_severity=ecoregion_drought_severity(i,:);        
        ploty_temp=ecoregion_drought_resistance(i,:);
                
        s1=scatter(drought_severity,ploty_temp,scatter_size,'MarkerFaceColor',pts_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]); %%%% ,symbols{i} ,'filled'
        s1.MarkerFaceAlpha = scatter_alpha;
        s1.MarkerEdgeAlpha = scatter_alpha;
        
   
        hold on

end

line([-2 16], [-1 -1],'LineStyle','--','color','k') 

box on

set(gca,'fontsize',9)
ylabel('Resistance','fontsize',10)


set(gca,'xlim',[-0.8 16])
set(gca,'xtick',0:5:15)
set(gca,'xticklabel',{})


set(gca,'ylim',[-1.1 -0.3])
set(gca,'ytick',[-1:0.2:-0.4])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% severity, resilience %%%%%%%%%%%%%%%%%  

pos2 = [0.05+0.08 0.48+0.0865 0.23 0.21];
% pos2 = [0.05+0.08 0.34+0.0865 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1])

for i=[1 4 2 3] 
    
    
        drought_severity=ecoregion_drought_severity(i,:);
        
        ploty_temp=ecoregion_drought_resilience(i,:);
        
        s1=scatter(drought_severity,ploty_temp,scatter_size,'MarkerFaceColor',pts_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]); %%% %%%% ,symbols{i}'filled'
        s1.MarkerFaceAlpha = scatter_alpha;
        s1.MarkerEdgeAlpha = scatter_alpha;
        

        hold on

end

box on
set(gca,'fontsize',9)

ylabel('Resilience','fontsize',10)
xlabel('Drought severity','fontsize',10)


set(gca,'xlim',[-0.8 16])
set(gca,'xtick',0:5:15)

set(gca,'ylim',[-0.18 0.22])
set(gca,'ytick',[-0.1:0.1:0.1])


line([-2 16], [0 0],'LineStyle','--','color','k') 
hold on




%% Africa. Pixel and regional levels drought resistance and resilience


load Africa_tropical_belt_fishnet
load Africa_tropical_belt_LUCC
load Africa_tropical_belt_GLEAM_ET_19922018 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Hansen ITR mask 
% load Africa_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<5 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal


%%%%% Vancutsem ITR mask
load Af_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
disp(sum(ind_1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load Africa_tropical_belt_Ecoregion
Olson_ID=cell_dominantEcoregion;

load Africa_tropical_belt_CHIRPS_19902018

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Pixel and regional levels drought resistance and resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Africa_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1);
Olson_ID=Olson_ID(ind_1);
CHIRPS_WD_allmonth=CHIRPS_WD_allmonth(ind_1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Three regions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecoregion_id_africa{1,1}=30124;
ecoregion_id_africa{2,1}=[30104; 30129];
ecoregion_id_africa{3,1}=30126;


if regional_or_continental
    
    ecoregion_id_africa{1,1}=[ecoregion_id_africa{1,1};ecoregion_id_africa{2,1};ecoregion_id_africa{3,1}];
    ecoregion_id_africa{2,1}=[];
    ecoregion_id_africa{3,1}=[];

    
end



ecoregion_drought_severity=[];
ecoregion_drought_resistance=[];
ecoregion_drought_resilience=[];


for i_ecoregion=1:length(ecoregion_id_africa)
    

    all_pixel_all_drought_severity=nan(size(Final_time_series_FullRadar,2),length(1992:2018));                          
    all_pixel_all_drought_radar_resistance=nan(size(Final_time_series_FullRadar,2),length(1992:2018));
    all_pixel_all_drought_radar_resilience=nan(size(Final_time_series_FullRadar,2),length(1992:2018));

    all_slope_resilience=nan(size(Final_time_series_FullRadar,2),1);
    all_slope_resistaance=nan(size(Final_time_series_FullRadar,2),1);

    valid_pixel_count_all=0;
    
    for i=1:1:size(Final_time_series_FullRadar,2)

        if ~ismember(Olson_ID(i),ecoregion_id_africa{i_ecoregion,1})
            continue
        end


        valid_pixel_count_all=valid_pixel_count_all+1;

        ecoregion_fullRadar=Final_time_series_FullRadar(:,i);
        
        
        if detrend_or_not
            
            %%%%%%%%%%%%%%% Detrend won't take NaN %%%%%%%%%%%%%%
            ecoregion_fullRadar_detrend=detrend_nan(ecoregion_fullRadar);
            ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar_detrend,1,100);
        
        else
            
            ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar,1,100);
        
        end
        
        wd_one=CHIRPS_WD_allmonth(i,:);
        wd_std=(CHIRPS_WD_allmonth(i,:)-nanmean(CHIRPS_WD_allmonth(i,:)))/nanstd(CHIRPS_WD_allmonth(i,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get each drought period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
        measurements = regionprops(logical(wd_std<wd_std_threhold), 'Area');
        theLengths = [measurements.Area];

        position_wd_std=find(wd_std<wd_std_threhold);
        each_drought_positin = mat2cell(position_wd_std,1,theLengths);
        

        this_pixel_all_drought_severity=[];
        this_pixel_all_drought_radar_decrease=[];
        this_pixel_all_drought_time=[];

        this_pixel_all_drought_resistance=[];
        this_pixel_all_drought_resilience=[];

        localmax_01ind_radar=islocalmax(ecoregion_fullRadar);

        for iiidrought=1:length(each_drought_positin)

            drought_index_in_rain=each_drought_positin{iiidrought};
            drought_center_time=mean(CHIRPS_date_decimal(drought_index_in_rain));

            if floor(drought_center_time)>=2019
                continue
            end


            one_drought_severity=abs(nansum(wd_std(drought_index_in_rain))); %nanmin nansum
            

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


    %%%%%% 1992 and 2018 are the first and last year, respectively. Resistance and resilience can't be calculated.
    all_pixel_all_drought_radar_resistance(:,[1  27])=NaN;
    all_pixel_all_drought_radar_resilience(:,[1  27])=NaN;
    all_pixel_all_drought_severity(:,[1 27])=NaN;

    

    drought_severity=nanmedian(all_pixel_all_drought_severity);     %%%% nanmean
    ecoregion_drought_severity=[ecoregion_drought_severity;drought_severity];
    ecoregion_drought_resilience=[ecoregion_drought_resilience;nanmedian(all_pixel_all_drought_radar_resilience)];
    ecoregion_drought_resistance=[ecoregion_drought_resistance;nanmedian(all_pixel_all_drought_radar_resistance)];
    
end

%% Figure Africa

scatter_alpha=0.7;
scatter_size=30;

greens = cbrewer('seq','Greens',10);
colors_africa=greens([3 6 9],:);

pts_colors=colors_africa;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% severity, resistance %%%%%%%%%%%%%%%%%

pos2 = [0.32+0.08 0.69+0.0865 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 


for i=1:size(ecoregion_drought_severity,1)
        
        drought_severity=ecoregion_drought_severity(i,:);
        
        ploty_temp=ecoregion_drought_resistance(i,:);

        s1=scatter(drought_severity,ploty_temp,scatter_size,'MarkerFaceColor',pts_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;

        
        hold on

end

line([-2 12], [-1 -1],'LineStyle','--','color','k') 

box on


set(gca,'xtick',0:5:10)
set(gca,'xlim',[-0.8 12])
set(gca,'ylim',[-1.1 -0.3])
set(gca,'ytick',[-1:0.2:-0.4])
set(gca,'yticklabel',{})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% severity, resilience %%%%%%%%%%%%%%%%%

pos2 = [0.32+0.08 0.48+0.0865 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1])

for i=1:size(ecoregion_drought_severity,1)
        
        drought_severity=ecoregion_drought_severity(i,:);
        
        ploty_temp=ecoregion_drought_resilience(i,:);
 
        s1=scatter(drought_severity,ploty_temp,scatter_size,'MarkerFaceColor',pts_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;
        
        
        hold on

end

box on


set(gca,'yticklabel',{})
set(gca,'xtick',0:5:10)
set(gca,'xlim',[-0.8 12])

set(gca,'ylim',[-0.18 0.22])
set(gca,'ytick',[-0.1:0.1:0.1])

set(gca,'fontsize',9)

hold on
line([-2 12], [0 0],'LineStyle','--','color','k') 

xlabel('Drought severity','fontsize',10)




%% Asia. Pixel and regional levels drought resistance and resilience

load Asia_tropical_belt_fishnet
load Asia_tropical_belt_LUCC
load Asia_tropical_belt_GLEAM_ET_19922018


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Hansen ITR mask 
% load Asia_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<5 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal




%%%%% Vancutsem ITR mask
load Asia_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
disp(sum(ind_1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



load Asia_tropical_belt_Ecoregion
Olson_ID=cell_dominantEcoregion;

load Asia_tropical_belt_CHIRPS_19902018

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Pixel and regional levels drought resistance and resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Asia_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1);
Olson_ID=Olson_ID(ind_1);
CHIRPS_WD_allmonth=CHIRPS_WD_allmonth(ind_1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecoregion_id_asia{1,1}=[40103;40102]; % Borneo
ecoregion_id_asia{2,1}=[10115;10105;10120;10121;10122;10128;10127]; %  New Guinea


if regional_or_continental
    
    ecoregion_id_asia{1,1}=[ecoregion_id_asia{1,1};ecoregion_id_asia{2,1}];
    ecoregion_id_asia{2,1}=[];
    
end


ecoregion_drought_severity=[];
ecoregion_drought_resistance=[];
ecoregion_drought_resilience=[];


for i_ecoregion=1:length(ecoregion_id_asia)
    

    all_pixel_all_drought_severity=nan(size(Final_time_series_FullRadar,2),length(1992:2018));                          
    all_pixel_all_drought_radar_resistance=nan(size(Final_time_series_FullRadar,2),length(1992:2018));
    all_pixel_all_drought_radar_resilience=nan(size(Final_time_series_FullRadar,2),length(1992:2018));


    all_slope_resilience=nan(size(Final_time_series_FullRadar,2),1);
    all_slope_resistaance=nan(size(Final_time_series_FullRadar,2),1);


    valid_pixel_count_all=0;
    
    for i=1:1:size(Final_time_series_FullRadar,2)

        if ~ismember(Olson_ID(i),ecoregion_id_asia{i_ecoregion,1})
            continue
        end


        valid_pixel_count_all=valid_pixel_count_all+1;
        ecoregion_fullRadar=Final_time_series_FullRadar(:,i);
        
        
        if detrend_or_not
            
            %%%%%%%%%%%%%%% Detrend won't take NaN %%%%%%%%%%%%%%
            ecoregion_fullRadar_detrend=detrend_nan(ecoregion_fullRadar);
            ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar_detrend,1,100);
        
        else
            
            ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar,1,100);
        
        end


        wd_one=CHIRPS_WD_allmonth(i,:);
        wd_std=(CHIRPS_WD_allmonth(i,:)-nanmean(CHIRPS_WD_allmonth(i,:)))/nanstd(CHIRPS_WD_allmonth(i,:));
        
        

        measurements = regionprops(logical(wd_std<wd_std_threhold), 'Area');
        theLengths = [measurements.Area];

        position_wd_std=find(wd_std<wd_std_threhold);

        each_drought_positin = mat2cell(position_wd_std,1,theLengths);

    

        this_pixel_all_drought_severity=[];
        this_pixel_all_drought_radar_decrease=[];
        this_pixel_all_drought_time=[];

        this_pixel_all_drought_resistance=[];
        this_pixel_all_drought_resilience=[];

        localmax_01ind_radar=islocalmax(ecoregion_fullRadar);

        for iiidrought=1:length(each_drought_positin)

            drought_index_in_rain=each_drought_positin{iiidrought};

            drought_center_time=mean(CHIRPS_date_decimal(drought_index_in_rain));

            if floor(drought_center_time)>=2019
                continue
            end

            one_drought_severity=abs(nansum(wd_std(drought_index_in_rain))); %nanmin nansum

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

                resistance_one_drought_temp=((Lowest_radar_value_in_drought)-(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01))))/(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));

                if isempty(resistance_one_drought_temp)
    %                 disp('resilience empty....')
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
    %                 disp('resilience empty....')
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


    
   
    all_pixel_all_drought_radar_resistance(:,[1  27])=NaN;
    all_pixel_all_drought_radar_resilience(:,[1  27])=NaN;
    all_pixel_all_drought_severity(:,[1  27])=NaN;
    
  
    drought_severity=nanmedian(all_pixel_all_drought_severity); 

    ecoregion_drought_severity=[ecoregion_drought_severity;drought_severity];
    ecoregion_drought_resilience=[ecoregion_drought_resilience;nanmedian(all_pixel_all_drought_radar_resilience)];
    ecoregion_drought_resistance=[ecoregion_drought_resistance;nanmedian(all_pixel_all_drought_radar_resistance)];
    
end

%% Figure Asia

scatter_alpha=0.7;
scatter_size=30;


blues = cbrewer('seq','Blues',9);
colors_asia=blues([4 7],:); 

pts_colors=colors_asia;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% severity, resistance %%%%%%%%%%%%%%%%%

pos2 = [0.59+0.08 0.69+0.0865 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 


for i=1:size(ecoregion_drought_severity,1)
        
        drought_severity=ecoregion_drought_severity(i,:);
        
        ploty_temp=ecoregion_drought_resistance(i,:);

        s1=scatter(drought_severity,ploty_temp,scatter_size,'MarkerFaceColor',pts_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;          

        
        hold on

end
box on

line([-2 40], [-1 -1],'LineStyle','--','color','k') 

set(gca,'ylim',[-1.1 -0.3])
set(gca,'ytick',[-1:0.2:-0.4])
set(gca,'yticklabel',{})
set(gca,'xlim',[-2 37])

set(gca,'xtick',0:10:40)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% severity, resilience %%%%%%%%%%%%%%%%%
 
pos2 = [0.59+0.08 0.48+0.0865 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1])

for i=1:size(ecoregion_drought_severity,1)
        
        drought_severity=ecoregion_drought_severity(i,:);
        
        ploty_temp=ecoregion_drought_resilience(i,:);

        s1=scatter(drought_severity,ploty_temp,scatter_size,'MarkerFaceColor',pts_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;

        hold on

end
box on

set(gca,'fontsize',9)

set(gca,'xlim',[-2 37])
set(gca,'xtick',0:10:40)

set(gca,'ylim',[-0.18 0.22])
set(gca,'ytick',[-0.1:0.1:0.1])
set(gca,'yticklabel',{})
 
hold on
line([-2 40], [0 0],'LineStyle','--','color','k') 

xlabel('Drought severity','fontsize',10)



%%
tightfig
% print(gcf,'-dtiff','-r300',strcat('./Fig3.tif'))
