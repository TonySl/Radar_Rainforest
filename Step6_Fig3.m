%% 
%
% This script will create Fig. 3 in Tao et al. 

% Author: Shengli Tao
% Email: sltao1990@gmail.com
% Date: First version in 08.2018. Formatted in 02.2021.

% clc
clear

addpath('./util')
warning off


figure
set(gcf,'position',[46  46   950   550],...
                             'color','w','paperpositionmode','auto')


before_length=12*2;  %%% Maximum radar value 2 year before a drought, to account for the drought legecy effect. Also tried 1 year and 3 years
after_length=12*2;   %%% Maximum radar value 2 year after a drought, to account for the drought legecy effect. Also tried 1 year and 3 years

findmin_length=3;  %%% Minimum radar value with +- 3 months of central drought time



%% Americas. Pixel and regional levels drought resistance and resilience

load ./mat/Americas_tropical_belt_fishnet %%% Pixel ID et al...
load ./mat/Americas_tropic_belt_HansenCover90_non_forest_count  %%% For each pixel, we calculated its non-forest ratio using Hansen et al.' dataset.
load ./mat/Americas_tropical_belt_LUCC  %%% For each pixel, we calculated its dominant LUCC using ESA CCI landcover images.
load ./mat/Americas_tropical_belt_GLEAM_ET_19922018   %%% For each pixel, we calculated its monthly ET using GLEAM data

%%%% Labels for intact tropical rainforest 
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal

load ./mat/Americas_tropical_belt_Ecoregion  %%% For each pixel, we labeled its belonging ecoregion.
Olson_ID=cell_dominantEcoregion;

load ./mat/Americas_tropical_belt_CHIRPS_19902018 %%% For each pixel, we calculated its monthly precipitation using CHIRPS data

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



wd_std_threhold=-1.0; %%%%% threshold to define drought


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Pixel and regional levels drought resistance and resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

load ./mat/Americas_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1); %%%%%%%%%%%% subset to intact evergreen forests
Olson_ID=Olson_ID(ind_1); %%%%%%%%%%%% subset to intact evergreen forests
CHIRPS_WD_allmonth=CHIRPS_WD_allmonth(ind_1,:); %%%%%%%%%%%% subset to intact evergreen forests

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subregions were defined by aggregating ecoregions %%%%%%%%%%%%%%%%%%%%%%%%% 

ecoregion_id_americas{1,1}=[60142;  60163;  60107; 60132; 60158; 60143]; %%%%  Northwest Amazonia
ecoregion_id_americas{2,1}=[60125; 60182; 60124;60169;60173]; %%%% East Amazonia (Guiana Shield)
ecoregion_id_americas{3,1}=[60135; 60168; 60180; 60170; 60157; 60140; 60180; 60170; 60138; 60212]; %%%% South Amazonia
ecoregion_id_americas{4,1}=[60166;  60133; 60157]; %%%% Southwest Amazonia


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
        
        %%%%%%%%%%%%%%% Detrend won't take NaN %%%%%%%%%%%%%%
%         ecoregion_fullRadar_1_100_normalize=detrend_nan(ecoregion_fullRadar); %%%% also tried detrending the radar time series
        ecoregion_fullRadar_1_100_normalize=ecoregion_fullRadar;
        ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar_1_100_normalize,1,100); %%%% rescale the radar time series
        
        wd_one=CHIRPS_WD_allmonth(i,:);

        wd_std=(wd_one-nanmean(wd_one))/nanstd(wd_one);
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get each drought %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

            if floor(drought_center_time)>=2019 %%%% our dataset ends in 2018.xx
                continue
            end
                         
    
            one_drought_severity=abs(nansum(wd_std(drought_index_in_rain)));


            this_pixel_all_drought_severity=[this_pixel_all_drought_severity;one_drought_severity];
            this_pixel_all_drought_time=[this_pixel_all_drought_time;floor(drought_center_time)];


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

                if find(lowest_ind_in_radar_time_series==1)<=before_length
                    lowest_Before6mon_ind_in_rain_time_series=1:find(lowest_ind_in_radar_time_series==1); 
                else
                    lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1)-before_length:find(lowest_ind_in_radar_time_series==1); 
                end

                temp_ind2=zeros(length(localmax_01ind_radar),1);
                temp_ind2(lowest_Before6mon_ind_in_rain_time_series)=1;
                final_Before6mon_ind_01=localmax_01ind_radar & temp_ind2;


                [before_drought_max,~]=(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));

                local_max_ind_before_drought=find(final_Before6mon_ind_01==1);
                

                resistance_one_drought_temp=((Lowest_radar_value_in_drought)-(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01))))/(nanmax(ecoregion_fullRadar_1_100_normalize(final_Before6mon_ind_01)));

                if isempty(resistance_one_drought_temp)
    %                 disp('resilience empty....')
                    resistance_one_drought=NaN;
                    
                else
                    resistance_one_drought=resistance_one_drought_temp;

                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if length(Final_time_monthly)-find(lowest_ind_in_radar_time_series==1)<=after_length
                    lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):length(Final_time_monthly); 
                else
                    lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):find(lowest_ind_in_radar_time_series==1)+after_length; 
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


        if ~isempty(this_pixel_all_drought_time) %%%% if drought happened to this pixel... Some pixels didn't suffer from droughts at all
    

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


%% Figure America


scatter_size=30;
scatter_alpha=0.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mann-Kendall test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taub_resistance_allecoregions=[];
for i=1:size(ecoregion_drought_resistance,1)

    data1=ecoregion_drought_resistance(i,:);
    datain=[(1:length(data1))' data1'];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
    disp('Americas')
    disp([taub sig])
    
    taub_resistance_allecoregions=[taub_resistance_allecoregions taub];

end

taub_resilience_allecoregions=[];
for i=1:size(ecoregion_drought_resilience,1)

    data1=ecoregion_drought_resilience(i,:);
    datain=[(1:length(data1))' data1'];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
    disp('Americas')
    disp([taub sig])

    taub_resilience_allecoregions=[taub_resilience_allecoregions taub];

end

reds = cbrewer('seq','OrRd',10);
colors=reds([3 5 8 10],:); 

line_colors=colors;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resistance %%%%%%%%%%%%%%%%%

pos2 = [0.05+0.08 0.32+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

for i=1:size(ecoregion_drought_resistance,1)
        plot_time=1:length(1992:2018);
        plot_temp345=ecoregion_drought_resistance(i,:);

        idx3=isnan(plot_temp345); 

        s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,'MarkerFaceColor',line_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;
        hold on

        set(gca,'xtick',(1992:3:2019)-1992+1)

        set(gca,'xticklabel',{})

        hold on

end

box on
set(gca,'xlim',[0 28])
set(gca,'ytick',[-0.9 -0.7 -0.5 -0.3])
set(gca,'ylim',[-1.1 -0.2])
set(gca,'fontsize',9)
ylabel('Resistance','fontsize',10)


text(4,-0.3,sprintf('%.2f',taub_resistance_allecoregions(1)),'color',line_colors(1,:),'fontsize',9,'fontweight','bold')  
text(9.7,-0.3,sprintf('%.2f',taub_resistance_allecoregions(2)),'color',line_colors(2,:),'fontsize',9,'fontweight','bold')  
text(15.4,-0.3,sprintf('%.2f',taub_resistance_allecoregions(3)),'color',line_colors(3,:),'fontsize',9,'fontweight','bold') 
text(21,-0.3,sprintf('%.2f',taub_resistance_allecoregions(4)),'color',line_colors(4,:),'fontsize',9,'fontweight','bold')  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resilience %%%%%%%%%%%%%%%%%

pos2 = [0.05+0.08 0.11+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 


for i=1:size(ecoregion_drought_resilience,1)
        plot_time=1:length(1992:2018);
        plot_temp345=ecoregion_drought_resilience(i,:);
     
        idx3=isnan(plot_temp345); 
         
        s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,'MarkerFaceColor',line_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;
        hold on
 
        
        set(gca,'xtick',(1992:3:2019)-1992+1)
        set(gca,'xticklabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016','2019'},'fontsize',9)

        set(gca,'XTickLabelRotation',90)

        
        
        hold on

end

box on
set(gca,'fontsize',9)
ylabel('Resilience','fontsize',10);
set(gca,'xlim',[0 28])
set(gca,'ylim',[-0.19 0.22])
set(gca,'ytick',-0.1:0.1:0.1)


text(4,0.17,sprintf('%.2f',taub_resilience_allecoregions(1)),'color',line_colors(1,:),'fontsize',9,'fontweight','bold')  
text(9.7,0.17,sprintf('%.2f',taub_resilience_allecoregions(2)),'color',line_colors(2,:),'fontsize',9,'fontweight','bold')  
text(15.4,0.17,sprintf('%.2f',taub_resilience_allecoregions(3)),'color',line_colors(3,:),'fontsize',9,'fontweight','bold')  
text(21,0.17,sprintf('%.2f',taub_resilience_allecoregions(4)),'color',line_colors(4,:),'fontsize',9,'fontweight','bold')  



%% Africa. Pixel and regional levels drought resistance and resilience

load ./mat/Africa_tropical_belt_fishnet %%%%% See above section (Americas. Pixel and regional levels drought resistance and resilience) for explanations. Same for codes below.
load ./mat/Africa_tropic_belt_HansenCover90_non_forest_count
load ./mat/Africa_tropical_belt_LUCC
load ./mat/Africa_tropical_belt_GLEAM_ET_19922018  

ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal


load ./mat/Africa_tropical_belt_Ecoregion
Olson_ID=cell_dominantEcoregion;

load ./mat/Africa_tropical_belt_CHIRPS_19902018

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


wd_std_threhold=-1.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Pixel and regional levels drought resistance and resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ./mat/Africa_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1);
Olson_ID=Olson_ID(ind_1);
CHIRPS_WD_allmonth=CHIRPS_WD_allmonth(ind_1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Three subregions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecoregion_id_africa{1,1}=30124;
ecoregion_id_africa{2,1}=[30104 30129];
ecoregion_id_africa{3,1}=30126;

ecoregion_name_africa{1,1}={'Northeastern'; 'Congo'};
ecoregion_name_africa{2,1}='Central Congo';
ecoregion_name_africa{3,1}='Northwestern Congo';

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
        
        
        %%%%%%%%%%%%%%% Detrend won't take NaN %%%%%%%%%%%%%%
%         ecoregion_fullRadar_1_100_normalize=detrend_nan(ecoregion_fullRadar);
        ecoregion_fullRadar_1_100_normalize=ecoregion_fullRadar;
        ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar_1_100_normalize,1,100);

        wd_one=CHIRPS_WD_allmonth(i,:);
        wd_std=(CHIRPS_WD_allmonth(i,:)-nanmean(CHIRPS_WD_allmonth(i,:)))/nanstd(CHIRPS_WD_allmonth(i,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get each drought %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

                if find(lowest_ind_in_radar_time_series==1)<=before_length
                    lowest_Before6mon_ind_in_rain_time_series=1:find(lowest_ind_in_radar_time_series==1); 
                else
                    lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1)-before_length:find(lowest_ind_in_radar_time_series==1); 
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
                if length(Final_time_monthly)-find(lowest_ind_in_radar_time_series==1)<=after_length
                    lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):length(Final_time_monthly); 
                else
                    lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):find(lowest_ind_in_radar_time_series==1)+after_length; 
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


    %%%%%% 1992 and 2018 are the first and last year, respectively. Resistance and resilience can't be calculated.
    all_pixel_all_drought_radar_resistance(:,[1 27])=NaN;
    all_pixel_all_drought_radar_resilience(:,[1 27])=NaN;
    all_pixel_all_drought_severity(:,[1 27])=NaN;    

    drought_severity=nanmedian(all_pixel_all_drought_severity);  

    ecoregion_drought_severity=[ecoregion_drought_severity;drought_severity];    
    ecoregion_drought_resilience=[ecoregion_drought_resilience;nanmedian(all_pixel_all_drought_radar_resilience)];
    ecoregion_drought_resistance=[ecoregion_drought_resistance;nanmedian(all_pixel_all_drought_radar_resistance)];
    
end


%% Figure Africa

scatter_size=30;
scatter_alpha=0.7;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mann-Kendall test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taub_resistance_allecoregions=[];
for i=1:size(ecoregion_drought_resistance,1)

    data1=ecoregion_drought_resistance(i,:);
    datain=[(1:length(data1))' data1'];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
    disp('Africa')
    disp([taub sig])
    
    taub_resistance_allecoregions=[taub_resistance_allecoregions taub];

end

taub_resilience_allecoregions=[];
for i=1:size(ecoregion_drought_resilience,1)

    data1=ecoregion_drought_resilience(i,:);
    datain=[(1:length(data1))' data1'];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
    disp('Africa')
    disp([taub sig])

    taub_resilience_allecoregions=[taub_resilience_allecoregions taub];

end


set1 = cbrewer('qual','Set2',9);
greens = cbrewer('seq','Greens',10);
colors_africa=greens([3 6 9],:); 

line_colors=colors_africa;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resistance %%%%%%%%%%%%%%%%%

pos2 = [0.32+0.08 0.32+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

for i=1:size(ecoregion_drought_resistance,1)
        plot_time=1:length(1992:2018);
        plot_temp345=ecoregion_drought_resistance(i,:);

        idx3=isnan(plot_temp345); 

        s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,'MarkerFaceColor',line_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;
        hold on


        set(gca,'xtick',(1992:3:2019)-1992+1)
        set(gca,'xticklabel',{})
        set(gca,'XTickLabelRotation',90)
        hold on

end

set(gca,'xlim',[0 28])
box on
set(gca,'fontsize',9)
set(gca,'ytick',[-0.9 -0.7 -0.5 -0.3])
set(gca,'ylim',[-1.1 -0.2])
set(gca,'yticklabel',{})


text(7.5,-0.3,sprintf('%.2f',taub_resistance_allecoregions(1)),'color',line_colors(1,:),'fontsize',9,'fontweight','bold')  
text(13,-0.3,sprintf('%.2f',taub_resistance_allecoregions(2)),'color',line_colors(2,:),'fontsize',9,'fontweight','bold')  
text(18.5,-0.3,sprintf('%.2f',taub_resistance_allecoregions(3)),'color',line_colors(3,:),'fontsize',9,'fontweight','bold')  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resilience %%%%%%%%%%%%%%%%%

pos2 = [0.32+0.08 0.11+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

for i=1:size(ecoregion_drought_resilience,1)
        plot_time=1:length(1992:2018);
        plot_temp345=ecoregion_drought_resilience(i,:);
     
        idx3=isnan(plot_temp345); 
         
        s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,'MarkerFaceColor',line_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;
        hold on
   
        set(gca,'xtick',(1992:3:2019)-1992+1)
        set(gca,'xticklabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016','2019'},'fontsize',9)
        set(gca,'XTickLabelRotation',90)

        hold on

end
box on

set(gca,'xlim',[0 28])
set(gca,'fontsize',9)
set(gca,'ylim',[-0.19 0.22])
set(gca,'ytick',-0.1:0.1:0.1)

set(gca,'yticklabel',{})


text(7.5,0.17,sprintf('%.2f',taub_resilience_allecoregions(1)),'color',line_colors(1,:),'fontsize',9,'fontweight','bold')  
text(13,0.17,sprintf('%.2f',taub_resilience_allecoregions(2)),'color',line_colors(2,:),'fontsize',9,'fontweight','bold')  
text(18.5,0.17,sprintf('%.2f',taub_resilience_allecoregions(3)),'color',line_colors(3,:),'fontsize',9,'fontweight','bold')  




%% Asia. Pixel and regional levels drought resistance and resilience

load ./mat/Asia_tropical_belt_fishnet %%%%% See above section (Americas. Pixel and regional levels drought resistance and resilience) for explanations. Same for codes below.
load ./mat/Asia_tropic_belt_HansenCover90_non_forest_count
load ./mat/Asia_tropical_belt_LUCC
load ./mat/Asia_tropical_belt_GLEAM_ET_19922018


ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal

load ./mat/Asia_tropical_belt_Ecoregion
Olson_ID=cell_dominantEcoregion;

load ./mat/Asia_tropical_belt_CHIRPS_19902018

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


wd_std_threhold=-1.0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Pixel and regional levels drought resistance and resilience %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ./mat/Asia_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff
Final_time_series_FullRadar=Final_time_series_FullRadar(:,ind_1);
Olson_ID=Olson_ID(ind_1);
CHIRPS_WD_allmonth=CHIRPS_WD_allmonth(ind_1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Two subregions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecoregion_id_asia{1,1}=[40103 40102]; % Borneo
ecoregion_id_asia{2,1}=[10115 10105 10120 10121 10122 10128 10127]; %  New Guinea


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
        
        
         %%%%%%%%%%%%%%% Detrend won't take NaN %%%%%%%%%%%%%%
%          ecoregion_fullRadar_1_100_normalize=detrend_nan(ecoregion_fullRadar);
        ecoregion_fullRadar_1_100_normalize=ecoregion_fullRadar;
        ecoregion_fullRadar_1_100_normalize=rescale(ecoregion_fullRadar_1_100_normalize,1,100);
        
        wd_one=CHIRPS_WD_allmonth(i,:);
        wd_std=(CHIRPS_WD_allmonth(i,:)-nanmean(CHIRPS_WD_allmonth(i,:)))/nanstd(CHIRPS_WD_allmonth(i,:));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get each drought %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

                if find(lowest_ind_in_radar_time_series==1)<=before_length
                    lowest_Before6mon_ind_in_rain_time_series=1:find(lowest_ind_in_radar_time_series==1); 
                else
                    lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1)-before_length:find(lowest_ind_in_radar_time_series==1); 
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
                if length(Final_time_monthly)-find(lowest_ind_in_radar_time_series==1)<=after_length
                    lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):length(Final_time_monthly); 
                else
                    lowest_Before6mon_ind_in_rain_time_series=find(lowest_ind_in_radar_time_series==1):find(lowest_ind_in_radar_time_series==1)+after_length; 
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


    %%%%%% 1992 and 2018 are the first and last year, respectively. Resistance and resilience can't be calculated.
    all_pixel_all_drought_radar_resistance(:,[1 27])=NaN;
    all_pixel_all_drought_radar_resilience(:,[1 27])=NaN;
    all_pixel_all_drought_severity(:,[1 27])=NaN;   
    
    drought_severity=nanmedian(all_pixel_all_drought_severity);  
    
    
    ecoregion_drought_severity=[ecoregion_drought_severity;drought_severity];    
    ecoregion_drought_resilience=[ecoregion_drought_resilience;nanmedian(all_pixel_all_drought_radar_resilience)];
    ecoregion_drought_resistance=[ecoregion_drought_resistance;nanmedian(all_pixel_all_drought_radar_resistance)];
    
end

%% Figure Asia

scatter_size=30;
scatter_alpha=0.7;


taub_resistance_allecoregions=[];
for i=1:size(ecoregion_drought_resistance,1)

    data1=ecoregion_drought_resistance(i,:);
    datain=[(1:length(data1))' data1'];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
    disp('Asia')
    disp([taub sig])
    
    taub_resistance_allecoregions=[taub_resistance_allecoregions taub];

end

taub_resilience_allecoregions=[];
for i=1:size(ecoregion_drought_resilience,1)

    data1=ecoregion_drought_resilience(i,:);
    datain=[(1:length(data1))' data1'];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
    disp('Asia')
    disp([taub sig])

    taub_resilience_allecoregions=[taub_resilience_allecoregions taub];

end

blues = cbrewer('seq','Blues',9);
colors_asia=blues([4 7],:); 
line_colors=colors_asia;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resistance %%%%%%%%%%%%%%%%%

pos2 = [0.59+0.08 0.32+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

for i=1:size(ecoregion_drought_resistance,1)
        plot_time=1:length(1992:2018);
        plot_temp345=ecoregion_drought_resistance(i,:);

        idx3=isnan(plot_temp345); 

        s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,'MarkerFaceColor',line_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;
        hold on

        set(gca,'xtick',(1992:3:2019)-1992+1)
        set(gca,'xticklabel',{})
        set(gca,'XTickLabelRotation',90)

        hold on

end


box on
set(gca,'fontsize',9)
set(gca,'ytick',[-0.9 -0.7 -0.5 -0.3])
set(gca,'ylim',[-1.1 -0.2])
set(gca,'yticklabel',{})
set(gca,'xlim',[0 28])

text(10,-0.3,sprintf('%.2f',taub_resistance_allecoregions(1)),'color',line_colors(1,:),'fontsize',9,'fontweight','bold')
text(16,-0.3,sprintf('%.2f',taub_resistance_allecoregions(2)),'color',line_colors(2,:),'fontsize',9,'fontweight','bold')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% time series, resilience %%%%%%%%%%%%%%%%%

pos2 = [0.59+0.08 0.11+0.4 0.23 0.21];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

for i=1:size(ecoregion_drought_resilience,1)
        plot_time=1:length(1992:2018);
        plot_temp345=ecoregion_drought_resilience(i,:);
     
        idx3=isnan(plot_temp345); 
         
        s1=scatter(plot_time(~idx3),plot_temp345(~idx3),scatter_size,'MarkerFaceColor',line_colors(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
        s1.MarkerFaceAlpha = scatter_alpha;
        hold on
        
end

box on
set(gca,'fontsize',9)
set(gca,'ylim',[-0.19 0.22])
set(gca,'ytick',-0.1:0.1:0.1)
set(gca,'yticklabel',{})
set(gca,'xlim',[0 28])
set(gca,'xtick',(1992:3:2019)-1992+1)
set(gca,'xticklabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016','2019'},'fontsize',9)
set(gca,'XTickLabelRotation',90)


text(10,0.17,sprintf('%.2f',taub_resilience_allecoregions(1)),'color',line_colors(1,:),'fontsize',9,'fontweight','bold')
text(16,0.17,sprintf('%.2f',taub_resilience_allecoregions(2)),'color',line_colors(2,:),'fontsize',9,'fontweight','bold')


%%
tightfig
% print(gcf,'-dtiff','-r600',strcat('./Figure3_2020_GLEAM_ET_2yr_legacy_resubmit_nodetrend.tif'))
