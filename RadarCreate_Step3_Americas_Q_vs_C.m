%% 
%
% All ERS and QSCAT data have been scaled against ASCAT.
% This script will further address the monthly difference between C-band (ERS, ASCAT) and Ku-band (QSCAT) data, per pixel (25km * 25km)
% The input files can be downloaded at: https://doi.org/10.6084/m9.figshare.14061428.v3

% Author: Shengli Tao
% Email: sltao@pku.edu.cn
% Date: Created in 08.2018. Formatted in 01.2022.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

warning off
addpath('./util')


%% First, read in the pixel locations and IDs of each radar pixel 
%%% Here I included all tropical pixesl (~23N - 23S). Will subset to intact rainforest when drawing the figures.
load ./mat/Americas_tropical_belt_fishnet
cell_locations=[Am_cell_ID Am_cell_centers];%%% Am_cell_ID: ID of each radar pixel;  Am_cell_centers: x y center of each radar pixel.  Each row is a pixel 

%% TRMM data

load ./mat/Americas_tropical_belt_TRMM_19982018  %%% TRMM precipitation for each radar pixel. Seasonality of precipication and water deficit are also provided.
%%% each row is a pixel (22738), each column represents a month

%% Load Radar data 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QSCAT scaled against ASCAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('./mat/Americas_tropic_belt_Qscat_scaled_ASCAT_NoStripC20_nocorrection_withfirstMon_Linear.mat')

Q_data_scaled_no_corrected=scaled_QSCAT';
% Q_data_original_data=orginal_QSCAT';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERS scaled against scaled_QSCAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('./mat/Americas_tropic_belt_ERS2_SimpleDifferscaled_scaledQscat_ASCATNostripC20_nocorrection_withfirstmon_LinearTao','ERS2_date_final','ERS2_date_decimal','orginal_ERS2','scaled_ERS2')
ERS2_scaled_scaledQSCAT=scaled_ERS2'; clear scaled_ERS2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ASCAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('./mat/Americas_tropic_belt_Qscat_scaled_ASCAT_NoStripC20_nocorrection_withfirstMon_Linear.mat','A_data_final','A_date_final','A_date_decimal')

%% Figure ERS2 QSCAT ASCAT. All scaled against ASCAT. No correction applied
if 1
figure
set(gcf,'position',[248   226   950   450],...
                             'color','w','paperpositionmode','auto')

plot(ERS2_date_decimal,nanmean(ERS2_scaled_scaledQSCAT,2),'-*')


hold on
plot(Q_date_decimal,nanmean(Q_data_scaled_no_corrected,2),'-*')

hold on
a_mean=nanmean(A_data_final,2);
a_mean(a_mean>-4)=nan; % one outlier
plot(A_date_decimal,a_mean,'-*')
set(gca,'xtick',[1990 1993 1995 1997 2000 2003 2005 2007 2010 2013 2015 2017 2019])
legend({'ERS','QSCAT','ASCAT'},'location','southwest')

title('Americas-Before correction','FontSize',12,'FontWeight','bold')
xlabel('Year','FontSize',12,'FontWeight','bold')
ylabel('Radar signal (dB)','FontSize',12,'FontWeight','bold')


end

%% Ready to address the monthly differences between C- and Ku-band Radar data

QSCAT_final_comme_Cbands=nan(size(Q_data_scaled_no_corrected,2),size(Q_data_scaled_no_corrected,1));

good_count=0;
bad_count=0;
bad_corrected_count=0;

large_count=0;

R_thre = 0.3;

for coli=1:size(cell_locations,1)
%         disp(coli)

        QSCAT_pixel_data=Q_data_scaled_no_corrected(:,coli);         
        ERS2_pixel_data=ERS2_scaled_scaledQSCAT(:,coli);        
        ASCAT_pixel_data=A_data_final(:,coli); 
     
        
        [c,d]=ismember(ERS2_date_final,Q_date_final);
        [c2,d2]=ismember(Q_date_final,A_date_final);

        if sum(c)>6
                                 
            if (sum(isnan(ASCAT_pixel_data(d2(c2))))<sum(c2)/3) && (sum(isnan(QSCAT_pixel_data(c2)))<sum(c2)/3)
                
                large_count=large_count+1;
                
                %%%%%%%%%%%%% get E Q overlap  %%%%%%%%%%%%%%%%
                ERS2_data_benchmark= ERS2_pixel_data(c);               
                QSCAT_data_overlap_ERS2= QSCAT_pixel_data(d(c));
                
                %%%%%%%%%%%%% get Q A overlap %%%%%%%%%%%%%%%%
                QSCAT_data_overlap_ASCAT= QSCAT_pixel_data(c2);
                ASCAT_data_benchmark= ASCAT_pixel_data(d2(c2));

                %%%%%%%%%%%%%%%% combine the two overlapped periods. The variable name needs to be refined (but lazy) %%%%%%%%%%%
                QSCAT_tocdf=[QSCAT_data_overlap_ERS2;QSCAT_data_overlap_ASCAT];
                C_band_benchmark=[ERS2_data_benchmark;ASCAT_data_benchmark];
         
                %%%%%%%%%%%  Calculate r after removing NaN days %%%%%%%%%%%%%%%%%
                
                QSCAT_tocdf2=QSCAT_tocdf;
                C_band_benchmark2=C_band_benchmark;
                QSCAT_tocdf2(isnan(C_band_benchmark))=[];
                C_band_benchmark2(isnan(C_band_benchmark))=[];
                
                [corres2,p] = corrcoef(QSCAT_tocdf2,C_band_benchmark2);
                corre=corres2(1,2);
                p=p(1,2);
                
                clear  QSCAT_tocdf2  C_band_benchmark2
                

              %% counter           
                if  corre>=R_thre %p<0.05 &&                   
                    good_count=good_count+1;
                end
                
                if  corre<R_thre %p<0.05 &&                   
                    bad_count=bad_count+1;
                end
                                                       
               %% Address the monthly differences using TRMM, for pixels where C- and Ku- data have r < R_thre
                if corre<R_thre %corre<0.5 && p<0.05 
                    
                    AQ_differ=C_band_benchmark-QSCAT_tocdf; % after scaling
                   
                    Q_date_overlapped=ERS2_date_final(c);                   
                    [ccc,~]=ismember(TRMM_times,Q_date_overlapped); 
                    trmm_overlap=TRMM_allcell_allmonth(coli,ccc)';

                    Q_date_overlapped=Q_date_final(c2);                   
                    [ccc,~]=ismember(TRMM_times,Q_date_overlapped); 
                    trmm_overlap=[trmm_overlap;TRMM_allcell_allmonth(coli,ccc)'];
                    
                    %%%%%%%% linear regression correction
                    % stats: r2,F, p, error
                    [b,bint,r,rint,stats] = regress(AQ_differ,[trmm_overlap ones(length(trmm_overlap),1)]);
                    cor_slope=b(1);
                    cor_intercept=b(2);
                    QSCAT_tocdf_linear_correct=(trmm_overlap*cor_slope+cor_intercept)+QSCAT_tocdf;
                                        
                    %%%%%%% remove nan before corrcoef %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    C_band_benchmark2=C_band_benchmark;
                    SATCDF_ok_correct22=QSCAT_tocdf_linear_correct;
                    SATCDF_ok_correct22(isnan(C_band_benchmark))=[];
                    C_band_benchmark2(isnan(C_band_benchmark))=[];
                    
                    [corre_aftercorrection,p_aftercorrection] = corrcoef(C_band_benchmark2,SATCDF_ok_correct22);
                    RMSE_aftercorrection = sqrt(mean((C_band_benchmark2-SATCDF_ok_correct22).^2));
                    if corre_aftercorrection(1,2) > R_thre
                        bad_corrected_count=bad_corrected_count+1;
                                               
                        [cc2,dd2]=ismember(Q_date_final,TRMM_times);
                        QSCAT_all_data_Linear_corrected=(TRMM_allcell_allmonth(coli,dd2(cc2))*cor_slope+cor_intercept)'+QSCAT_pixel_data;
    %                     QSCAT_all_data_Linear_corrected=predict(tree,TRMM_allcell_allmonth(coli,dd2(cc2))')+QSCAT_pixel_data;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        QSCAT_final_comme_Cbands(coli,:)=QSCAT_all_data_Linear_corrected';
                        
                    else  %%%%% decision tree correction
                        tree = fitrtree(trmm_overlap,AQ_differ); 
                        QSCAT_tocdf_linear_correct = predict(tree,trmm_overlap)+QSCAT_tocdf;
                        
                        [cc2,dd2]=ismember(Q_date_final,TRMM_times);
%                         QSCAT_all_data_Linear_corrected=(TRMM_allcell_allmonth(coli,dd2(cc2))*cor_slope+cor_intercept)'+QSCAT_pixel_data;
                        QSCAT_all_data_Linear_corrected=predict(tree,TRMM_allcell_allmonth(coli,dd2(cc2))')+QSCAT_pixel_data;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        QSCAT_final_comme_Cbands(coli,:)=QSCAT_all_data_Linear_corrected';
                        

                        
                    end
                                               
                end
                
               %% Check for pixels where C- and Ku- data have r >= R_thre, whether r can be further improved 
                if corre>=R_thre %corre<0.5 && p<0.05 
                    
                    AQ_differ=C_band_benchmark-QSCAT_tocdf; 
                   
                    Q_date_overlapped=ERS2_date_final(c);                   
                    [ccc,~]=ismember(TRMM_times,Q_date_overlapped); 
                    trmm_overlap=TRMM_allcell_allmonth(coli,ccc)';

                    Q_date_overlapped=Q_date_final(c2);                   
                    [ccc,~]=ismember(TRMM_times,Q_date_overlapped); 
                    trmm_overlap=[trmm_overlap;TRMM_allcell_allmonth(coli,ccc)'];
                    
                    %%%%%%%% linear regression correction
                    % stats: r2, F, p, error
                    [b,bint,r,rint,stats] = regress(AQ_differ,[trmm_overlap ones(length(trmm_overlap),1)]);
                    cor_slope=b(1);
                    cor_intercept=b(2);
                    QSCAT_tocdf_linear_correct=(trmm_overlap*cor_slope+cor_intercept)+QSCAT_tocdf;
                    
                    
                    %%%%%%% remove nan before corrcoef %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    C_band_benchmark2=C_band_benchmark;
                    SATCDF_ok_correct22=QSCAT_tocdf_linear_correct;
                    SATCDF_ok_correct22(isnan(C_band_benchmark))=[];
                    C_band_benchmark2(isnan(C_band_benchmark))=[];
                    
                    [corre_aftercorrection,p_aftercorrection] = corrcoef(C_band_benchmark2,SATCDF_ok_correct22);
                    RMSE_aftercorrection = sqrt(mean((C_band_benchmark2-SATCDF_ok_correct22).^2));
                    
                    if corre_aftercorrection(1,2) > corre  %%% if after linear correction, r increased
                        
                        [cc2,dd2]=ismember(Q_date_final,TRMM_times);
                        QSCAT_all_data_Linear_corrected=(TRMM_allcell_allmonth(coli,dd2(cc2))*cor_slope+cor_intercept)'+QSCAT_pixel_data;
  
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        QSCAT_final_comme_Cbands(coli,:)=QSCAT_all_data_Linear_corrected';
                        
                    else  %%%%%% if after linear correction, r increased, use the original  file.

                        QSCAT_final_comme_Cbands(coli,:)=QSCAT_pixel_data';
                        
                    end
                                        
                    
                end
            end
    
        end

end

% save('./mat/Americas_tropic_belt_Final_time_series_LinearCDF_Noutlier_ERSSimpleDiff','ERS2_date_decimal','ERS2_date_final','ERS2_scaled_scaledQSCAT','Q_date_final','Q_date_decimal','QSCAT_final_comme_Cbands','A_date_final','A_date_decimal','A_data_final')

%% Figure ERS2 QSCAT ASCAT. All scaled against ASCAT. Correction applied
if 1

CT=cbrewer('qual', 'Paired', 12);

% CT_back = lbmap(10,'RedBlue'); %RedBlue BlueGray BrownBlue Blue
CT_back=cbrewer('seq', 'Greys', 12);

figure
set(gcf,'position',[248   226   950   450],...
                             'color','w','paperpositionmode','auto')
                         

plot(ERS2_date_decimal,nanmean(ERS2_scaled_scaledQSCAT,2),'-*','color',CT(2,:))


hold on
plot(Q_date_decimal,nanmean(QSCAT_final_comme_Cbands,1),'-*','color',CT(8,:))

hold on
a_mean=nanmean(A_data_final,2);
plot(A_date_decimal,a_mean,'-*','color',CT(4,:))

set(gca,'xtick',[1990 1993 1995 1997 2000 2003 2005 2007 2010 2013 2015 2017 2019])
legend({'ERS','QSCAT','ASCAT'},'location','southwest')

title('Americas-After correction','FontSize',12,'FontWeight','bold')
xlabel('Year','FontSize',12,'FontWeight','bold')
ylabel('Radar signal (dB)','FontSize',12,'FontWeight','bold')

end

%% Merge all sensor data into one continous time series at monthly time-step 
%%%% Data in the overlapped time period will be averaged across sensors.
if 1
   
    temp_time=num2str(ERS2_date_final(1));
    start_year=str2num(temp_time(1:4));
    start_mon=str2num(temp_time(5:6));
    temp_time=num2str(A_date_final(end));
    end_year=str2num(temp_time(1:4));
    end_mon=str2num(temp_time(5:6));
    
    %%%%% create a full time series from start date to end data at monthly time-step
    Final_time_monthly=[];
    Final_time_monthly_decimal=[];
    
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
            Final_time_monthly=[Final_time_monthly;temp_time];
            temp_time2=num2str(temp_time);
            temp_time_decimal=str2num(temp_time2(1:4))+str2num(temp_time2(5:6))/12;
            Final_time_monthly_decimal=[Final_time_monthly_decimal;temp_time_decimal];
        
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% full time series of three dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%  
    [c,d]=ismember(Q_date_final,Final_time_monthly);
    Final_time_series_QSCAT=nan(length(Final_time_monthly),size(cell_locations,1));
    Final_time_series_QSCAT(d(c),:)=QSCAT_final_comme_Cbands';
       
    [c,d]=ismember(ERS2_date_final,Final_time_monthly);
    Final_time_series_ERS=nan(length(Final_time_monthly),size(cell_locations,1));
    Final_time_series_ERS(d(c),:)=ERS2_scaled_scaledQSCAT;
    
    [c,d]=ismember(A_date_final,Final_time_monthly);
    Final_time_series_ASCAT=nan(length(Final_time_monthly),size(cell_locations,1));
    Final_time_series_ASCAT(d(c),:)=A_data_final;
    
    %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Final_time_series_FullRadar(:,:,1)=Final_time_series_ERS;
    Final_time_series_FullRadar(:,:,2)=Final_time_series_QSCAT;
    Final_time_series_FullRadar(:,:,3)=Final_time_series_ASCAT;
    
    
    Final_time_series_FullRadar=nanmean(Final_time_series_FullRadar,3);

    figure   
    plot(Final_time_monthly_decimal,nanmean(Final_time_series_FullRadar,2),'-','color',[0.5 0.5 0.5])

%     save('./mat/Americas_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff','Final_time_series_FullRadar','Final_time_monthly_decimal','Final_time_monthly','-v7.3')


end
