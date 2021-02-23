%% Project Radar_Rainforest
% Until 2021.02, there has been no well-calibrated C-band microwave (either active or passive) data set covering a long time span. 
% There have been passive microwave data (C-band included) from AMSR-E and AMSR2 covering a relatively long time span (since 2002), but the two didn't overlap temporally. 
% Hence, merging them with a full calibration has not been possible (Du et al. 2017. Earth Syst. Sci. Data 9,791–808; Moesinger et al. 2020. Earth Syst. Sci. Data 12, 177–196). 
% The aim of Project Radar_Rainforest is to fill this gap, by providing a well calibrated, long-term (since 1992) C-band radar data set for global land areas, especially for tropical rainforests.

%%
% This script will scale QSCAT against ASCAT data, per pixel (25km * 25km)
% The input (and output) files can be downloaded at:  https://filesender.renater.fr/?s=download&token=c0cb8523-76ab-4190-96f3-8cc260c2a67b
% Tropical Americas was used as an example here, but all input (and output) data for tropical Africa and Asia were also provided.

% Author: Shengli Tao
% Email: sltao1990@gmail.com
% Date: Created in 08.2018. Formatted in 02.2021.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

warning off
addpath('./util')



%% First, read in the pixel locations and IDs of each radar pixel 
%%% Here I included all tropical pixesl (~23N - 23S). Will subset to intact rainforest when drawing the figures.
load ./mat/Americas_tropical_belt_fishnet
cell_locations=[Am_cell_ID Am_cell_centers]; %%% Am_cell_ID: ID of each radar pixel;  Am_cell_centers: x y center of each radar pixel.  Each row is a pixel 

%% Load QSCAT and ASCAT data
load ./mat/Americas_tropical_belt_Q_data  %%%% This mat file contains the monthly QSCAT data for each radar pixel. 
%%% Q_data_final: monthly QSCAT data, 125 months, 22738 pixels; Q_date_final & Q_date_decimal: month labels.


load ./mat/Americas_tropical_belt_A_data_NoStripC20 %%%% This mat file contains the monthly ASCAT data for each radar pixel. Strips in ASCAT images have been removed by thresholding the number of observations.
%%% A_data_final: monthly ASCAT data, 141 months, 22738 pixels; A_date_final & A_date_decimal: month labels.

%% Quality control. Check outliers first. Very few pixels have outliers (visual inspection). Could be skipped
if 0

%%%%%%%%%%%%% ASCAT data are of high quality, but these two months seem to be outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_data_final(57,:)=(A_data_final(56,:)+A_data_final(58,:))/2; %%% replace these two months as averaged value of two adjacent months;
A_data_final(75,:)=(A_data_final(74,:)+A_data_final(76,:))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


outlier_count=0;
outlier_FID=[];
Q_data_final_noOutlier=[];

for i=1:size(Q_data_final,2)
    one_data=Q_data_final(:,i);
    
    %%%%%%%%%%%%%% hampel method, too loose %%%%%%%%%%%%%%%%%%%%%%%
%     [y,j] = hampel(one_data);
%     A_data_final2(:,i)=y;
%     if sum(j)>0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    outlier_ind=Isoutlier_3stdMean(one_data);
    if outlier_ind(13)==1 
        
        before_outlier_remove=one_data(outlier_ind);%%%% data before outlier interpolation, for figure purpose
        
        %%interpolate the outlier with adjacent values 

          left_neighbor=one_data(12);
          right_neighbor=one_data(14);

          one_data(13)=(left_neighbor+right_neighbor)/2;
              
  

        %%%%%%%% check the efficace of outlier interopolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            figure('visible','off')
%            plot(Q_date_decimal,one_data)
%            hold on
%            scatter(Q_date_decimal(outlier_ind),before_outlier_remove,30,'r','filled') %%%% old data
%            hold on
%            scatter(Q_date_decimal(outlier_ind),one_data(outlier_ind),20,'b','filled') %%%% new data
%            title(num2str(Q_date_final(outlier_ind)))          
%            saveas(gca,strcat('./jpg_Q_outlier_ecoregion60133/Pixel_',num2str(cell_locations(i,1)),'.jpg'))
%            close          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q_data_final_noOutlier(:,i)=one_data;

    else
        Q_data_final_noOutlier(:,i)=one_data;


    end
    
end

Q_data_final=Q_data_final_noOutlier;
clear Q_data_final_noOutlier

end

%% Ready to scale QSCAT data against ASCAT
if 1

orginal_QSCAT=nan(size(cell_locations,1),length(Q_date_final));
scaled_QSCAT=nan(size(cell_locations,1),length(Q_date_final));


good_count=0;
bad_count=0;

large_count=0;

R_thre=0.30;


for coli=1:size(cell_locations,1)

        QSCAT_pixel_data=Q_data_final(:,coli);           
        ASCAT_pixel_data=A_data_final(:,coli);  
        
       %%%%%%%%%%%% movemean didn't help much %%%%%%%%%%%%       
%         QSCAT_pixel_data=movmean(QSCAT_pixel_data,3,'omitnan');
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [c,d]=ismember(Q_date_final,A_date_final);

        if sum(c)>6 %%%%% make sure there are enough observations overlapped.

            if (sum(isnan(ASCAT_pixel_data(d(c))))<sum(c)/3) && (sum(isnan(QSCAT_pixel_data(c)))<sum(c)/3) %%%%% make sure there are enough observations overlapped.
                
                large_count=large_count+1;

                QSCAT_data_tocdf= QSCAT_pixel_data(c);
                ASCAT_data_benchmark= ASCAT_pixel_data(d(c));
              
               %% percentile cdf
                if 0
                percentiles=[0 30 50 70 100]; %10  90
%                 percentiles=[0 20 40 60 80 100]; %10  90
                SATmont_percentiles = prctile(QSCAT_data_tocdf,percentiles);
                [COEFFs,RMSE,corre,p] = QMAPP_piecewise(ASCAT_data_benchmark,QSCAT_data_tocdf,percentiles);

                SATCDF_ok=nan(length(QSCAT_data_tocdf),1);
                QSCAT_pixel_scaled = nan(length(QSCAT_pixel_data),1);


                for j=1:length(COEFFs)

                    ind1=(QSCAT_pixel_data>=SATmont_percentiles(j)) & (QSCAT_pixel_data<=SATmont_percentiles(j+1));
                    QSCAT_pixel_scaled(ind1)=QSCAT_pixel_data(ind1)*COEFFs(j,1)+COEFFs(j,2);

                    ind1=(QSCAT_data_tocdf>=SATmont_percentiles(j)) & (QSCAT_data_tocdf<=SATmont_percentiles(j+1));
                    SATCDF_ok(ind1)=QSCAT_data_tocdf(ind1)*COEFFs(j,1)+COEFFs(j,2);


                end

                ind1=(QSCAT_pixel_data<SATmont_percentiles(1));
                QSCAT_pixel_scaled(ind1)=QSCAT_pixel_data(ind1)*COEFFs(1,1)+COEFFs(1,2);

                ind1=(QSCAT_pixel_data>SATmont_percentiles(end));
                QSCAT_pixel_scaled(ind1)=QSCAT_pixel_data(ind1)*COEFFs(end,1)+COEFFs(end,2);

                
                end
                
               %%
                %%%%% linear scaling by Tao %%%%%%
                if 1
                    SATCDF_ok1 = scaledata(QSCAT_data_tocdf,nanmin(ASCAT_data_benchmark),nanmax(ASCAT_data_benchmark));
                    [COEFFs,RMSE,~,~] = QMAPP_linear(SATCDF_ok1,QSCAT_data_tocdf); 
                    
                    [corres,p] = corrcoef(SATCDF_ok1,ASCAT_data_benchmark);
                    corre=corres(1,2);
                    p=p(1,2);

                    SATCDF_ok=nan(length(QSCAT_data_tocdf),1);
                    QSCAT_pixel_scaled = nan(length(QSCAT_pixel_data),1);

                    SATCDF_ok=QSCAT_data_tocdf*COEFFs(1)+COEFFs(2);
                    QSCAT_pixel_scaled=QSCAT_pixel_data*COEFFs(1)+COEFFs(2);
                end
                
               %% add a mean difference as scale
                if 0
                    differ11=nanmean(QSCAT_data_tocdf-ASCAT_data_benchmark);
                    SATCDF_ok=QSCAT_data_tocdf-differ11;

                    %%%%%%%%%%% remove bad days, they are nan
                    SATCDF_ok2=SATCDF_ok;
                    ASCAT_data_benchmark2=ASCAT_data_benchmark;
                    SATCDF_ok2(isnan(SATCDF_ok))=[];
                    ASCAT_data_benchmark2(isnan(SATCDF_ok))=[];

                    [corres2,p] = corrcoef(ASCAT_data_benchmark2,SATCDF_ok2);
                    corre=corres2(1,2);
                    p=p(1,2);

                    clear  SATCDF_ok2  ASCAT_data_benchmark2

                    QSCAT_pixel_scaled=QSCAT_pixel_data-differ11;
                end
                               
               %% Some pixels (close to savana) have very low values <-10, need to check
                if sum(ASCAT_data_benchmark<-10)>0
%                     disp(cell_locations(coli,1))
%                     figure
%                     subplot(2,1,1)
%                     plot(ASCAT_data_benchmark)
%                     ylabel('ASCAT')
%                     subplot(2,1,2)
%                     plot(SATCDF_ok)
%                     ylabel('QSCAT')
                end    
               %% count good or bad
                if  corre>=R_thre %p<0.05 &&
                    good_count=good_count+1;
                end
                
                if  corre<R_thre %p<0.05 &&                   
                    bad_count=bad_count+1;
                end
                

               %%

                scaled_QSCAT(coli,:)=QSCAT_pixel_scaled';
                orginal_QSCAT(coli,:)=QSCAT_pixel_data';
                    
            end
        end


end

% save('./mat/Americas_tropic_belt_Qscat_scaled_ASCAT_NoStripC20_nocorrection_withfirstMon_Linear','Q_date_decimal','Q_date_final','orginal_QSCAT','scaled_QSCAT','A_data_final','A_date_final','A_date_decimal')

end

%% Check the results. Draw a continous time series at monthly level
if 1
   
figure

hold on
plot(Q_date_decimal,nanmean(scaled_QSCAT,1),'-*')

hold on
a_mean=nanmean(A_data_final,2);
plot(A_date_decimal,a_mean,'-*')

ylabel('Radar signal (dB)')

legend({'QSCAT','ASCAT'})


end
