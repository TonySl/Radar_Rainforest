%% Matlab code to scale QSCAT against ASCAT monthly radar observations for tropical Americas (25N-25S).

%%%    This code requires several mat files as input, which can be downloaded at figshare: xxxx.
%%%    Specially, 
%%%    1. We provided mat files for monthly radar data of ERS, QSCAT and ASCAT covering an extent from 25N to 25S. Raw radar images were downloaded from eumetsat and BYU.
%%%       We processed the images into a monthly time step and a resolution of 25km. 
%%%    2. The center latitudes, longitudes, ids, and area of all 25km sized pixels were also provided.
%%%    3. The land cover type (based on ESA CCI 2015 map), deforestation ratio (based on Hansen's maps), TRMM monthly precipitation, and CHIRPS monthly precipitation of each pixel were also calcualted and provided as mat files.

%%%    In the codes, all the mat files were put into a folder named as "mat"

%%% Author: Shengli Tao. 
%%% Email; sltao1990@gmail.com

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

warning off
addpath('D:\OneDrive\Drought')


load ./mat/Americas_tropical_belt_fishnet
cell_locations=[Am_cell_ID Am_cell_centers];


%% TRMM data
load ./mat/Americas_tropical_belt_TRMM_19982018

%% load Radar data
load ./mat/Americas_tropical_belt_Q_data
load ./mat/Americas_tropical_belt_A_data_NoStripC20

%% %%%%%%%%%%%%%%%%%%%%%%% Check outliers first %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%% Visual insepection found, ASCAT images in 201109 (row 57) and 201303 (row 75) have obvious outliers in Americas.
%%%%%%%%%%%%% Replace the values in these two months with averaged value of two adjacent months 
A_data_final(57,:)=(A_data_final(56,:)+A_data_final(58,:))/2; 
A_data_final(75,:)=(A_data_final(74,:)+A_data_final(76,:))/2;

%%%%%%%%%%%%%%%%%% QSCAT images in 2000.07 (row 13) has outliers in Americas, but not in every pixel
%%%%%%%%%%%%%%%%%% Use mean(3std) to find them
if 1

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
%            saveas(gca,strcat('./jpg_Q_outlier/Pixel_',num2str(cell_locations(i,1)),'.jpg'))
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

%% Ready to rescale QSCAT data against ASCAT
if 1

template1=ones(size(cell_locations,1),1);

orginal_QSCAT=nan(size(cell_locations,1),length(Q_date_final));
scaled_QSCAT=nan(size(cell_locations,1),length(Q_date_final));


Final_merged=[];
Final_RMSE=[];

good_count=0;
bad_count=0;


R_thre=0.30;


for coli=1:size(cell_locations,1)

        QSCAT_pixel_data=Q_data_final(:,coli);           
        ASCAT_pixel_data=A_data_final(:,coli);  
        
        
        [c,d]=ismember(Q_date_final,A_date_final);

        if sum(c)>6 %%%%% make sure there are enough observations in the overlaping period

            if (sum(isnan(ASCAT_pixel_data(d(c))))<sum(c)/3) && (sum(isnan(QSCAT_pixel_data(c)))<sum(c)/3)
                

                QSCAT_data_tocdf= QSCAT_pixel_data(c);
                ASCAT_data_benchmark= ASCAT_pixel_data(d(c));
                
                %%%%%%%%%%%%% Three scalling method tested below %%%%%%%%%%%%%
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
                %%%%% linear scaling %%%%%%
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
                                 
               %% count number of pixels with good or bad match quality
                if  corre>=R_thre %p<0.05 &&
                    good_count=good_count+1;
                end
                
                if  corre<R_thre %p<0.05 &&                   
                    bad_count=bad_count+1;
                end
                
              %%
                Final_merged(coli,:)=[QSCAT_pixel_scaled' ASCAT_pixel_data'];
                Final_RMSE(coli,:)=[corre p];

                scaled_QSCAT(coli,:)=QSCAT_pixel_scaled';
                orginal_QSCAT(coli,:)=QSCAT_pixel_data';
                    
            end
        end


end

save('./mat/Americas_tropic_belt_Qscat_scaled_ASCAT_NoStripC20_Linear','Q_date_decimal','Q_date_final','orginal_QSCAT','scaled_QSCAT','A_data_final','A_date_final','A_date_decimal')

end

%% draw a continous time series at monthly level
if 1
   
figure

hold on
plot(Q_date_decimal,nanmean(scaled_QSCAT,1),'-*')

hold on
a_mean=nanmean(A_data_final,2);
plot(A_date_decimal,a_mean,'-*')


end
