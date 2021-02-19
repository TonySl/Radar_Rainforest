
%% 
%
% This script will scale ERS against QSCAT data, per pixel (25km * 25km)
% The input (and output) files can be downloaded at: xxxxx
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

%% Load ERS data

load ./mat/Americas_tropical_belt_ERS1_data %%%% This mat file contains the monthly ERS1 data for each radar pixel. 
%%% ERS1_data_final: monthly ERS1 data, 36 months, 22738 pixels; ERS1_date_final & ERS1_date_decimal: month labels.

load ./mat/Americas_tropical_belt_ERS2_data %%%% This mat file contains the monthly ERS2 data for each radar pixel. 
%%% ERS2_data_final: monthly ERS2 data, 48 months, 22738 pixels; ERS2_date_final & ERS2_date_decimal: month labels.

%% Load scaled QSCAT data

load ./mat/Americas_tropic_belt_Qscat_scaled_ASCAT_NoStripC20_nocorrection_withfirstMon_Linear %%%% This mat file contains the new monthly QSCAT data that 
%%% have been scaled against ASCAT previously

Q_data_scaled_no_corrected=scaled_QSCAT';

%%%%%%%%%%%%%%%%%%%%%% Combine ERS 1 & 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ERS2_data_final=[ERS1_data_final;ERS2_data_final]; % A better variable name should be ERS_data_final...
ERS2_date_decimal=[ERS1_date_decimal; ERS2_date_decimal];
ERS2_date_final=[ERS1_date_final;ERS2_date_final];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Quality control. Check outliers first. Very few pixels have outliers (visual inspection). Could be skipped
%%%% 1997.06 (row 49), north amazonia is contaminated. Visually confirmed %%%%%%% 
if 0

% figure
% histogram(ERS2_data_final(49,:))
% hold on
% histogram(ERS2_data_final(48,:))

%%%% how many pixels have 1997.06 as outlier?
outlier_count=0;
outlier_FID=[];
ERS2_data_final_noOutlier=[];

for i=1:size(ERS2_data_final,2)
    one_data=ERS2_data_final(:,i);
    
    %%%%%%%%%%%%%% hampel method, too loose %%%%%%%%%%%%%%%%%%%%%%%
%     [y,j] = hampel(one_data);
%     A_data_final2(:,i)=y;
%     if sum(j)>0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    outlier_ind=Isoutlier_3stdMean_America_ERS199706(one_data);
    if outlier_ind(49)==1
        outlier_count=outlier_count+1;
        %%%% in which ecoregion %%%%%
        outlier_FID=[outlier_FID; cell_locations(i,1)];
        
%        figure('visible','off')
%        plot(ERS2_date_decimal,one_data)
%        hold on
%        scatter(ERS2_date_decimal(outlier_ind),one_data(outlier_ind),20,'r','filled') %%%% new data
%        title(num2str(ERS2_date_final(outlier_ind)))          
%        saveas(gca,strcat('./jpg_ERS2_outlier_3stdMean_199706/Pixel_',num2str(cell_locations(i,1)),'.jpg'))
%        close    
        
        
        %%%%%%%%%% correct 1997.06 by reducing 10% %%%%%%%%%%
        one_data(49)=one_data(49)+one_data(49)*0.1;
        ERS2_data_final_noOutlier(:,i)=one_data;

    else
        ERS2_data_final_noOutlier(:,i)=one_data;


    end
    
end

ERS2_data_final=ERS2_data_final_noOutlier;
clear ERS2_data_final_noOutlier



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outlier_count=0;
outlier_time=[];

ERS2_data_final_noOutlier=[];

for i=1:size(ERS2_data_final,2)
    one_data=ERS2_data_final(:,i);
    
    %%%%%%%%%%%%%% hampel method, too loose %%%%%%%%%%%%%%%%%%%%%%%
%     [y,j] = hampel(one_data);
%     A_data_final2(:,i)=y;
%     if sum(j)>0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    outlier_ind=Isoutlier_4stdMean_America_ERS(one_data);
    if sum(outlier_ind)>0
           outlier_count=outlier_count+1;
           outlier_time=[outlier_time;ERS2_date_final(outlier_ind)];
           
           before_outlier_remove=one_data(outlier_ind);%%%% data before outlier interpolation, for figure purpose
           
           %%interpolate the outlier with adjacent values 
           outlier_id=find(outlier_ind==1);
    
           for j2=1:length(outlier_id)
               
               if outlier_id(j2)==1 % if first value is outlier
                  left_neighbor=one_data(outlier_id(j2)+1);
                  right_neighbor=one_data(outlier_id(j2)+1);
               elseif outlier_id(j2)==length(one_data) % if last value is outlier
                  left_neighbor=one_data(outlier_id(j2)-1);
                  right_neighbor=one_data(outlier_id(j2)-1);
               else
                  left_neighbor=one_data(outlier_id(j2)-1);
                  right_neighbor=one_data(outlier_id(j2)+1);
               end    
                 
               
               if sum(isnan([left_neighbor right_neighbor]))>0
                   disp('neighbors nan.......')
               end              
               one_data(outlier_id(j2))=(left_neighbor+right_neighbor)/2;           
           end   
            
           %%%%%%%% check the efficace of outlier interopolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            figure('visible','off')
%            plot(ERS2_date_decimal,one_data)
%            hold on
%            scatter(ERS2_date_decimal(outlier_ind),before_outlier_remove,30,'r','filled') %%%% old data
%            hold on
%            scatter(ERS2_date_decimal(outlier_ind),one_data(outlier_ind),20,'r','filled') %%%% new data
%            title(num2str(ERS2_date_final(outlier_ind)))          
%            saveas(gca,strcat('./jpg_ERS2_outlier_4stdMean/Pixel_',num2str(cell_locations(i,1)),'.jpg'))
%            close          
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           ERS2_data_final_noOutlier(:,i)=one_data;
    
    else % if no outlier
        
        ERS2_data_final_noOutlier(:,i)=one_data;
        
    end


end

% disp(unique(outlier_time))

ERS2_data_final=ERS2_data_final_noOutlier;
clear ERS2_data_final_noOutlier

end

%% Ready to scale ERS against scaled_QSCAT


orginal_ERS2=nan(size(cell_locations,1),length(ERS2_date_final));
scaled_ERS2=nan(size(cell_locations,1),length(ERS2_date_final));


good_count=0;
bad_count=0;

large_count=0;
R_thre=0.30;


for coli=1:size(cell_locations,1)
%         disp(coli)

        QSCAT_pixel_data=Q_data_scaled_no_corrected(:,coli);
        ERS2_pixel_data=ERS2_data_final(:,coli); 
        

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [c,d]=ismember(ERS2_date_final,Q_date_final);

        if sum(c)>6
                       
          
            if (sum(isnan(QSCAT_pixel_data(d(c))))<sum(c)/3) && (sum(isnan(ERS2_pixel_data(c)))<sum(c)/3)
                
                large_count=large_count+1;
                
                ERS_data_tocdf= ERS2_pixel_data(c);
                QSCAT_data_benchmark= QSCAT_pixel_data(d(c));
                
              %% Add a mean difference as scale
                if 1
                differ11=nanmean(ERS_data_tocdf-QSCAT_data_benchmark);
                SATCDF_ok=ERS_data_tocdf-differ11;
  
                %%%%%%%%%%% remove nan
                SATCDF_ok2=SATCDF_ok;
                QSCAT_data_benchmark2=QSCAT_data_benchmark;
                SATCDF_ok2(isnan(SATCDF_ok))=[];
                QSCAT_data_benchmark2(isnan(SATCDF_ok))=[];
                
                [corres2,p] = corrcoef(QSCAT_data_benchmark2,SATCDF_ok2);
                corre=corres2(1,2);
                p=p(1,2);
                
                clear  SATCDF_ok2  QSCAT_data_benchmark2
                               
                ERS2_pixel_scaled=ERS2_pixel_data-differ11;
                end           
              %% linear cdf
                if 0
                    SATCDF_ok1 = scaledata(ERS_data_tocdf,nanmin(QSCAT_data_benchmark),nanmax(QSCAT_data_benchmark));
                    [COEFFs,RMSE,~,~] = QMAPP_linear(SATCDF_ok1,ERS_data_tocdf); 
                    
                    [corres,p] = corrcoef(SATCDF_ok1,QSCAT_data_benchmark);
                    corre=corres(1,2);
                    p=p(1,2);

                    SATCDF_ok=nan(length(ERS_data_tocdf),1);
                    ERS2_pixel_scaled = nan(length(ERS2_pixel_data),1);

                    SATCDF_ok=ERS_data_tocdf*COEFFs(1)+COEFFs(2);
                    ERS2_pixel_scaled=ERS2_pixel_data*COEFFs(1)+COEFFs(2);
                end
              %% percentile cdf
                if 0
                percentiles=[0 30 50 70 100]; %10  90
%                 percentiles=[0 20 40 60 80 100]; %10  90
                SATmont_percentiles = prctile(ERS_data_tocdf,percentiles);
                [COEFFs,RMSE,corre,p] = QMAPP_piecewise(QSCAT_data_benchmark,ERS_data_tocdf,percentiles);

                SATCDF_ok=nan(length(ERS_data_tocdf),1);
                ERS2_pixel_scaled = nan(length(ERS2_pixel_data),1);

                for j=1:length(COEFFs)

                    ind1=(ERS2_pixel_data>=SATmont_percentiles(j)) & (ERS2_pixel_data<=SATmont_percentiles(j+1));
                    ERS2_pixel_scaled(ind1)=ERS2_pixel_data(ind1)*COEFFs(j,1)+COEFFs(j,2);

                    ind1=(ERS_data_tocdf>=SATmont_percentiles(j)) & (ERS_data_tocdf<=SATmont_percentiles(j+1));
                    SATCDF_ok(ind1)=ERS_data_tocdf(ind1)*COEFFs(j,1)+COEFFs(j,2);


                end

                ind1=(ERS2_pixel_data<SATmont_percentiles(1));
                ERS2_pixel_scaled(ind1)=ERS2_pixel_data(ind1)*COEFFs(1,1)+COEFFs(1,2);

                ind1=(ERS2_pixel_data>SATmont_percentiles(end));
                ERS2_pixel_scaled(ind1)=ERS2_pixel_data(ind1)*COEFFs(end,1)+COEFFs(end,2);
                         
                
                %%%%%%%%%%% remove bad days, they are nan
                SATCDF_ok2=SATCDF_ok;
                QSCAT_data_benchmark2=QSCAT_data_benchmark;
                SATCDF_ok2(isnan(SATCDF_ok))=[];
                QSCAT_data_benchmark2(isnan(SATCDF_ok))=[];
                
                [corres2,p] = corrcoef(QSCAT_data_benchmark2,SATCDF_ok2);
                corre=corres2(1,2);
                p=p(1,2);
                
                clear  SATCDF_ok2  QSCAT_data_benchmark2
                

                end
              %%                  
                if  corre>=R_thre %p<0.05 &&                   
                    good_count=good_count+1;
                end
                
                if  corre<R_thre %p<0.05 &&                   
                    bad_count=bad_count+1;
                end
                             

                scaled_ERS2(coli,:)=ERS2_pixel_scaled';
                orginal_ERS2(coli,:)=ERS2_pixel_data';
                    
                                   
            end
            
      
        end

end

% save('./mat/Americas_tropic_belt_ERS2_SimpleDifferscaled_scaledQscat_ASCATNostripC20_nocorrection_withfirstmon_LinearTao','ERS2_date_final','ERS2_date_decimal','orginal_ERS2','scaled_ERS2')

%% Figure ERS2 QSCAT ASCAT. All scaled against ASCAT. No correction applied
if 1
figure
set(gcf,'position',[248   226   1450   450],...
                             'color','w','paperpositionmode','auto')

plot(ERS2_date_decimal,nanmean(scaled_ERS2,1),'-*')

hold on
plot(Q_date_decimal,nanmean(scaled_QSCAT,1),'-*')

hold on
a_mean=nanmean(A_data_final,2);
a_mean(a_mean>-4)=nan; % catch one outlier here :)
plot(A_date_decimal,a_mean,'-*')
set(gca,'xtick',[1990 1993 1995 1997 2000 2003 2005 2007 2010 2013 2015 2017 2019])
legend({'ERS','QSCAT','ASCAT'},'location','southwest')

title('Americas','FontSize',12,'FontWeight','bold')
xlabel('Year','FontSize',12,'FontWeight','bold')
ylabel('Radar signal (dB)','FontSize',12,'FontWeight','bold')
% xlim([1998 2002])

% print(gcf,'-dtiff','-r300',strcat('./Americas_90_19_Radar.tif'))

end
