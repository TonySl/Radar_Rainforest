%% 
% This script will create Fig. 2 in Tao et al. 


% Author: Shengli Tao, assisstant professor in Peking University
% Email: sltao@pku.edu.cn
% Date: First version in 08.2018. Formatted in May.2022.


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

%%%%% Color set by Wenwen %%%%%
green_wenwen=[155 213 216;95 197 202;31 159 164; 4 105 109]/255;
blue_wenwen=[186 220 249;121 189 243;3 118 206;0 85 153]/255; 
warm_wenwen=[253 196 178;249 161 130; 240 114 69;226 65 8; 170 47 3]/255;
yellow_wenwen=[248 192 17]/255;

colors=[1 1 1; flipud(warm_wenwen);green_wenwen(1,:);green_wenwen(2,:)]; 



%% 
% figure
% set(gcf,'position',[46  46   950   600],...
%                              'color','w','paperpositionmode','auto')

                         
%% Americas Pixel trends
figure
set(gcf,'position',[46  46   950   600],...
                             'color','w','paperpositionmode','auto')
                         
axes1=axes('position',[.08 .62 .27-0.02 .3]);


load Americas_tropical_belt_fishnet
load Americas_tropical_belt_LUCC


%%%%%% define mask coordinate %%%%%%
xll=-85.091241018;
yll=-16.2614185455;

xright=-46.091241018;
yup=15.4885814545;

nrows=127;
ncols=156;
cellsize=0.25;

R_map=[1/0.25 yup xll];


%%%%% radar signal trend calculation %%%%%
load Americas_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff 
all_trend_matlab=size(Final_time_series_FullRadar,2);
for i=1:size(Final_time_series_FullRadar,2)
    signal_onepixel=Final_time_series_FullRadar(:,i);
    ind_nan=(~isnan(signal_onepixel));
    P = polyfit(Final_time_monthly_decimal(ind_nan),signal_onepixel(ind_nan),1); 
    slope = P(1);
    all_trend_matlab(i)=slope;   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Hansen ITR mask 
%load Americas_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal


%%%%% Vancutsem ITR mask
load Am_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal



all_scores=all_trend_matlab(ind_1);
Am_cell_centers=Am_cell_centers(ind_1,:);



display('decreasing ratio')
100*sum(all_scores<0)/length(all_scores)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% catogorize trend for figureing %%%%%%%%%%%%%%%%%%%%%%
P=nan(nrows,ncols);
[row, col] = setpostn(P, R_map, Am_cell_centers(:,2), Am_cell_centers(:,1));
row=double(nrows)-double(row)+1;

P(sub2ind(size(P),row,col))=all_scores;

P_new=nan(size(P));
P_new(isnan(P))=0;
P_new(P<-1.3*0.001*12)=1;
P_new(P>=-1.3*0.001*12 & P<-1*0.001*12)=2;
P_new(P>=-1*0.001*12 & P<-0.6*0.001*12)=3;
P_new(P>=-0.6*0.001*12 & P<-0.3*0.001*12)=4;
P_new(P>=-0.3*0.001*12 & P<0*12)=5;
P_new(P>=0 & P<0.3*0.001*12)=6;
P_new(P>=0.3*0.001*12)=7;

P_new_america=P_new;

latitudes = yll:cellsize:yup; 
latitudes=flipud(latitudes');  
longitudes = xll:cellsize:xright;
[latGrid, lonGrid] = meshgrat(latitudes,longitudes);
geoshow(latGrid,lonGrid,double(P_new),'DisplayType','texturemap');   


caxis([0 8]);
colormap(gca,colors);

load ('world_country_new','world_country_new')
S = geoshape(world_country_new);
hold on
geoshow(S, 'DefaultFaceColor', [1 1 1],'DefaultEdgeColor', 'k','Linewidth',0.5,'FaceAlpha',0)
xlim([xll-9 xright+1])
ylim([yll-1 yup+1])
box off
axis off
% tightmap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Legend
chandle = colorbar('Location','NorthOutside','FontSize',6,'FontWeight','normal');  
set(get(chandle,'ylabel'),'String','Trend in radar signal ({10}^{-3} dB/yr)','FontSize',46,'FontWeight','normal');  

set(chandle, 'ylim', [1 8])


lab = [{'-16'};{'-12'};{'-8'};{'-4'};{'0'};{'4'}];
set(chandle,'YTick',[2:1:7],'yticklabel',lab,'fontsize',9);
set(chandle, 'TickLength', [0 0]);

x=get(chandle,'Position');
x(1)=0.31;
x(2)=0.61;
x(3)=0.37;
x(4)=.016;
set(chandle,'Position',x)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% histogram inset
axes('position',[.09 .68 .05 .11],'XColor', [1,1,1], 'YColor', [0,0,0])
hold on 

P_new_america=reshape(P_new_america,size(P_new_america,1)*size(P_new_america,2),1);
P_new_america(P_new_america==0)=[];

x = 1:7;
for i = x
hbar = bar(i,sum(P_new_america==i)/length(P_new_america)*100,...
    'facecolor',colors(i+1,:),'barwidth',1,...
    'edgecolor','none');
end


set(gca,'xlim',[0.5 7.5],...
    'ylim',[0 45],...
    'fontname','arial',...
    'fontsize',8,...
    'xtick', 1.5:1:7.5,...
    'xticklabel',{},...
    'ytick',0:20:40,...
    'yticklabel',{'','20%','40%'},...
    'tickdir','out',...
    'color','none')
set(gca, 'TickLength', [0.03 0.03]);


% print(gcf,'-dtiff','-r300',strcat('./Fig2_americas_vancutsem95_v2.tif')) 

%% %%%%%% Africa Pixel trends 
figure
set(gcf,'position',[46  46   950   600],...
                             'color','w','paperpositionmode','auto')
                         
axes('position',[.36 .62 .27-0.02 .3]);

load Africa_tropical_belt_fishnet
load Africa_tropical_belt_LUCC


xll=-9.09128937231;
yll=-5.26141318182;

xright=29.9087106277;
yup=6.23858681818;

nrows=46;
ncols=156;
cellsize=0.25;

R_map=[1/0.25 yup xll];



load Africa_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff 
all_trend_matlab=size(Final_time_series_FullRadar,2);
for i=1:size(Final_time_series_FullRadar,2)
    signal_onepixel=Final_time_series_FullRadar(:,i);
    ind_nan=(~isnan(signal_onepixel));
    P = polyfit(Final_time_monthly_decimal(ind_nan),signal_onepixel(ind_nan),1); 
    slope = P(1);
    all_trend_matlab(i)=slope;   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Hansen ITR mask 
% load Africa_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal


%%%%% Vancutsem ITR mask 
load Af_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Af_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal



all_scores=all_trend_matlab(ind_1);
Af_cell_centers=Af_cell_centers(ind_1,:);



display('decreasing ratio')
100*sum(all_scores<0)/length(all_scores)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% catogorize the trends %%%%%%%%%%%%%%%%%%%%%%
P=nan(nrows,ncols);
[row, col] = setpostn(P, R_map, Af_cell_centers(:,2), Af_cell_centers(:,1));
row=double(nrows)-double(row)+1;

P(sub2ind(size(P),row,col))=all_scores;


P_new=nan(size(P));
P_new(isnan(P))=0;

P_new(P<-1.3*0.001*12)=1;
P_new(P>=-1.3*0.001*12 & P<-1*0.001*12)=2;
P_new(P>=-1*0.001*12 & P<-0.6*0.001*12)=3;
P_new(P>=-0.6*0.001*12 & P<-0.3*0.001*12)=4;
P_new(P>=-0.3*0.001*12 & P<0*12)=5;
P_new(P>=0 & P<0.3*0.001*12)=6;
P_new(P>=0.3*0.001*12)=7;

P_new_africa=P_new;


latitudes = yll:cellsize:yup; 
latitudes=flipud(latitudes');  
longitudes = xll:cellsize:xright;
[latGrid, lonGrid] = meshgrat(latitudes,longitudes);
geoshow(latGrid,lonGrid,double(P_new),'DisplayType','texturemap');

caxis([0 8]);  

colormap(gca,colors);



load('world_country_new','world_country_new')
S = geoshape(world_country_new);
hold on
geoshow(S, 'DefaultFaceColor', [1 1 1],'DefaultEdgeColor', 'k','Linewidth',0.5,'FaceAlpha',0)
xlim([xll+9 xright+3])
ylim([yll-3 yup+3])
hold on
% tightmap
box off
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% histogram inset
axes('position',[.37 .68 .05 .11],'XColor', [1,1,1], 'YColor', [0,0,0])
hold on 

P_new_africa=reshape(P_new_africa,size(P_new_africa,1)*size(P_new_africa,2),1);
P_new_africa(P_new_africa==0)=[];

x = 1:7;
for i = x
hbar = bar(i,sum(P_new_africa==i)/length(P_new_africa)*100,...
    'facecolor',colors(i+1,:),'barwidth',1,...
    'edgecolor','none');
end


set(gca,'xlim',[0.5 7.5],...
    'ylim',[0 45],...
    'fontname','arial',...
    'fontsize',8,...
    'xtick', 1.5:1:7.5,...
    'xticklabel',{},...
    'ytick',0:20:40,...
    'yticklabel',{'','20%','40%'},...
    'tickdir','out',...
    'color','none')
set(gca, 'TickLength', [0.03 0.03]);

% print(gcf,'-dtiff','-r300',strcat('./Fig2_africa_vancutsem95_v2.tif')) 

%% %%%%%% Asia Pixel trends 

figure
set(gcf,'position',[46  46   950   600],...
                             'color','w','paperpositionmode','auto')


axes('position',[.7 .62 .27-0.02 .3]);

load Asia_tropical_belt_fishnet
load Asia_tropical_belt_LUCC


xll=95.4087106277;
yll=-10.7614588941;

xright=162.158710628;
yup=5.23854110592;

nrows=64;
ncols=267;
cellsize=0.25;

R_map=[1/0.25 yup xll];



load Asia_tropic_belt_Final_time_series_FullRadar_LinearCDF_Noutlier_ERSSimpleDiff 
all_trend_matlab=size(Final_time_series_FullRadar,2);
for i=1:size(Final_time_series_FullRadar,2)
    signal_onepixel=Final_time_series_FullRadar(:,i);
    ind_nan=(~isnan(signal_onepixel));
    P = polyfit(Final_time_monthly_decimal(ind_nan),signal_onepixel(ind_nan),1); 
    slope = P(1);
    all_trend_matlab(i)=slope;   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Hansen ITR mask 
% load Asia_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal


%%%%% Vancutsem ITR mask 
load Asia_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
% size2=(0.25/0.00025)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Asia_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal



all_scores=all_trend_matlab(ind_1);
Asia_cell_centers=Asia_cell_centers(ind_1,:);



display('decreasing ratio')
100*sum(all_scores<0)/length(all_scores)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% catogorize the trends %%%%%%%%%%%%%%%%%%%%%%
P=nan(nrows,ncols);
[row, col] = setpostn(P, R_map, Asia_cell_centers(:,2), Asia_cell_centers(:,1));
row=double(nrows)-double(row)+1;

P(sub2ind(size(P),row,col))=all_scores;


P_new=nan(size(P));
P_new(isnan(P))=0;

P_new(P<-1.3*0.001*12)=1;
P_new(P>=-1.3*0.001*12 & P<-1*0.001*12)=2;
P_new(P>=-1*0.001*12 & P<-0.6*0.001*12)=3;
P_new(P>=-0.6*0.001*12 & P<-0.3*0.001*12)=4;
P_new(P>=-0.3*0.001*12 & P<0*12)=5;
P_new(P>=0 & P<0.3*0.001*12)=6;
P_new(P>=0.3*0.001*12)=7;

P_new_asia=P_new;



latitudes = yll:cellsize:yup;
latitudes=flipud(latitudes'); 
longitudes = xll:cellsize:xright;
[latGrid, lonGrid] = meshgrat(latitudes,longitudes);
geoshow(latGrid,lonGrid,double(P_new),'DisplayType','texturemap');

caxis([0 8]); 

colormap(gca,colors);


load('world_country_new','world_country_new')
S = geoshape(world_country_new);
hold on
geoshow(S, 'DefaultFaceColor', [1 1 1],'DefaultEdgeColor', 'k','Linewidth',0.5,'FaceAlpha',0)
xlim([xll-1 xright+1])
ylim([yll-1 yup+2])
hold off 
% tightmap
box off
axis off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% histogram inset
axes('position',[.67 .68 .05 .11],'XColor', [1,1,1], 'YColor', [0,0,0])
hold on 

P_new_asia=reshape(P_new_asia,size(P_new_asia,1)*size(P_new_asia,2),1);
P_new_asia(P_new_asia==0)=[];

x = 1:7;
for i = x
hbar = bar(i,sum(P_new_asia==i)/length(P_new_asia)*100,...
    'facecolor',colors(i+1,:),'barwidth',1,...
    'edgecolor','none');
end

set(gca,'xlim',[0.5 7.5],...
    'ylim',[0 45],...
    'fontname','arial',...
    'fontsize',8,...
    'xtick', 1.5:1:7.5,...
    'xticklabel',{},...
    'ytick',0:20:40,...
    'yticklabel',{'','20%','40%'},...
    'tickdir','out',...
    'color','none')
set(gca, 'TickLength', [0.03 0.03]);

% print(gcf,'-dtiff','-r300',strcat('./Fig2_asia_vancutsem95_v2.tif')) 
