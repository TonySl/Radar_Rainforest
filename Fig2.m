%% 
%
% This script will create Fig. 2 in Tao et al. 
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

                         
%% Americas Pixel trends
axes1=axes('position',[.08 .62 .27-0.02 .3]);


load Americas_tropical_belt_fishnet
load Americas_tropical_belt_LUCC

xll=min(Am_cell_centers(Am_cell_area>0.062,1));
yll=min(Am_cell_centers(Am_cell_area>0.062,2));

xright=max(Am_cell_centers(Am_cell_area>0.062,1));
yup=max(Am_cell_centers(Am_cell_area>0.062,2));

nrows=int32((yup-yll)/0.25+1);
ncols=int32((xright-xll)/0.25+1);
cellsize=0.25;

R_map=[1/0.25 yup xll];


load DBEST_Americas_TrendMatlab_Tropical_belt %%%% DBEST trend of each pixel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Hansen ITR mask 
%load Americas_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<10 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal


%%%%% Vancutsem ITR mask
load Am_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
% size2=(0.25/0.00025)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal



all_scores=all_trend_matlab(ind_1);
Am_cell_centers=Am_cell_centers(ind_1,:);

display('decreasing ratio')
100*sum(all_scores<0)/length(all_scores)
sum(all_scores<0)*30*30/1000000 %%%in area, million km2
length(all_scores)*30*30/1000000 %%%in area, million km2


% dlmwrite('cell_center_americas.csv',Am_cell_centers,'delimiter',',')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% catogorize trend for figureing %%%%%%%%%%%%%%%%%%%%%%
P=nan(nrows,ncols);
[row, col] = setpostn(P, R_map, Am_cell_centers(:,2), Am_cell_centers(:,1));
row=double(nrows)-double(row);

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
% geotiffwrite('pixel_trend_radar.tif',flipud(P_new),R_map)
% geotiffwrite('pixel_trend_radar_original.tif',flipud(P),R_map)
 
% axesm('MapProjection','mercator','maplatlimit',[-19 15],'maplonlimit',[-100 172],...  %%%%  [-100 -33] mercator  eqdcylin 'maplatlimit',[-24 10],'maplonlimit',[-80 -44] [-50 50],'maplonlimit',[0 360]
% 'ParallelLabel','off','PlabelMeridian','west','MeridianLabel','off','MLabelParallel','south',...
% 'FontSize',6,'FontWeight','normal','PLineLocation',20,'MLineLocation',20);
 

latitudes = yll:cellsize:yup; 
latitudes=flipud(latitudes');  
longitudes = xll:cellsize:xright;
[latGrid, lonGrid] = meshgrat(latitudes,longitudes);
geoshow(latGrid,lonGrid,double(P_new),'DisplayType','texturemap');  %%%% texturemap


caxis([0 8]);
colormap(gca,colors);

load ('world_country_new','world_country_new')
S = geoshape(world_country_new);
hold on
geoshow(S, 'DefaultFaceColor', [1 1 1],'DefaultEdgeColor', 'k','Linewidth',0.5,'FaceAlpha',0)
xlim([xll+16 xright-10])
ylim([yll+11 yup-11])
box off
axis off
% tightmap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Legend
chandle = colorbar('Location','NorthOutside','FontSize',6,'FontWeight','normal'); % This line places the colorbar
set(get(chandle,'ylabel'),'String','Trend in radar signal ({10}^{-3} dB/yr)','FontSize',46,'FontWeight','normal'); % Set the colorbar's label

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

%'xticklabel',{'','10','20','30','40','50'},...

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
% xlabel('Trend in radar signal');
% set(gca, 'YAxisLocation', 'right')



%% %%%%%% Africa Pixel trends 
% figure
% set(gcf,'position',[46  46   950   600],...
%                              'color','w','paperpositionmode','auto')
                         
axes('position',[.36 .62 .27-0.02 .3]);

load Africa_tropical_belt_fishnet
load Africa_tropical_belt_LUCC

xll=min(Af_cell_centers(Af_cell_area>0.062,1));
yll=min(Af_cell_centers(Af_cell_area>0.062,2));

xright=max(Af_cell_centers(Af_cell_area>0.062,1));
yup=max(Af_cell_centers(Af_cell_area>0.062,2));

nrows=int32((yup-yll)/0.25+1);
ncols=int32((xright-xll)/0.25+1);
cellsize=0.25;

R_map=[1/0.25 yup xll];


load DBEST_Africa_TrendMatlab_Tropical_belt

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
sum(all_scores<0)*30*30/1000000
length(all_scores)*30*30/1000000
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
xlim([xll+17 xright-17])
ylim([yll+15 yup-17])
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
% xlabel('Trend in radar signal');
% set(gca, 'YAxisLocation', 'right')


%% %%%%%% Asia Pixel trends 

% figure
% set(gcf,'position',[46  46   950   600],...
%                              'color','w','paperpositionmode','auto')


axes('position',[.7 .62 .27-0.02 .3]);

load Asia_tropical_belt_fishnet
load Asia_tropical_belt_LUCC


xll=min(Asia_cell_centers(Asia_cell_area>0.062,1));
yll=min(Asia_cell_centers(Asia_cell_area>0.062,2));

xright=max(Asia_cell_centers(Asia_cell_area>0.062,1));
yup=max(Asia_cell_centers(Asia_cell_area>0.062,2));

nrows=int32((yup-yll)/0.25+1);
ncols=int32((xright-xll)/0.25+1);
cellsize=0.25;

R_map=[1/0.25 yup xll];


load DBEST_Asia_TrendMatlab_Tropical_belt

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
sum(all_scores<0)*30*30/1000000
length(all_scores)*30*30/1000000

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
xlim([xll+24 xright-11])
ylim([yll+15 yup-12])
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
% xlabel('Trend in radar signal');
% set(gca, 'YAxisLocation', 'right')


%%
% print(gcf,'-dtiff','-r300',strcat('./Fig2.tif')) 

