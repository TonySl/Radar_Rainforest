%% 
%
% This script will create Fig. 3 in Tao et al. 
% All input (and output) files can be downloaded at: https://doi.org/10.6084/m9.figshare.14061428.v3

% Author: Shengli Tao
% Email: sltao1990@gmail.com, sltao@pku.edu.cn
% Date: First version in 08.2018. Formatted in 01.2022.

clc
clear

addpath('./util')
addpath('./mat')
warning off

scatter_alpha=0.3;
scatter_size=10;


%%
figure
set(gcf,'position',[46  46   950   550],...
                             'color','w','paperpositionmode','auto')
   

%% %%%%%% Americas  
pos2 = [0.05+0.08 0.69+0.0865 0.23 0.24];

axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

load AR1_America_lineardetrend_pixellevel
load Americas_tropical_belt_LUCC
load Americas_tropical_belt_fishnet


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subset to intact tropical rainforest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Hansen ITR mask % load Americas_tropic_belt_HansenCover90_non_forest_count
% ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & All_tile_allcell_non_forest_counter_ratio<5 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal


%%%%% Vancutsem ITR mask
load Am_intact_vancutsem_count
size2=(0.25/0.000269495)^2;
allcell_intact_ratio=100*cell2mat(allcell_intact_count)/size2;
ind_1=cell_dominantLUCC_pixelratio>50 & ismember(cell_dominantLUCC,[50 160 170]) & Am_cell_area>0.062 & allcell_intact_ratio>=95 & cell_water_pixelratio<5  & cell_dominantLUCC~=200; %%% derset removal
disp(sum(ind_1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Final_time_series_FullRadar=all_scores(ind_1,:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Final_time_series_FullRadar_allpixel_average=nanmean(Final_time_series_FullRadar,1);

ax = gca;
ax.YColor = 'k';

load Americas_tropical_belt_Ecoregion
Olson_ID=cell_dominantEcoregion;
ecoregion_id_americas{1,1}=[60142;  60163;  60107; 60132; 60158; 60143]; 
ecoregion_id_americas{2,1}=[60125; 60182; 60124;60169;60173]; 
ecoregion_id_americas{3,1}=[60135; 60168; 60180; 60170; 60157; 60140; 60180; 60170; 60138; 60212]; 
ecoregion_id_americas{4,1}=[60166;  60133; 60157];

Olson_ID=Olson_ID(ind_1);

ecoregion1_pixelAR1_ind=ismember(Olson_ID,ecoregion_id_americas{1,1});
ecoregion2_pixelAR1_ind=ismember(Olson_ID,ecoregion_id_americas{2,1});
ecoregion3_pixelAR1_ind=ismember(Olson_ID,ecoregion_id_americas{3,1});
ecoregion4_pixelAR1_ind=ismember(Olson_ID,ecoregion_id_americas{4,1});
ecoregion1_pixelAR1=Final_time_series_FullRadar(ecoregion1_pixelAR1_ind,:);
ecoregion2_pixelAR1=Final_time_series_FullRadar(ecoregion2_pixelAR1_ind,:);
ecoregion3_pixelAR1=Final_time_series_FullRadar(ecoregion3_pixelAR1_ind,:);
ecoregion4_pixelAR1=Final_time_series_FullRadar(ecoregion4_pixelAR1_ind,:);

ecoregions_pixelAR1=[nanmean(ecoregion1_pixelAR1,1);nanmean(ecoregion2_pixelAR1,1);nanmean(ecoregion3_pixelAR1,1);nanmean(ecoregion4_pixelAR1,1)];

reds = cbrewer('seq','OrRd',10);
colors=reds([3 5 8 10],:); 
line_colors=colors;


taub_allecoregions=[];

for i=[1 4 2 3] 
      
        s1=scatter(ar_time,ecoregions_pixelAR1(i,:),scatter_size,'MarkerFaceColor',line_colors(i,:),'MarkerEdgeColor','none'); %%%% ,symbols{i} ,'filled'
        s1.MarkerFaceAlpha = scatter_alpha;
        s1.MarkerEdgeAlpha = scatter_alpha;
       
        hold on
               
        data1=ecoregions_pixelAR1(i,:);
        datain=[(1:length(data1))' data1'];
        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
        disp('Americas taub sig')
        disp([taub sig])
        taub_allecoregions=[taub_allecoregions taub];

end

box on


dim = [0.03+0.112-0.05 0.32+0.55 0.04 0.075];
str = strcat(sprintf('%.2f',taub_allecoregions(1)),'*');  
annotation('textbox',dim,'String',str,'BackgroundColor',line_colors(1,:),'edgecolor','none','FaceAlpha',0.7,'fontsize',9,'HorizontalAlignment','center','Margin',0.5);

dim = [0.03+0.165-0.05 0.32+0.55 0.04 0.075];
str = strcat(sprintf('%.2f',taub_allecoregions(3)),'*');  
annotation('textbox',dim,'String',str,'BackgroundColor',line_colors(2,:),'edgecolor','none','FaceAlpha',0.7,'fontsize',9,'HorizontalAlignment','center','Margin',0.5);

dim = [0.03+0.218-0.05 0.32+0.55 0.04 0.075];
str = strcat(sprintf('%.2f',taub_allecoregions(4)),'*');  
annotation('textbox',dim,'String',str,'BackgroundColor',line_colors(3,:),'edgecolor','none','FaceAlpha',0.7,'fontsize',9,'HorizontalAlignment','center','Margin',0.5);

dim = [0.03+0.271-0.05 0.32+0.55 0.04 0.075];
str = strcat(sprintf('%.2f',taub_allecoregions(2)),'*');  
annotation('textbox',dim,'String',str,'BackgroundColor',line_colors(4,:),'edgecolor','none','FaceAlpha',0.7,'fontsize',9,'HorizontalAlignment','center','Margin',0.5);
    

set(gca,'xtick',[1995+1/12:5:2015+1/12])

set(gca,'xticklabel',{'1995', '2000', '2005', '2010', '2015'})


set(gca,'ylim',[0 0.93])
set(gca,'ytick',0.3:0.2:0.7)


ylabel_pos=ylabel('AR(1)','FontSize',10); %,'FontWeight','bold   ,'color',red_dark_final
% set(ylabel_pos, 'Units', 'Normalized', 'Position', [-0.116, 0.44, 0]);

box on


set(gca,'xlim',[1991 2020])

set(gca,'fontsize',10)

ylim1=ylim;
xlim1=xlim;
hold on



% grid on
h=gca;
h.XRuler.TickLength=[0.01 0.01];

%% %%%%%% Africa  

pos2 = [0.32+0.08 0.69+0.0865 0.23 0.24];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

load AR1_Africa_lineardetrend_pixellevel


load Africa_tropical_belt_LUCC
load Africa_tropical_belt_fishnet


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


Final_time_series_FullRadar=all_scores(ind_1,:); 


ax = gca;
ax.YColor = 'k';

load Africa_tropical_belt_Ecoregion
Olson_ID=cell_dominantEcoregion;
ecoregion_id_africa{1,1}=30124;
ecoregion_id_africa{2,1}=[30104 30129];
ecoregion_id_africa{3,1}=30126;

Olson_ID=Olson_ID(ind_1);

ecoregion1_pixelAR1_ind=ismember(Olson_ID,ecoregion_id_africa{1,1});
ecoregion2_pixelAR1_ind=ismember(Olson_ID,ecoregion_id_africa{2,1});
ecoregion3_pixelAR1_ind=ismember(Olson_ID,ecoregion_id_africa{3,1});

ecoregion1_pixelAR1=Final_time_series_FullRadar(ecoregion1_pixelAR1_ind,:);
ecoregion2_pixelAR1=Final_time_series_FullRadar(ecoregion2_pixelAR1_ind,:);
ecoregion3_pixelAR1=Final_time_series_FullRadar(ecoregion3_pixelAR1_ind,:);


ecoregions_pixelAR1=[nanmean(ecoregion1_pixelAR1,1);nanmean(ecoregion2_pixelAR1,1);nanmean(ecoregion3_pixelAR1,1)];

greens = cbrewer('seq','Greens',10);
colors=greens([3 6 9],:); 
line_colors=colors;


for i=1:3
      
        s1=scatter(ar_time,ecoregions_pixelAR1(i,:),scatter_size,'MarkerFaceColor',line_colors(i,:),'MarkerEdgeColor','none'); %%%% ,symbols{i} ,'filled'  'none'
        s1.MarkerFaceAlpha = scatter_alpha;
        s1.MarkerEdgeAlpha = scatter_alpha;
       
        hold on
        
        data1=ecoregions_pixelAR1(i,:);
        datain=[(1:length(data1))' data1'];
        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
        disp('Africa taub sig')
        disp([taub sig])
        taub_allecoregions=[taub_allecoregions taub];

end

box on
  

dim = [0.414+0.112-0.07 0.32+0.55 0.04 0.075];
str = strcat(sprintf('%.2f',taub_allecoregions(5)),'*');  
annotation('textbox',dim,'String',str,'BackgroundColor',line_colors(1,:),'edgecolor','none','FaceAlpha',0.7,'fontsize',9,'HorizontalAlignment','center','Margin',0.5);

dim = [0.414+0.165-0.07 0.32+0.55 0.04 0.075];
str = strcat(sprintf('%.2f',taub_allecoregions(6)),'*');  
annotation('textbox',dim,'String',str,'BackgroundColor',line_colors(2,:),'edgecolor','none','FaceAlpha',0.7,'fontsize',9,'HorizontalAlignment','center','Margin',0.5);

dim = [0.414+0.218-0.07 0.32+0.55 0.04 0.075];
str = strcat(sprintf('%.2f',taub_allecoregions(7)),'*');  
annotation('textbox',dim,'String',str,'BackgroundColor',line_colors(3,:),'edgecolor','none','FaceAlpha',0.7,'fontsize',9,'HorizontalAlignment','center','Margin',0.5);


set(gca,'ylim',[0 0.93])
set(gca,'ytick',0.3:0.2:0.7)
set(gca,'yticklabel',{}) 

set(gca,'xlim',[1991 2020])


set(gca,'xtick',[1995+1/12:5:2015+1/12 ])
% set(gca,'xticklabel','')
set(gca,'xticklabel',{'1995', '2000', '2005', '2010', '2015'})


set(gca,'FontSize',10)
box on

ylim1=ylim;
xlim1=xlim;
hold on


% grid on
h=gca;
h.XRuler.TickLength=[0.01 0.01];

%% %%%%%% Asia  

pos2 = [0.59+0.08 0.69+0.0865 0.23 0.24];
axes('position',pos2,'XColor', [1,1,1], 'YColor', [1,1,1]) 

load AR1_Asia_lineardetrend_pixellevel


load Asia_tropical_belt_LUCC
load Asia_tropical_belt_fishnet

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



Final_time_series_FullRadar=all_scores(ind_1,:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax = gca;
ax.YColor = 'k';


load Asia_tropical_belt_Ecoregion
Olson_ID=cell_dominantEcoregion;
ecoregion_id_asia{1,1}=[40103 40102]; 
ecoregion_id_asia{2,1}=[10115 10105 10120 10121 10122 10128 10127]; 

Olson_ID=Olson_ID(ind_1);

ecoregion1_pixelAR1_ind=ismember(Olson_ID,ecoregion_id_asia{1,1});
ecoregion2_pixelAR1_ind=ismember(Olson_ID,ecoregion_id_asia{2,1});


ecoregion1_pixelAR1=Final_time_series_FullRadar(ecoregion1_pixelAR1_ind,:);
ecoregion2_pixelAR1=Final_time_series_FullRadar(ecoregion2_pixelAR1_ind,:);



ecoregions_pixelAR1=[nanmean(ecoregion1_pixelAR1,1);nanmean(ecoregion2_pixelAR1,1)];

blues = cbrewer('seq','Blues',9);
colors=blues([4 7],:); 
line_colors=colors;


for i=1:2 
      
        s1=scatter(ar_time,ecoregions_pixelAR1(i,:),scatter_size,'MarkerFaceColor',line_colors(i,:),'MarkerEdgeColor','none'); %%%% ,symbols{i} ,'filled'
        s1.MarkerFaceAlpha = scatter_alpha;
        s1.MarkerEdgeAlpha = scatter_alpha;
       
        hold on
        
        
        data1=ecoregions_pixelAR1(i,:);
        datain=[(1:length(data1))' data1'];
        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.05,0);
        disp('Asia taub sig')
        disp([taub sig])
        taub_allecoregions=[taub_allecoregions taub];

end

box on


dim = [0.764+0.112-0.07 0.32+0.55 0.04 0.075];
str = strcat(sprintf('%.2f',taub_allecoregions(8)));  
annotation('textbox',dim,'String',str,'BackgroundColor',line_colors(1,:),'edgecolor','none','FaceAlpha',0.7,'fontsize',9,'HorizontalAlignment','center','Margin',0.5);

dim = [0.764+0.165-0.07 0.32+0.55 0.04 0.075];
str = strcat(sprintf('%.2f',taub_allecoregions(9)),'*');  
annotation('textbox',dim,'String',str,'BackgroundColor',line_colors(2,:),'edgecolor','none','FaceAlpha',0.7,'fontsize',9,'HorizontalAlignment','center','Margin',0.5);


set(gca,'xtick',[1995+1/12:5:2015+1/12])
% set(gca,'xticklabel','')
set(gca,'xticklabel',{'1995', '2000', '2005', '2010', '2015'})


set(gca,'FontSize',10)

box on

set(gca,'xlim',[1991 2020])

set(gca,'ylim',[0 0.93])
set(gca,'ytick',0.3:0.2:0.7)
set(gca,'yticklabel',{}) 



ylim1=ylim;
xlim1=xlim;
hold on

% grid on
h=gca;
h.XRuler.TickLength=[0.01 0.01];

tightfig
%%
% print(gcf,'-dtiff','-r300',strcat('./Fig5.tif')) 
