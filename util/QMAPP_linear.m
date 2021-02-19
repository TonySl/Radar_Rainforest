function [COEFFs,RMSE,corre2,p2]=QMAPP_linear(OBSmont,SATmont)


    [b,~,~,~,~] = regress(OBSmont,[SATmont,ones(length(SATmont),1)]); % R2, the F, p value, error variance.
    COEFFs=b';    

    
    RMSE = sqrt(mean((OBSmont-SATmont).^2));
    [corres,p] = corrcoef(OBSmont,SATmont);
    corre2=corres(1,2);
    p2=p(1,2);
    
    % Interpolation of input data (SAT) to the same percentiles of
    % benchmark data (OBS)
%     SATint= interp1(PSAT,sort(SATmont),sort(POBS),'linear','extrap');
%     
%     % Computation of the differences between the CDFs of OBS and SAT data
%     DIFF=sort(OBSmont)-SATint;
    
    % Fitting of a polynomial curve to DIFF
%     COEFF= polyfit(SATint,DIFF, 5);
    
    % Evaluation of the polynomial curve to SAT data 
%     SATCDF= polyval(COEFF,SATmont)+SATmont;
%     SATCDF_ok(ID_SATmont) = SATCDF;
%     
    % Comparison of cdf curves estimated for: benchmark data (OBSmont), 
    % data to modify (SATmont), corrected data (SATCDF)
%     set(gcf,'position',[ 530, 190, 1111, 794])
%     subplot(3,4,i)
%     plot( sort(OBSmont),(1:length(OBSmont))/(length(OBSmont)+1),'Color',0.7*[1,1,1], 'linewidth',7)
%     hold on
%     plot(sort(SATmont),(1:length(SATmont))/(length(SATmont)+1), 'b-','linewidth',4)
%     plot( sort(SATCDF),(1:length(SATCDF))/(length(SATCDF)+1), 'r--', 'linewidth',2)
%     xlabel('data'), ylabel('Cumulative Density Function')
%     monthTitle =['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
%     title(monthTitle(i,:),'fontweight','bold','fontsize',10), grid on
%     if i==12, legend ('Reference data','Original biased data','Corrected data','Location','southeast'), end
%     
%     
%     M_STAT_OBS(i,1)= nanmean(OBSmont); V_STAT_OBS(i,1) = nanvar(OBSmont);
%     M_STAT_ST(i,1)= nanmean(SATCDF);   V_STAT_ST(i,1) = nanvar(SATCDF);
    
end

