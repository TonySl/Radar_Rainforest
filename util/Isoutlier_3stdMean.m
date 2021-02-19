 function [av]=Isoutlier_3stdMean(FF)
 

%   mk=median(FF);
%   M_d=mad(FF,0);
%   c=-1/(sqrt(2)*erfcinv(3/2));
%   smad=c*M_d;
%   tsmad=3*smad; 
%   tsmad=3*M_d;
%   av=(abs(FF-mk)>=tsmad);
  
  mk=nanmean(FF);
  
  M_d=nanstd(FF);
  tsmad=3*M_d;
  av=(abs(FF-mk)>=tsmad);
  
 end
