function [LowerLim, UpperLim,s] = limites(Vel)

% MaxVe = max(max(Vel));
% MaxVe = round(MaxVe*10)/10;
% MinVe = min(min(Vel));
% MinVe = round(MinVe*10)/10;
% TFVE = 0 == isnan(Vel);
% NumberOfSamplesVel = sum(sum(TFVE));
% 
% while NumOfSampleMax <= NumberOfSamplesVel*.95 || NumOfSampleMin <= NumberOfSamplesVel*.95
%     if NumOfSampleMax <= NumberOfSamplesVel*.95
%         NumOfSampleMax = sum(sum(Vel <= MaxVe))
%         MaxVe = MaxVe - .10;
%         
%     end
%     if NumOfSampleMin <= NumberOfVel*.95
%         
%     end
%     
%     
% 
% end
[m,s] = nanmean3(Vel);
LowerLim = floor((m - 2*s)*10)/10;
UpperLim = ceil((m + 2*s)*10)/10;