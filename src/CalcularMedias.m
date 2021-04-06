function [AVERAGES] = CalcularMedias(xx,zz,AverageVelocityVector,AverageVelocityMag,VelExtap)
zz = zz(2:end);
xxdiff = diff(xx);
Ytop = VelExtap.Ytop;
diffYtop = diff(Ytop);
Ybot = VelExtap.Ybot;
diffYbot = diff(Ybot);
TopVel = VelExtap.VelocityTopExtrp;
BotVel = VelExtap.VelocityBotExtrp;
s = size(diffYbot);
xxdiff1=repmat(xxdiff,s(1),1);
for i=1:s
    TopVelMean(i,:) = mean(TopVel(i:i+1,:));
    BotVelMean(i,:) = mean(BotVel(i:i+1,:));
end
  
AverageQBody = (nansum([xxdiff'.*zz'.*AverageVelocityVector(:,1) xxdiff'.*zz'.*AverageVelocityVector(:,2)...
    xxdiff'.*zz'.*AverageVelocityVector(:,3)]));
TotalAverageVelocityVector = AverageQBody/nansum(sum(zz'.*xxdiff'));
TotalAverageVelocity2 = nansum(xxdiff.*zz.*AverageVelocityMag)/nansum(nansum(zz'.*xxdiff'));
TotalAverageVelocity = sqrt(nansum(TotalAverageVelocityVector.^2));
AverageBothAverages = mean([TotalAverageVelocity,TotalAverageVelocity2]);
QMiddle = nansum(xxdiff'.*zz'.*sqrt(AverageVelocityVector(:,1).^2+...
    AverageVelocityVector(:,2).^2+AverageVelocityVector(:,3).^2));

AVERAGES.QTop = nansum(nansum(diffYtop.*xxdiff1.*TopVelMean));
AVERAGES.VelTop = nansum(nansum(diffYtop.*xxdiff1.*TopVelMean))/nansum(nansum(diffYtop.*xxdiff1));
AVERAGES.QBot = nansum(nansum(diffYbot.*xxdiff1.*BotVelMean));
AVERAGES.VelBot = nansum(nansum(diffYtop.*xxdiff1.*BotVelMean))/nansum(nansum(diffYbot.*xxdiff1));
AVERAGES.QMiddle = QMiddle;

AVERAGES.AverageQBody = AverageQBody;
AVERAGES.TotalAverageVelocityVector = TotalAverageVelocityVector;
AVERAGES.TotalAverageVelocity2 = TotalAverageVelocity2;
AVERAGES.TotalAverageVelocity = TotalAverageVelocity;
AVERAGES.AverageBothAverages = AverageBothAverages;
AVERAGES.AreaMeasured = nansum(sum(zz'.*xxdiff'));
end