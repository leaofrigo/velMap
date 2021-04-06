function [mean,std] = nanmean3(x)
x = reshape(x,1,[]);
TF = 0 == isnan(x);
mean = sum(sum(TF).*nanmean(x))/sum(sum(TF));
std = sqrt(sum(sum((x(TF)-mean).^2))/(sum(sum(TF))-1));