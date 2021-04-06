c= -TotalAverageVelocityVector(2)/TotalAverageVelocityVector(1);
d=-AveLong*c + AveLat;
t = AveLong-.00022:.000001:AveLong+.00032;
r= c*t + d;
plot(t,r)
hold on
scalefactor = .0005;
quiver(AveLong,AveLat,TotalAverageVelocityVector(2)*scalefactor,TotalAverageVelocityVector(1)*scalefactor)
xlim([Xlim(1)-.005 Xlim(2)+.005]);
ylim([Ylim(1)-.005 Ylim(2)+.005]);
theta = atan(c);
plot_google_map