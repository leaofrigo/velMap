function [Perimeter] = FindPerimeterRiver(depthnew1,tracknew1)
p = tracknew1(end)-tracknew1(1);
p1 = sum(abs(diff(depthnew1)));
Perimeter = 2*p+p1;