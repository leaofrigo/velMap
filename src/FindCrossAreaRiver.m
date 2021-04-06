function [Area] = FindCrossAreaRiver(depth,xx)
xx = xx(2:end);
a = zeros(length(xx)-1,1);
for n = 1:length(xx)-1
    a(n) = (xx(n+1) - xx(n))*(depth(n)+depth(n+1))*.5;
end
Area = sum(a);