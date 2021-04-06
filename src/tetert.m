clear Face
g=1;
r=1;
% for nn = 1:length(Sumnum2)-1
for nn = 1:3
    k = Sumnum2(nn);
    j = Sumnum2(nn+1);
    
    for ii = r:(j-k-2)+r
        
        Face(ii,1:4) = g:g+3;
        g = g+2;
        
    end
    g= g + 2;
    r= (j-k-1)+r;
    
end