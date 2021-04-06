VerticeYY=zeros((Sumnum(2)-1)*4,1);
VerticeXX=zeros((Sumnum(2)-1)*4,1);
w=1;
ww=1;
Face = zeros((Sumnum(2)-1),4);
YY=zeros(Sumnum(end),1);
cx = 0;
for n = 1:length(Sumnum2)-1
    k = Sumnum2(n);
    j = Sumnum2(n+1);
    kk = Sumnum(n);
    jj = Sumnum(n+1);
    for i = k:j-1
        VerticeYY(k+2*(i-1)+1-k:k+2*(i-1)-k+2) = ones(2,1)*cellDepths(i);
    end
    for ii = kk:jj-1        
        YY(ii) = mean(cellDepths(ii+cx:ii+cx+1));
    end
    cx=cx+1;
    VerticeXX(w:2:w+(j-k-1)*2) = ones(i-k+1,1)*track(n);
    VerticeXX(w+1:2:w+(j-k-1)*2+1) = ones(i-k+1,1)*track(n+1);
    XX(ww:ww+(j-k-2),1) = ones(i-k,1)*mean(track(n:n+1));
    w = (j-k)*2+ w;
    ww = (j-k)+ ww -1;
end
YY = YY(1:end-1);
VerticeXXYY = [VerticeXX, VerticeYY];
g=1;
r=1;
for nn = 1:length(Sumnum2)-1
% for nn = 1:2
    k = Sumnum2(nn);
    j = Sumnum2(nn+1);
    
    for ii = r:(j-k-2)+r
        
        Face(ii,1:4) = g:g+3;
        g = g+2;
        
    end
    g= g+2;
    r= (j-k-1)+r;
    
end
Face =[Face(:,1) Face(:,2) Face(:,4) Face(:,3)];
p = patch('Faces',Face,'Vertices',VerticeXXYY);
colorbar
set(gca,'CLim',[-1 1])
cdata = Quiver3VectorU;
set(p,'FaceColor','flat',...
'FaceVertexCData',cdata,...
'CDataMapping','scaled','EdgeColor','none')
% Face = [1 2 3 4; ...
%         3 4 5 6; ...
%         5 6 7 8; ...
%         9 10 11 12; ...
%         11 12 13 14; ...
%         13 14 15 16; ...
%         15 16 17 18; ...
%         19 20 21 22; ...
%         21 22 23 24];