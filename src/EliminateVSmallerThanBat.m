function [Vel] = EliminateVSmallerThanBat(X1,Y1,Z1,Vel,XrR,YrR,FrR)
Space = size(X1);
for j=1:Space(1)
    a=Y1(j,1,1);
    tmp = abs(XrR(1,:)-a);
    [~, idx]=min(tmp);
    TFXrR=XrR==XrR(1,idx);
    y=YrR(TFXrR);
    for i=1:Space(2)
       b=X1(1,i,1);
       tmpy = abs(y-b);
       [~, idy]=min(tmpy);
       TFYrR=YrR==YrR(idy,1);
       z = FrR(and(TFXrR,TFYrR));
        for k=1:Space(3)
           c = Z1(1,1,k);
           if c<z;
               Vel(j,i,k)=nan;
           end           
        end
    end
end