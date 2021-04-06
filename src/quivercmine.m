function quivercmine(x,y,u,v,c)
if length(x)~=length(y)
    error('All vector lengths must be the same');
elseif length(x)~=length(u)
    error('All vector lengths must be the same');
elseif length(x)~=length(v)
    error('All vector lengths must be the same');
elseif length(x)~=length(c)
    error('All vector lengths must be the same');
end

for i=1:length(x)
    quiver(x(i),y(i),u(i),v(i),'Color',c(i))
end