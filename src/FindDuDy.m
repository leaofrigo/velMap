function MAX = FindDuDy(uprime1,yloc)
MAX = 0;
for i= 1:length(uprime1)
    for j = i+1:length(uprime1)
        VALUE = abs(uprime1(i)-uprime1(j))/abs(yloc(i)-yloc(j));
        if VALUE > MAX
            MAX = VALUE;
        end
    end
end