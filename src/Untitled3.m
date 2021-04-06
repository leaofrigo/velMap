close all;clear all; clc
[c1,h1]=contourf(peaks,8);
hleve=get(h1,'LevelList');
hchil = get(h1,'Children');
cont_level = zeros(1,length(hchil));
legend_entries = cell(1,length(hleve)-1);
TFtot = true(1,length(hchil));
for dd= 1:length(hchil)
    cont_level(dd) = get(hchil(dd),'userdata');
    TF = cont_level(dd)==cont_level;
    if sum(TF(:)) == 1
        TFtot(dd) = true;
    else
        TFtot(dd) = false;
    end
    if hleve(1)==cont_level(dd)
        TFtot(dd) = false;
    end
end
for kg = 1:length(hleve)-1
    legend_entries{kg} = [num2str(hleve(kg),3) '<x<' num2str(hleve(kg+1),3)];
end
[B,I]=sort(cont_level(TFtot));
B1 = hchil(TFtot);
legend(B1(I),legend_entries)