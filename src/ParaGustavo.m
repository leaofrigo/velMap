clear all;close all;clc
directory = 'C:\Users\Roberta\Desktop\Ricardo\Gustavo\';
listing = dir(directory);
a = cell(length(listing),1);
tf = zeros(length(listing),1);
k = a;
for i = 1:length(listing)
    a{i} = listing(i,1).name;
    k{i} = strfind(a{i}, '.sum');
    tf(i) = ~isempty(k{i});

end
OpenFiles = a(logical(tf));

delete([directory 'AllFiles.xls'])
for n = 1:length(OpenFiles)
    A(n,1).name = ['File_' OpenFiles{n}(1:end-4)];
    A(n,1).data = importdata([directory OpenFiles{n}]);
    xlswrite([directory 'AllFiles.xls'],[{'File Name:' OpenFiles{n}(1:end-4)}]...
        ,'main',[char((n*5)-4+'A'-1) '1']);
    xlswrite([directory 'AllFiles.xls'],[{'Depth (m)' 'Latitude (deg)' 'Longitude (deg)'}]...
        ,'main',[char((n*5)-3+'A'-1) '2']);
    xlswrite([directory 'AllFiles.xls'],[A(n,1).data.data(:,1:3)],'main',[char((n*5)-3+'A'-1) '3']);
    xlswrite([directory 'AllFiles.xls'],[A(n,1).data.textdata(:,2)]...
        ,'main',[char((n*5)-4+'A'-1) '2']);
        
end


    



