function OpenFiles = DirectoryOpenFiles(Directory)
listing = dir(Directory);
a = cell(length(listing),1);
tf = zeros(length(listing),1);
for i = 1:length(listing)
    a{i} = listing(i,1).name;
    k{i} = strfind(a{i}, '.mat');
    tf(i) = ~isempty(k{i});

end
OpenFiles = a(logical(tf));
