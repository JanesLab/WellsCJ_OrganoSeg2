%load files from directory
files = dir('Gold Standard Masks');

for k = length(files):-1:1
    % remove non-folders
    
    if files(k).isdir 
        files(k) = [ ];
        continue
    end

    if ~contains(files(k).name, '.tif')
        files(k) = [ ];
        continue
    end
    
    % remove folders starting with .
    fname = files(k).name;
    if fname(1) == '.'
        files(k) = [ ];
    end
end

cd 'Gold Standard Masks'\


fullMask = [];
k = 0;
fullMasks = cell(1,k);
goldCCs = cell(1,k);
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,'_0.tif'))
        
        if ~isempty(fullMask) %store most recent full mask/ cc
            
            fullMasks{k} = fullMask;
            goldCCs{k} = cc;
        end
        k = k + 1;
        fullMask = imread(files(i).name);
        cc = bwconncomp(fullMask);
    else
        tempMask = imread(files(i).name);
        tempcc = bwconncomp(tempMask); 
        cc.PixelIdxList(end+1) = tempcc.PixelIdxList;
        cc.NumObjects = cc.NumObjects + 1;
        fullMask = tempMask + fullMask;
    end

end
fullMasks{k} = fullMask;
goldCCs{k} = cc;
cd ..\