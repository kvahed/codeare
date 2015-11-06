function [] = load_ch5 (fname)
%LOAD_CH5 Reads data from HDF5 files.
%    
%   usage: load_ch5 (filename)
%   
%   Kaveh Vahedipour - NYU School of Medicine, 2015
    
    hinfo = h5info(fname);
    nsets = numel(hinfo.Datasets);
    dsets = hinfo.Datasets;
    names = cell(nsets);
    
    for i = 1:nsets
        names{i} = dsets(i).Name;
        cplx = 0;
        if numel(dsets(i).Attributes)>0
            for j=1:numel(dsets(i).Attributes);
                if (strcmp(dsets(i).Attributes(j).Name,'complex'))
                    cplx = dsets(i).Attributes(j).Value;
                end
            end
        end 
        tmp = h5read (fname,['/' names{i}]);
        tsz = size(tmp);
        tn  = numel(tmp);
        if (tsz(1) == 2 && cplx == 1)
            tmp = tmp(:);
            tmp = tmp(1:2:end)+1i*tmp(2:2:end);
            if (numel(tsz)==2)
                tsz = [tsz 1];
            end
            tmp = reshape(tmp,tsz(2:end));
        end
        assignin ('caller',names{i},tmp);
    end
end