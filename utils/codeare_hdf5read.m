function [] = codeare_hdf5read (fname)
%    
%   Read codeare hdf5 output.
%   
%   [data,names] = codeare_hdf5read (fname);
%
%   IN: 
%     fname: File name
%   OUT: 
%     data:  Cell array of matrices
%     names: The named under which the data was saved in HDF5 file
%   
%   Kaveh Vahedipour - Research Centre Juelich, 2013
    
    hinfo = hdf5info(fname);
    nsets = numel(hinfo.GroupHierarchy.Datasets);
    dsets = hinfo.GroupHierarchy.Datasets;
    data  = cell(nsets);
    names = cell(nsets);
    
    %display (sprintf('   reading %d datasets ...', nsets));
    
    for i = 1:nsets
        names{i} = dsets(i).Name(2:end);
        %display (sprintf('     %s', names{i}));
        tmp = hdf5read (dsets(i));
        tsz = size(tmp);
        tn  = numel(tmp);
        if (tsz(1) == 2)
            tmp = tmp(:);
            tmp = tmp(1:2:end)+1i*tmp(2:2:end);
            if (numel(tsz)==2)
                tsz = [tsz 1];
            end
            tmp = reshape(tmp,tsz(2:end));
        end
        assignin ('caller',names{i},tmp);
    end

    clear tmp names tsz tn data dsets nsets hinfo cmd i fname;
    whos
    %display ('   done.');
    
end