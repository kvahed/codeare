function [] = save_ch5 (varargin) 
%SAVE_CH5 Writes data into codeare HDF5 files.
%    
%   usage: save_ch5 (filename, data1, data2, ...)
%   
%   Kaveh Vahedipour - NYU School of Medicine, 2015

    delf = sprintf('rm %s', fname);
    [~,~] = system(delf);
    for i=2:nargin
        vname = sprintf('/%s',inputname(i));
        cplx  = iscomplex(varargin{i});
        if (cplx)
            varargin{i} = reshape([real(varargin{i}(:))';imag(varargin{i}(:))'],[2 size(varargin{i})]);
        end
        h5create (fname,vname,size(varargin{i}),'Datatype',class(varargin{i}));
        h5write (fname,vname,varargin{i});
        if (cplx)
            h5writeatt(fname,vname,'complex',1);
        end
    end
end