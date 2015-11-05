function [] = codeare_hdf5write (varargin) 
    fname = varargin{1};
    display(sprintf('filename = %s',fname));
    delf = sprintf('rm %s', fname);
    [~,~] = system(delf);
    for i=2:nargin
        vname = sprintf('/%s',inputname(i));
        if (iscomplex(varargin{i}))
            varargin{i} = reshape([real(varargin{i}(:))';imag(varargin{i}(:))'],[2 size(varargin{i})]);
        end
        h5create (fname,vname,size(varargin{i}));
        h5write (fname,vname,varargin{i});
    end
end