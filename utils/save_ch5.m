function [] = save_ch5 (varargin) 
    fname = varargin{1};
    display(sprintf('filename = %s',fname));
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