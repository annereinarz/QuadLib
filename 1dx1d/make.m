function make(varargin)
if(size(varargin,2) < 1)
    target = 'all';
else
    target = varargin{1};
end

switch (target)
    case 'all'
        make 1d
    case '1d'
        make squad1d
    case 'squad1d'
        mex -o squad1d             squad1d_gateway.c        ...
            squad1d.c                ...
            sing_identical1d.c       ...
            sing_commonvertex1d.c    ...
            Ref2PhyS1d_affine.c      ...
            Ref2PhyS1d_circle.c      ...
            quadpoints.c
end

end

