% t :: nxd matrix (n zeilen, d spalten)
function [t,wt] = import_quadpoints(fname)
    f = fopen(fname, 'rb');
    n = fread(f, 1, 'int');
    d = fread(f, 1, 'int');
    t = fread(f, [n,d], 'double');
    wt = fread(f, n, 'double');
    fclose(f);
end
