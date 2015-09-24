% t :: nxd matrix (n zeilen, d spalten)
function export_quadpoints(t, wt, fname)
    [n,d] = size(t);
    f = fopen(fname, 'wb');
    fwrite(f, n, 'int');
    fwrite(f, d, 'int');
    fwrite(f, t, 'double');
    fwrite(f, wt, 'double');
    fclose(f);
end
