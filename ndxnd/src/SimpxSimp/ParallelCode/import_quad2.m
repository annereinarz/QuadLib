% t :: nxd matrix (n zeilen, d spalten)
function [u,v,zh,wt] = import_quad2(fname)
    f = fopen(fname, 'rb');

    luv = fread(f, 1, 'int');
    ruv = fread(f, 1, 'int');

    u = fread(f, [ruv,luv], 'double');
    v = fread(f, [ruv,luv], 'double');

    lzh = fread(f, 1, 'int');
    rzh = fread(f, 1, 'int');

    zh = fread(f, [rzh,lzh], 'double');

    wt = fread(f, ruv, 'double');

    fclose(f);
end
