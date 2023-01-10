function Jp = calcMaxLvlDWT(Y,M)

    N = size(Y,1);

    Jp = fix(log2(N/(M-1)));
end