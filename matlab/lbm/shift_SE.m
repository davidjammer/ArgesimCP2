function [out] = shift_SE(f,rank,RANK)
E = 2;  S = 3;
W = 4;  N = 5;

out = f;

out(2:end,2:end) = f(1:end-1,1:end-1);

%East Part
[dest,src]=cart2D_shift(RANK,rank,[0,1],1);

if dest >= 1
    data = f(1:end-1,end);
    labSend(data, dest, E);
end

if src >= 1
    data = labReceive(src, E);
    out(2:end,1) = data;
end

%South Part
[dest,src]=cart2D_shift(RANK,rank,[1,0],1);

if dest >= 1
    data = f(end,1:end-1);
    labSend(data, dest, S);
end

if src >= 1
    data = labReceive(src, S);
    out(1,2:end) = data;
end


%South-East Part
[dest,src]=cart2D_shift(RANK,rank,[1,1],1);

if dest >= 1
    data = f(end,end);
    labSend(data, dest, S+E);
end

if src >= 1
    data = labReceive(src, S+E);
    out(1,1) = data;
end

end

