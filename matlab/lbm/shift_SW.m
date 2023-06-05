function [out] = shift_SW(f,rank,RANK)
E = 2;  S = 3;
W = 4;  N = 5;

out = f;

out(2:end,1:end-1) = f(1:end-1,2:end);

%West Part
[dest,src]=cart2D_shift(RANK,rank,[0,-1],1);

if dest >= 1
    data = f(1:end-1,1);
    labSend(data, dest, W);
end

if src >= 1
    data = labReceive(src, W);
    out(2:end,end) = data;
end

%South Part
[dest,src]=cart2D_shift(RANK,rank,[1,0],1);

if dest >= 1
    data = f(end,2:end);
    labSend(data, dest, S);
end

if src >= 1
    data = labReceive(src, S);
    out(1,1:end-1) = data;
end


%South-West Part
[dest,src]=cart2D_shift(RANK,rank,[1,-1],1);

if dest >= 1
    data = f(end,1);
    labSend(data, dest, S+W);
end

if src >= 1
    data = labReceive(src, S+W);
    out(1,end) = data;
end

end

