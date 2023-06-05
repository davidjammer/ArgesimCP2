function [out] = shift_NW(f,rank,RANK)
E = 2;  S = 3;
W = 4;  N = 5;

out = f;

out(1:end-1,1:end-1) = f(2:end,2:end);

%West Part
[dest,src]=cart2D_shift(RANK,rank,[0,-1],1);

if dest >= 1
    data = f(2:end,1);
    labSend(data, dest, W);
end

if src >= 1
    data = labReceive(src, W);
    out(1:end-1,end) = data;
end

%North Part
[dest,src]=cart2D_shift(RANK,rank,[-1,0],1);

if dest >= 1
    data = f(1,2:end);
    labSend(data, dest, N);
end

if src >= 1
    data = labReceive(src, N);
    out(end,1:end-1) = data;
end


%North-West Part
[dest,src]=cart2D_shift(RANK,rank,[-1,-1],1);

if dest >= 1
    data = f(1,1);
    labSend(data, dest, N+W);
end

if src >= 1
    data = labReceive(src, N+W);
    out(end,end) = data;
end

end

