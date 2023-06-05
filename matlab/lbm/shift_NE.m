function [out] = shift_NE(f,rank,RANK)
E = 2;  S = 3;
W = 4;  N = 5;

out = f;

out(1:end-1,2:end) = f(2:end,1:end-1);

%East Part
[dest,src]=cart2D_shift(RANK,rank,[0,1],1);

if dest >= 1
    data = f(2:end,end);
    labSend(data, dest, E);
end

if src >= 1
    data = labReceive(src, E);
    out(1:end-1,1) = data;
end

%North Part
[dest,src]=cart2D_shift(RANK,rank,[-1,0],1);

if dest >= 1
    data = f(1,1:end-1);
    labSend(data, dest, N);
end

if src >= 1
    data = labReceive(src, N);
    out(end,2:end) = data;
end


%North-East Part
[dest,src]=cart2D_shift(RANK,rank,[-1,1],1);

if dest >= 1
    data = f(1,end);
    labSend(data, dest, N+E);
end

if src >= 1
    data = labReceive(src, N+E);
    out(end,1) = data;
end

end

