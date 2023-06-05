function [out] = shift_W(f,rank,RANK)
E = 2;  S = 3; 
W = 4;  N = 5;

	out = f;

    [dest,src]=cart2D_shift(RANK,rank,[0,-1],1);
    
    if dest >= 1
		data = f(:,1);
		labSend(data, dest, W);
	end
	out(:,1:end-1) = f(:,2:end);
	if src >= 1
		data = labReceive(src, W); 
		out(:,end) = data;
	end
end

