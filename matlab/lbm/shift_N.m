function [out] = shift_N(f,rank,RANK)
E = 2;  S = 3; 
W = 4;  N = 5;

	out = f;
    
    [dest,src]=cart2D_shift(RANK,rank,[-1,0],1);
    
	if dest >= 1
		data = f(1,:);
		labSend(data, dest, N);
	 end
	 out(1:end-1,:) = f(2:end,:);
	 if src >= 1
		data = labReceive(src, N); 
		out(end,:) = data;
	 end
end

