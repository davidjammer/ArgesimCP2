function [out] = shift_E(f,rank,RANK)
E = 2;  S = 3; 
W = 4;  N = 5;

	out = f;
    
    [dest,src]=cart2D_shift(RANK,rank,[0,1],1);
    
    if dest >= 1
	    data = f(:,end);
	    labSend(data, dest, E);
    end
    out(:,2:end) = f(:,1:end-1);
    if src >= 1
	   data = labReceive(src, E); 
	   out(:,1) = data;
    end
end

