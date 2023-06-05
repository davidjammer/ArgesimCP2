function [dest, source] = cart2D_shift(CART,rank,direction,disp)
%CART: Map of labindex
%rank: coordinate in CART
%direction: dimension of shift operation
%disp: displacement

[ny,nx] = size(CART);

dest_rank = rank;

dest_rank(1) = rank(1) + direction(1) * disp;
dest_rank(2) = rank(2) + direction(2) * disp;
if 1 <= dest_rank(1) && dest_rank(1) <= ny && 1 <= dest_rank(2) && dest_rank(2) <= nx
    dest = CART(dest_rank(1), dest_rank(2));
else
    dest = -1;
end

source_rank = rank;

source_rank(1) = rank(1) - direction(1) * disp;
source_rank(2) = rank(2) - direction(2) * disp;
if 1 <= source_rank(1) && source_rank(1) <= ny && 1 <= source_rank(2) && source_rank(2) <= nx
    source = CART(source_rank(1), source_rank(2));
else
    source = -1;
end

end

