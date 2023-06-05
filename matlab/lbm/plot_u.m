function plot_u(u,it)
% plot macroscopic velocity magnitude
  imagesc(u);   
  axis("off")
  axis("image")
  colormap("turbo")
  colorbar("FontSize",12); 

  if nargin == 2
    nx = size(u,1);
    sTitle = sprintf( ...
      "Macroscopic velocity |u|/u_{0} (nx=%d) after %d iterations", nx, it);
    title(sTitle)
  end
end
