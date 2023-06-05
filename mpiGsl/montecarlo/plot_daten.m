function nSteps = plot_daten(file)
	data = load(file);

  plot(data(:,1), data(:,2));

  grid on;
  xlabel('t');
  ylabel('y(t)');
  legend('Weg');
  
  nSteps = size(data, 1) - 1;
end
