function plot_u(file,rep)
	data = load(file);

	imagesc(data); % plot macroscopic velicity magnitude
	axis('equal'); % make display square
	colorbar; % show color index
	title(['Relative macroscopic velocity magnitude (u/u_0) after ',...
	num2str(rep),' iterations']); % show plot title
end