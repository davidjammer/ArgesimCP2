function plot_u(file)
	data = load(file);

	imagesc(data); % plot macroscopic velicity magnitude
	axis('equal'); % make display square
	colorbar; % show color index
	title(['Relative macroscopic velocity magnitude (u/u_0) after ',...
	num2str(2000),' iterations']); % show plot title
end