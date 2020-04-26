function plot_lvlh_error(tvec,x_error, SC, font_size)

figure
axes_labels = {'x', 'y', 'z'};

for ax = 1:3
subplot(3,1,ax)

plot(tvec/MINUTES, x_error(:,ax), 'k')
ylabel({axes_labels{ax} ' [m]'}, 'Interpreter', 'Latex','FontSize', font_size);

xlabel('$Time$ [min]', 'Interpreter', 'Latex','FontSize', font_size);
grid on

set(gca,'fontsize',font_size-2)
end

end