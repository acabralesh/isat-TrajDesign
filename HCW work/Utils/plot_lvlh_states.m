function plot_lvlh_states(x_lvlh,font_size)
%PLOT_LVLH_STATES Plots the position in lvlh frame

figure;
subplot(2,2,1)
plot(x_lvlh(:,1)/METERS, x_lvlh(:,2)/METERS, 'k-', 'linewidth', 1)
xlabel('$x$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
ylabel('$y$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
set(gca,'fontsize',font_size-2)
grid on
axis square
subplot(2,2,2)
plot(x_lvlh(:,2)/METERS, x_lvlh(:,3)/METERS, 'k-', 'linewidth', 1)
xlabel('$y$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
ylabel('$z$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
set(gca,'fontsize',font_size-2)
grid on
axis square
subplot(2,2,3)
plot(x_lvlh(:,1)/METERS, x_lvlh(:,3)/METERS, 'k-', 'linewidth', 1)
xlabel('$x$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
ylabel('$z$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
set(gca,'fontsize',font_size-2)
grid on
axis square
subplot(2,2,4)
plot3(x_lvlh(:,1)/METERS,x_lvlh(:,2)/METERS, x_lvlh(:,3)/METERS, 'k-', 'linewidth', 1)
xlabel('$x$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
ylabel('$y$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
zlabel('$z$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
set(gca,'fontsize',font_size-2)
grid on
axis square


end

