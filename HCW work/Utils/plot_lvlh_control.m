function plot_lvlh_control(tvec,xvec, SC, font_size)
% Plot the controls

figure

for ax = 1:3
subplot(3,1,ax)
if SC.CTRL.mode == SC.CTRL.impulsive
    
    stem(tvec/MINUTES,xvec.delta_v(:,ax), 'k') 
   ylabel('$\Delta v$ [m/s]', 'Interpreter', 'Latex','FontSize', font_size);
elseif  SC.CTRL.mode == SC.CTRL.zoh
    
    stairs(tvec/MINUTES,xvec.force_lvlh(:,ax), 'k')
    ylabel('$Force$ [N]', 'Interpreter', 'Latex','FontSize', font_size);
end

xlabel('$Time$ [min]', 'Interpreter', 'Latex','FontSize', font_size);
grid on
set(gca,'fontsize',font_size-2)

end


end