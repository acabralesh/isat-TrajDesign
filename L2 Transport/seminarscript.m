% seminar slides
close all;

font_size = 22;
if 1
figure; 

% plot3(1 - mu, 0, 0, '*b')
hold on
% plot3(1+ geo_orb(:,1), geo_orb(:,2), geo_orb(:,3), '--k')
% plot3(1+ moon_orb(:,1), moon_orb(:,2), moon_orb(:,3), '--b')
plot3(xnum, ynum, znum,'k-', 'linewidth', 2)
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% view([-16, 5])
axis equal
view(0,90) 

grab_ind = 5:5:35;
grab_ind = 1:10;

mani_plots = orb_sort_ind(grab_ind);



for mani = 1:length(mani_plots)
    
    
xnumvec = sol{mani_plots(mani)}.xnumvec;

plot3(xnumvec(1,1),xnumvec(1,2),xnumvec(1,3),'sk','MarkerFaceColor', [0 0 0], 'linewidth', 0.1)

plot3(xnumvec(:,1),xnumvec(:,2),xnumvec(:,3),'--','MarkerFaceColor', [0 0 0], 'linewidth', 0.1)

end
end

%%


% make gif
grab_ind = 1:11;

tindices = 0:-0.05:sol{1}.t(end);

mani_plots = orb_sort_ind(grab_ind);
h = figure;
axis equal % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
set(gcf, 'color', 'w')

v = VideoWriter('testvideo');
open(v);

hold on
% plot3(1+ geo_orb(:,1), geo_orb(:,2), geo_orb(:,3), '--k')
% plot3(1+ moon_orb(:,1), moon_orb(:,2), moon_orb(:,3), '--b')
plot3(xnum, ynum, znum,'k-', 'linewidth', 2)
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% view([-16, 5])
axis equal
view(0,90) 


for mani = 1:length(mani_plots)
% ii = 17;

xnumvec = sol{mani_plots(mani)}.xnumvec;
t = sol{mani_plots(mani)}.t;
plot3(xnumvec(1,1),xnumvec(1,2),xnumvec(1,3),'sk','MarkerFaceColor', [0 0 0], 'linewidth', 0.1)



start_plot = 0;
keep_plot = 0;
for td = 2:(length(tindices))
    
    ind2plot = length(t(t> tindices(td)));
    
    % Draw plot for y = x.^n
    plot3(xnumvec(1:ind2plot,1),xnumvec(1:ind2plot,2),xnumvec(1:ind2plot,3), 'k', 'linewidth', 0.1)
    postsat = plot3(xnumvec(ind2plot,1),xnumvec(ind2plot,2),xnumvec(ind2plot,3), 'ok', 'linewidth', 0.1);
    % check if need to plot earth
    
    
    
    if ~isempty(find(xnumvec(1:ind2plot,1) <= 1,1))
        start_plot = 1;
    end
    
    
    if start_plot == 1 || keep_plot == 1
        
        plot3(1-mu, 0, 0, 'ob', 'MarkerFaceColor', [0 0 1])
        keep_plot = 1;
    end
    
    drawnow
    
      % Capture the plot as an image 
      frame = getframe(h); 
      delete(postsat);
      writeVideo(v,frame);
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if td == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
end
end
close(v);


%% Plot closest trajectories
if 0 
grab_ind = 1:7;

mani_plots = orb_sort_ind(grab_ind);

figure;
hold on
plot3(1-mu, 0, 0, 'ob', 'MarkerFaceColor', [0 0 1])
plot3(xnum, ynum, znum,'k-', 'linewidth', 1)

for mani = 1:length(mani_plots)
    
    xnumvec = sol{mani_plots(mani)}.xnumvec;
    
    plot3(xnumvec(:,1),xnumvec(:,2),xnumvec(:,3), 'linewidth', 0.5)
    plot3(xnumvec(minind(mani_plots(mani)),1),xnumvec(minind(mani_plots(mani)),2),...
        xnumvec(minind(mani_plots(mani)),3), 'ks')
    
end
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% view([-16, 5])
axis equal
xlim([.99, 1.012]);
ylim([-7.2e-3 7.2e-3]);
zlim([-1e-3 1e-3]);


%% Plot the trajectories


num_radia = 3;

figure;
hold on
% --- earth ------
% plot3(1 - mu, 0, 0, '.b')
r_radius = 6371*KILOMETERS/a;
x_offset = (1-mu)*0;
plot_earth(1, x_offset)

% ----earth

plot3((xnumvec_orb(:,1)- (1-mu))/r_radius, xnumvec_orb(:,2)/r_radius, xnumvec_orb(:,3)/r_radius, 'k')
%  -----manifold------
ind_sol = orb_sort_ind(mani_ind);
xnumvec_mani = sol{ind_sol}.xnumvec;
plot3((xnumvec_mani(:,1)-(1-mu))/r_radius,xnumvec_mani(:,2)/r_radius,xnumvec_mani(:,3)/r_radius, '--b')
% plot the closest approach
plot3((xnumvec_mani(minind(ind_sol),1)-(1-mu))/r_radius,xnumvec_mani(minind(ind_sol),2)/r_radius,...
        xnumvec_mani(minind(ind_sol),3)/r_radius, 'k.')
    
% -----manifold end ----
axis equal
xlabel('$x$ [$R_E$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$R_E$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$R_E$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on

xlim([-num_radia, num_radia]);
ylim([-num_radia, num_radia]);
zlim([-num_radia, num_radia]);
leg = legend({'Earth', 'Parking Orb', 'Mainfold', 'Int. Point'}, 'Interpreter', 'Latex');
leg.FontSize = font_size;

end