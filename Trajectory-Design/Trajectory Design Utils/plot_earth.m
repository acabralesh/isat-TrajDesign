function plot_earth(r_radius, xoffset)
%plot_earth(r_radius, xoffset) plots a sphere on current figure
%
% Purpose: Automatically plots earth for a given radius and xoffset. Note
% units MUST match accros inputs
% 
% Inputs:
% r_radius  [1x1] radius of earth in desired units
% xoffset   [3x1] offset of earth, if input a [1x1] vector offest is
%           assumed to be only on x direction
%
% Outputs:
% None
%
% Co-dependencies
%   sphere
%
% Source: None
% 
% Author: Alejandro Cabrales H
%  
[xearth, yearth, zearth] = sphere;

switch length(xoffset)
    
    case 1
        
        hSurface=surf(xoffset + xearth*r_radius,yearth*r_radius,zearth*r_radius);
        set(hSurface,'FaceColor',[0 0 1], ...
            'FaceAlpha',0.5,'FaceLighting','gouraud')
        
    case 2
        hSurface=surf(xoffset(1) + xearth*r_radius,xoffset(2) +yearth*r_radius,zearth*r_radius);
        set(hSurface,'FaceColor',[0 0 1], ...
            'FaceAlpha',0.5,'FaceLighting','gouraud')
        
    case 3
        hSurface=surf(xoffset(1) + xearth*r_radius,xoffset(2) + ...
            yearth*r_radius, xoffset(3) + zearth*r_radius);
        set(hSurface,'FaceColor',[0 0 1], ...
            'FaceAlpha',0.5,'FaceLighting','gouraud')     
end
end