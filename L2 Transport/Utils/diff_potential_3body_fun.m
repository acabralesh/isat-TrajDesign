function [dU_vec] = diff_potential_3body_fun(mu,x,ax)
%DIFF_POTENTIAL_3BODY_FUN Computes the pseudo potential of a function
%
% Inputs:
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2
% x         [3x1] Non dimensional position vector
% ax        [1xn] desired axis for example 1:3 will give derivative for
%           each axis

len_ax = length(ax);
dU_vec = zeros(len_ax,1);

% Relative position vector from primary to desired vector
d = sqrt(sum((x(1)+mu).^2 + x(2).^2 + x(3).^2));

% Relative position vector from secondary to desired vector
r = sqrt( sum((x(1) - 1 + mu).^2  + x(2).^2 + x(3).^2 ));

for ii = 1:len_ax
    
    dU_temp = NaN;
    
    switch ax(ii)
        
        case 1  
            dU_temp = x(1) - (1 - mu).*(x(1) + mu)./d.^3 - ...
                mu*(x(1) - (1 - mu))./r.^3;
            
        case 2
            dU_temp = x(2) - (1 - mu).*x(2)./d.^3 - mu*x(2)./r.^3;
            
        case 3
            dU_temp = -(1 - mu)*x(3)./d.^3 - mu*x(3)./r^.3;
            
    end
    
    dU_vec(ii) = dU_temp;
end
          
end

