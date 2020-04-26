function Xdot = threeDimensional3BodyOrbitRotatingFrame(t,M)

    global params;
    
    % state to be integrated
    x = M(1);
    y = M(2);
    z = M(3);
    vx = M(4);
    vy = M(5);
    vz = M(6);
    
    % Distance from primary mass (Sun) and secondary mass (planet) to
    %   spacecraft - r1 and r2 respectively.
    % spacecraft - r1 and r2 respectively.
    r1 = sqrt((x+params.mu)^2 + y^2 + z^2);
    r2 = sqrt((x-1+params.mu)^2 + y^2 + z^2);
    
    % gradient of effective potential function
    U_barx = -x + (1-params.mu)*(x+params.mu)/(r1^3) + params.mu*(x-1+params.mu)/(r2^3);
    U_bary = -y + (1-params.mu)*y/(r1^3) + params.mu*y/(r2^3);
    U_barz = (1-params.mu)*z/(r1^3) + params.mu*z/(r2^3);
          
    Ubar_xx = -1 + (1-params.mu)*((r1^3-3*(x+params.mu)^2*r1)/r1^6) + params.mu*((r2^3-3*(x-1+params.mu)^2*r2)/r2^6);
    Ubar_yy = -1 + (1-params.mu)*((r1^3-3*y^2*r1)/r1^6) + params.mu*((r2^3-3*y^2*r2)/r2^6);
    Ubar_zz = (1-params.mu)*((r1^3-3*z^2*r1)/r1^6) + params.mu*((r2^3-3*z^2*r2)/r2^6);;
    Ubar_xy = -3*y*(((1-params.mu)*(x+params.mu)/r1^5) + (params.mu*(x-1+params.mu)/r2^5));
    Ubar_yx = Ubar_xy;
    Ubar_xz = -3*z*(((1-params.mu)*(x+params.mu)/r1^5) + (params.mu*(x-1+params.mu)/r2^5));
    Ubar_zx = Ubar_xz;
    Ubar_yz = -3*z*y*(((1-params.mu)/r1^5) + (params.mu/r2^5));
    Ubar_zy = Ubar_yz;
    
    % 3D CR3BP dynamics
    Xdot(1,1) = vx;
    Xdot(2,1) = vy;
    Xdot(3,1) = vz;
    Xdot(4,1) = 2*vy - U_barx;
    Xdot(5,1) = -2*vx - U_bary;
    Xdot(6,1) = -U_barz;    
end