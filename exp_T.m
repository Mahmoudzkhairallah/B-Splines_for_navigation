function T = exp_T(zeta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function returns the exponential of zeta velocity vector as explained 
% in the report  "A tutorial on SE(3) transformation parameterizations
% and on-manifold optimization"
% it takes:
% Tin: an 6-by-1 zeta vector (vx,vy,vz,wx,wy,wz)
% and it returns:
% qout: a 4-y-4 transformation matrix
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % validate that the values passed are correct 
%     validateattributes(zeta,{'numeric'},{'vector',6,'real','finite','nonnan'});
    
    % calculate componantes of the log vector
    v = zeta(1:3);
    w = zeta(4:6);
    W =[0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    theta = norm(w);
    if (theta == 0)
        expW = eye(3);
        V = eye(3);
    else
        expW = eye(3) + (sin(theta)/theta)*W + ((1 - cos(theta))/theta^2)*W*W;
        V = eye(3) + ((1 - cos(theta))/theta^2)*W + ((theta - sin(theta))/theta^3)*W*W;
    end
    T = [expW, V*v; 0 0 0 1];
end
