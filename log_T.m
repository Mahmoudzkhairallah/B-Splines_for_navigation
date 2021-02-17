function zeta = log_T(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function returns the log of transformation matrix as explained in 
% the report  "A tutorial on SE(3) transformation parameterizations
% and on-manifold optimization"
% it takes:
% Tin: an 4-by-4 transformation matrix
% and it returns:
% qout: a 6-by-1 zeta vector of the (vx,vy,vz,wx,wy,wz)
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % validate that the values passed are correct 
%     validateattributes(T,{'numeric'},{'square',4,'real','finite','nonnan'});
    
    % calculate componantes of the log vector
    R = T(1:3,1:3);
    t = T(1:3,4);
    theta = acos((trace(R)-1)/2);
    if (theta == 0)
        w = [0;0;0];
    else
        logR = theta/(2*sin(theta))*(R-R');
        w = [logR(3,2);logR(1,3);logR(2,1)];
    end
    
    if (norm(w)==0)
        V = eye(3);
    else
        V = eye(3) + ((1-cos(norm(w)))/(norm(w)^2))*logR + ((norm(w)-sin(norm(w)))/(norm(w)^3))*(logR*logR);
    end
    zeta = [V^-1*t;w];
end
