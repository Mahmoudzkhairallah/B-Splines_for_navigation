function [P,Tout,T] = SE3_cumul_BSpline(n,p,t,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function demonstrates the cumulative B spline curves for SE3 
% euclidean space (namely pose coordinates "position and orientation") 
% it takes:
% n: (n-1) defines the order of of the polynomial (4 for cubic splines)
% p: an n-by-7 vector of rx,ry,rz,qw,qx,qy and qz coordinates to be splined
% t: defines the knots needed for the splines
% m: the number of uniform splined elements to be returned
% and it returns:
% P: a m-by-7 vector of splined coordinates
% Tout: the set of transformations that correspond to each splined P
% T: the corresponding t value of the generated splines
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate the correctness of the passed variables
    validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
    validateattributes(t, {'numeric'}, {'real','vector'});
    validateattributes(p,{'numeric'},{'ncols',7,'real','finite','nonnan'});
    check = t(2:end) - t(1:end-1);
    assert(all(check>= 0 ), ...
        'knots are not correctly arranged in ascending order');
    n_knots = size(p,1) + n;
    assert(n_knots == length(t), ...
        'number of knots is not correctly assigned, t = n + p')
    assert(m > size(p,1), ...
        'number of points of approximation should be more than the number of control points')
    
    % calculation of spline values
%     
    T = linspace(t(n),t(end-n+1),m);
    Tr = transformation(p);
    Tout = zeros(4,4,m);
    P = zeros(m,7);
    for i = 1 :length(T)
        knot = T(i);
        B = zeros(1,size(p,1));
        if (i == 1) && (t(1)==t(2))
            Tout(:,:,1) = Tr(:,:,1);
            P(i,:) = p(1,:);
        else
            for j = 1:length(B)-(n-2)
                B(j) = cumul_spline_basis(j-1,n,t,knot);
            end
            Tout(:,:,i) = eye(4);
            for j = 1:length(B)
                Tout(:,:,i) = Tout(:,:,i)*exp_T(B(j)*log_T(Tr(:,:,j))); 
            end
            P(i,:) = [Tout(1:3,4,i)',rotm2quat(Tout(1:3,1:3,i))];
        end
    end
    
end

function Tr = transformation(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function creates the transformation matrix Tr between each given
% pose P while putting first pose as it is with no subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    len = size(P,1);
    Tr = zeros(4,4,len);
    Tr(:,:,1) = [quat2rotm(P(1,4:7)),P(1,1:3)'; 0 0 0 1];
    for i = 2:len
        Ti = [quat2rotm(P(i,4:7)), P(i,1:3)';0 0 0 1];
        Ti_1 = [quat2rotm(P(i-1,4:7)), P(i-1,1:3)';0 0 0 1];
        Tr(:,:,i) = inv_T(Ti_1)*Ti;
    end
end

function Tout = inv_T(Tin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function returns the inverse of a transformation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tout = [Tin(1:3,1:3)', -Tin(1:3,1:3)'*Tin(1:3,4); 0 0 0 1];
end
