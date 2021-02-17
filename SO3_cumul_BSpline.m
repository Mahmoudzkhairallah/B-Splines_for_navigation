function [Q,T] = SO3_cumul_BSpline(n,q,t,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function demonstrates the cumulative B spline curves for SO3 
% euclidean space (namely the quaternions) 
% it takes:
% n: (n-1) defines the order of of the polynomial (4 for cubic splines)
% q: an n-by-4 vector of w,x,y and z coordinates to be splined
% t: defines the knots needed for the splines
% m: the number of uniform splined elements to be returned
% and it returns:
% X: a 3-by-m vector of splined coordinates
% T: the corresponding t value of the generated splines
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate the correctness of the passed variables
    validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
    validateattributes(t, {'numeric'}, {'real','vector'});
    validateattributes(q,{'numeric'},{'ncols',4,'real','finite','nonnan'});
    check = t(2:end) - t(1:end-1);
    assert(all(check>= 0 ), ...
        'knots are not correctly arranged in ascending order');
    n_knots = size(q,1) + n;
    assert(n_knots == length(t), ...
        'number of knots is not correctly assigned, t = n + p')
    assert(m > size(q,1), ...
        'number of points of approximation should be more than the number of control points')
    
    % calculation of spline values
    
    T = linspace(t(n),t(end-n+1),m);
    q0 = q(1,:);
    w = [q0; quatmultiply(quatinv(q(1:end-1,:)),q(2:end,:))];
    Q = [ones(m,1),zeros(m,3)];
    for i = 1 :length(T)
        knot = T(i);
        B = zeros(1,size(q,1));
        if (i == 1) && (t(1)==t(2))
            Q(i,:) = q(i,:);
        else
            for j = 1:length(B)-(n-2)
                B(j) = cumul_spline_basis(j-1,n,t,knot);
            end
            for j = 1:length(B)
                Q(i,:) = quatmultiply(Q(i,:),exp_q(B(j)*log_q(w(j,:))));
            end
        end
    end
    
end
