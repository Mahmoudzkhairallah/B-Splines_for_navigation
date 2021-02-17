function [X,T] = R3_cumul_BSpline(n,x,t,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function demonstrates the cumulative B spline curves for R3 
% euclidean space 
% it takes:
% n: (n-1) defines the order of of the polynomial (4 for cubic splines)
% x: an n-by-3 vector of x,y and z coordinates to be splined
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
    validateattributes(x,{'numeric'},{'real','2d'});

    check = t(2:end) - t(1:end-1);
    assert(all(check>= 0 ), ...
        'knots are not correctly arranged in ascending order');
    n_knots = size(x,2) + n;
    assert(n_knots == length(t), ...
        'number of knots is not correctly assigned, t = n + p')
    assert(m > size(x,2), ...
        'number of points of approximation should be more than the number of control points')
    
    % calculation of spline values
    
    T = linspace(t(n),t(end-n+1),m);
    P0 = x(:,1);
    dP = x(:,2:end) - x(:,1:end-1);
    X = zeros(size(x,1),m);
    for i = 1 :length(T)
        knot = T(i);
        B = zeros(1,size(x,2));
        if (i == 1) && (t(1)==t(2))
            X(:,i) = P0;
        else
            for j = 1:length(B)-(n-2)
                B(j) = cumul_spline_basis(j-1,n,t,knot);
            end
            X(:,i) = P0*B(1) + sum((dP.*B(2:end))')';
        end
        
        xxx = X(:,i);
    end
    
end
    
    