function [B,x] = cumul_spline_basis(i,n,t,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function demonstrates the cumulative b-spline basis values and 
% returns values of t corrsponding to the provided knot
% it takes:
% i: the index of the active interval for the spline
% n: the order of the b-splines (not the polynomial degree)
% t: the knots vector
% and it returns:
% y: the values of b-spline basis
% x: the values where b0-spline basis are evaluated
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% validate the nature of passed values
    validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
    validateattributes(t, {'numeric'}, {'real','vector'});
    check = t(2:end) - t(1:end-1);
    assert(all(check>= 0 ), ...
        'knots are not correctly arranged in ascending order');
    assert(0 <= i && i < numel(t)-n, ...
    'Invalid interval index j = %d, expected 0 =< j < %d (0 =< j < numel(t)-n).', i, numel(t)-n);
    if nargin < 4
        x =  linspace(t(n),t(end-n+1),300);  % allocate points uniformly
        B = zeros(1,300);
    end
    
    % calculation of basis value
    if length(x)>1
        for j = 1:length(x)
            b = 0;
            knot = x(j);
            if knot <= t(i+1)
                b = 0;
            elseif knot < t((i+1)+(n)-1)
                for k = i : i+n-2
                    cumul = spline_basis(k,n,t,knot);
                    b = b + cumul;
                end
            else 
                b = 1;
            end
            B(j) = b;
            
        end
    else
        b = 0;
        knot = x;
        if knot <= t(i+1)
            b = 0;
        elseif knot < t((i+1)+(n)-1)
            for k = i : i+n-2
                cumul = spline_basis(k,n,t,knot);
                b = b + cumul;
            end
        else
            b = 1;
        end
        B = b;
    
    end
end
