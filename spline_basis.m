function [B,x] = spline_basis(i,n,t,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function demonstrates the b-spline basis values and returns values
% of t corrsponding to the provided knot
% it takes:
% i: the index of the active interval for the splinr
% n: the order of the b-splines (not hte polynomial degree)
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
        x = linspace(t(n), t(end-n+1), 300);  % allocate points uniformly
    end
    B = spline_basis_recurrence(i,n,t,x);
    
end

function y = spline_basis_recurrence(j,n,t,x)

    y = zeros(size(x));
    if n > 1
        b = spline_basis(j,n-1,t,x);
        num = x - t(j+1);
        den = t(j+n) - t(j+1);
        if den ~= 0  % indeterminate forms 0/0 are deemed to be zero
            y = y + b.*(num./den);
        end
        b = spline_basis(j+1,n-1,t,x);
        num = t(j+n+1) - x;
        den = t(j+n+1) - t(j+1+1);
        if den ~= 0
            y = y + b.*(num./den);
        end
    elseif t(j+2) < t(end)  % treat last element of knot vector as a special case
        y(t(j+1) <= x & x < t(j+2)) = 1;
    else
        y(t(j+1) <= x) = 1;
    end
end