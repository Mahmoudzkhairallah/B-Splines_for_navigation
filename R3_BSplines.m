function [X,T] = R3_BSplines(n,x,t,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function demonstrates the B spline curves for R3 euclidean space 
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
    knots = linspace(t(n),t(end-n+1),m);


    % prepare for the for loop to get the basis function
 
    deg = n-1;
    mul = sum((repmat(knots',1,length(t)) == repmat(t,length(knots),1))')';
    acive_knots = (repmat(knots',1,length(t)) >= repmat(t,length(knots),1))...
        &(repmat(knots',1,length(t)) < repmat([t(2:end) 2*t(end)],length(knots),1));
    [row,col] = find(acive_knots);
    [~,ind] = sort(row); 
    index = col(ind);
    rec_control = zeros(size(x,1),n,n); % container for recursion values
    alpha = zeros(n,n); % container for basis of recursion
    X = zeros(size(x,1), numel(knots)); % the approximated points variable to be returned

    % the for loop to get the splines
    for j = 1 : length(knots)
        num = mul(j);
        ind = index(j);
        knot = knots(j);
        rec_control(:,(ind-deg):(ind-num),1) = x(:,(ind-deg):(ind-num)); % fill n degree control points
        h = deg - num; % to check from where we can start the approximation
        if (h>0)
            % start the recursion process
            for r = 1:h
                l = ind-1;
                for i = l-deg+r:l-num % the active interval
                    alpha(i+1,r+1) = (knot-t(i+1)) / (t((i+1)+deg-r+1)-t(i+1));
                    rec_control(:,i+1,r+1) = (1-alpha(i+1,r+1)) * rec_control(:,i,r)...
                        + alpha(i+1,r+1) * rec_control(:,i+1,r);
                end
            end
           X(:,j) =rec_control(:,ind-num,deg-num+1);
        elseif ind == numel(t)
            X(:,j) = x(:,end); % the final point
        else
            X(:,j) = x(:,ind-deg);
        end

    end
end