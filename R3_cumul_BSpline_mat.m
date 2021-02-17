function [X,T] = R3_cumul_BSpline_mat(n,x,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function demonstrates the cumulative B spline curves for R3 
% euclidean space with the matrix representation as explained in "efficient
% derivative computation for cumulative B-splines on Lie-Groups"
% it takes:
% n: (n-1) defines the order of of the polynomial (4 for cubic splines)
% x: an n-by-3 vector of x,y and z coordinates to be splined
% num: the number of uniform splined elements to be returned in eache
% interval
% and it returns:
% X: a 3-by-m vector of splined coordinates
% T: the corresponding t value of the generated splines
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate the correctness of the passed variables
    validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
    validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
    validateattributes(x,{'numeric'},{'real','2d'});
       
    % create the basis matrix
    k = n;
    m = zeros(n);
    M = zeros(n);
    U = zeros(n,num);
    u = linspace(0,1,num);
    T = linspace(0,1,(num-1)*(size(x,2)-1));
    dx = x(:,2:end) - x(:,1:end-1);
    X = [];
    b = 0;
    for s = 0 : k-1 
        for n = 0 : k-1
            for l = s : k-1
                b = b + (-1)^(l-s)*nchoosek(k,l-s)*(k-1-l)^(k-1-n);
            end
            m(s+1,n+1) = (nchoosek(k-1,n)/(factorial(k-1))) * b; % Basis function matrix
            b = 0;
        end
    end
    % create the cumulative basis matrix
    for j = 1 : k % column
        for n = 1 : k % row
            M(j,n) = sum(m(j:k,n)); % Comulative Basis function matrix
        end
    end
    % initialize the u vector
    for i = 1:n
        U(i,:) = u.^(i-1);
    end
    % calculate the spline vlaues 
    for i = 1:size(x,2)-n
        B = [x(:,i),dx(:,i:i+k-2)]*M*U;
        X = [X B]; 
    end 
    T = linspace(0,1,size(X,2));

end
    
    