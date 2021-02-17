function [Q,T] = SO3_cumul_BSpline_mat(n,q,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function demonstrates the cumulative B spline curves for R3 
% euclidean space with the matrix representation as explained in "efficient
% derivative computation for cumulative B-splines on Lie-Groups"
% it takes:
% n: (n-1) defines the order of of the polynomial (4 for cubic splines)
% q: an n-by-4 vector of qw,qx,qy and qz coordinates to be splined
% num: the number of uniform splined elements to be returned in eache
% interval
% and it returns:
% Q: a 4 column vector of splined coordinates
% T: the corresponding t value of the generated splines
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate the correctness of the passed variables
    validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
    validateattributes(num, {'numeric'}, {'positive','integer','scalar'});
    validateattributes(q,{'numeric'},{'real','2d'});
       
    % create the basis matrix
    k = n;
    m = zeros(n);
    M = zeros(n);
    U = zeros(n,num);
    u = linspace(0,1,num);
    Q = [];
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
    Bt = M*U; % k x pr terms
    for i = 1 : size(q,1)-n
        qb = q(i,:);
        for pr = 1 : length(u)
            for j = 1 : n-1
                logq = log_q(quatmultiply(quatinv(q(i+j-1,:)),q(i+j,:)));
                qb  = quatmultiply(qb,exp_q(Bt(j+1,pr)*logq));
            end
            Q = [Q; qb]; % B-spline segment 
            qb = q(i,:);
        end
    end 
    T = linspace(0,1,size(Q,1));

end
    
    