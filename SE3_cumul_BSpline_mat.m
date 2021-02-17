function [P,T] = SE3_cumul_BSpline_mat(n,p,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function demonstrates the cumulative B spline curves for R3 
% euclidean space with the matrix representation as explained in "efficient
% derivative computation for cumulative B-splines on Lie-Groups"
% it takes:
% n: (n-1) defines the order of of the polynomial (4 for cubic splines)
% p: an n-by-7 vector of qw,qx,qy and qz coordinates to be splined
% num: the number of uniform splined elements to be returned in eache
% interval
% and it returns:
% P: a 7-column vector of splined coordinates
% T: the corresponding t value of the generated splines
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate the correctness of the passed variables
    validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
    validateattributes(num, {'numeric'}, {'positive','integer','scalar'});
    validateattributes(p,{'numeric'},{'real','2d'});
       
    % create the basis matrix
    k = n;
    m = zeros(n);
    M = zeros(n);
    U = zeros(n,num);
    u = linspace(0,1,num);
    P = [];
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
    for i = 1 : size(p,1)-n
        Tb = transformation(p(i,:));
        for pr = 1 : length(u)
            for j = 1 : n-1
                tr = transformation(p(i+j-1:i+j,:));
                logT = log_T(tr(:,:,2)) ;
                Tb  = Tb*exp_T(Bt(j+1,pr)*logT);
            end
            pb = [Tb(1:3,4)',rotm2quat(Tb(1:3,1:3))];
            P = [P; pb]; % B-spline segment 
            Tb = transformation(p(i,:));
        end
    end 
    T = linspace(0,1,size(P,1));

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
    
    