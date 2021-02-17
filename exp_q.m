function qout = exp_q(qin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function returns the exp of quaternion as explained in the report 
% "A tutorial on SE(3) transformation parameterizations
% and on-manifold optimization"
% it takes:
% qin: an n-by-4 vector of w,x,y and z
% and it returns:
% qout: a n-by-4 vector of the exponant quaternions
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % validate that the values passed are correct 
%     validateattributes(qin,{'numeric'},{'ncols',4,'real','finite','nonnan'});
    
    % make sure the passed quaternions are normalised 
    normth = arrayfun(@(i) norm(qin(i,2:4)),1:size(qin,1),'UniformOutput',true)';
    
    % initialize quaterion exponential
    qout = zeros(size(qin));
    
    % calculate the exponential of quaternion
    for i =1:size(qin,1)
        qout(i,:) = [cos(normth(i)), (sin(normth(i))/normth(i))*qin(i,2:4)];
        if (normth(i) ==0)
            qout(i,:) = [1 0 0 0];
        end
    end
end
