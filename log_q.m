function qout = log_q(qin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function returns the log of quaternion as explained in the report 
% "A tutorial on SE(3) transformation parameterizations
% and on-manifold optimization"
% it takes:
% qin: an n-by-4 vector of w,x,y and z
% and it returns:
% qout: a n-by-4 vector of the ogarithmic quaternions
% author: Mahmoud Z. KHAIRALLAH
% mail: mahmoud.khairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % validate that the values passed are correct 
%     validateattributes(qin,{'numeric'},{'ncols',4,'real','finite','nonnan'});
    
    % make sure the passed quaternions are normalised 
    qin = quatnormalize(qin);
    normth = arrayfun(@(i) norm(qin(i,2:4)),1:size(qin,1),'UniformOutput',true)';
    th = acos(qin(:,1));
    
    % initialize quaterion log
    qout = zeros(size(qin));
    
    % calculate the log of quaternion
    for i =1:size(qin,1)
        qout(i,2:4) = th(i)/normth(i,:) * qin(i,2:4);   
    end
end
