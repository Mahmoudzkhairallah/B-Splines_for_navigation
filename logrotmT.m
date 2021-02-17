function D  = logrotmT(T)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>
    R     = T(1:3,1:3);
    t     = T(1:3,4);
    theta = acos((trace(R)-1)/2);    
    d     = (theta/(2*sin(theta)))*(R-R');
    omega = [d(3,2);d(1,3);d(2,1)];
    
    % Singularity case
    if norm(omega) == 0 
        V = eye(3);
    else
        V = eye(3) + ((1-cos(norm(omega)))/(norm(omega)^2))*d + ((norm(omega)-sin(norm(omega)))/(norm(omega)^3))*(d*d);
    end
    
    D = [omega;(V^-1)*t];
    
end