function E  = exprotmT(zeta)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>
    omega   = zeta(1:3,1); % se3
    omega_s = [0,-omega(3,1),omega(2,1);omega(3,1),0,-omega(1,1);-omega(2,1),omega(1,1),0];     
    v       = zeta(4:6,1); % se3
    
    % Singularity case
    if norm(omega) == 0 
        V = eye(3);
        R = eye(3);
    else
        V = eye(3) + ((1-cos(norm(omega)))/(norm(omega)^2))*omega_s + ((norm(omega)-sin(norm(omega)))/(norm(omega)^3))*(omega_s*omega_s);
        
        
        
        R = eye(3) + (sin(norm(omega))/norm(omega))*omega_s + ((1-cos(norm(omega)))/(norm(omega)^2))*(omega_s*omega_s);
    end
    
    E = [R,V*v;0 0 0 1];
    
end