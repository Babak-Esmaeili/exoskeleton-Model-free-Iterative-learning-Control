function qddot = twoDOF_exoskeleton_dynamics_imp(q,qdot,u,dist)
        
    %% Parameters
    m1 = 2;                 % kg
    m2 = 0.85;              % kg
    
    L1 = 0.35;              % m
    L2 = 0.31;              % m
    
    L_c1 = L1/2;            % m
    L_c2 = L2/2;            % m
    
    I1 = (1/4)*m1*(L1^2);   % kg*(m^2)
    I2 = (1/4)*m2*(L2^2);   % kg*(m^2)

    g = 9.81;               % m/(s^2)
    
    %% Matrices
    M_11 = m1*(L_c1^2) + m2*((L_c2^2)+(L1^2)+2*L1*L_c2*cos(q(2))) + I1 + I2;
    M_12 = I2 + m2*((L_c2^2)+L1*L_c2*cos(q(2)));
    M_21 = M_12;
    M_22 = I2 + m2*(L_c2^2);
    
    M = [ M_11  M_12 
          M_21  M_22 ];
      
    C_11 = -m2*L1*L_c2*qdot(2)*sin(q(2));
    C_12 = -m2*L1*L_c2*(qdot(2)+qdot(1))*sin(q(2));
    C_21 = m2*L1*L_c2*qdot(1)*sin(q(2));
    C_22 = 0;
    
    C = [ C_11  C_12 
          C_21  C_22 ];
      
    G_11 = (m1*L_c2+m2*L1)*g*cos(q(1)) + m2*L_c2*g*cos(q(1)+q(2));
    G_21 = m2*L_c2*g*cos(q(1)+q(2));
    
    G = [ G_11
          G_21 ];
               
    %% Dynamics  
    qddot = (M^-1)*u - (M^-1)*C*qdot - (M^-1)*(G+dist);

end
