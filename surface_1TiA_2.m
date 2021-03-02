function [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,y1,y2,y3,y4,y5,y6] = surface_1TiA_2(ks1,ks2)
% clc;clear ;
%     ks1 = 0.0001;
%     ks2 = 0.005;
ep0 = 8.85E-12;
ev = 1.6e-19;
e =  1.6e-19;
TiA = 0.05*ev;
TeA = 3*ev;
TiB = 0.025*ev;
TeB = 3*ev;
mi = 6.6e-26;
md = 1.15e-18;
ni0A = 5e16;
ne0A = 5e16; 
ni0B = 4e16;
ne0B = 1.35e16;
nd0 = 1.325e12;
zd = 20000; 
nn0 = 1.14e10;
nu_in = 1e4;
nu_dn = 90;
alpha_d = 0.5;
D = 0.01;
K_i = 1e-4;


%%%% Calculated values %%%%%
    tau = TiA/TeA; 
    lambda_iB = sqrt((ep0*TiB)/(ni0B*e^2));
    lambda_eA = sqrt((ep0*TeA)/(ne0A*e^2));
    lambda_eB = sqrt((ep0*TeB)/(ne0B*e^2));
    lambda_d = sqrt((lambda_iB^2*lambda_eB^2)/(lambda_iB^2+lambda_eB^2));
    w_pd = sqrt((nd0*e^2*zd^2)/(ep0*md));
    w_piA = sqrt((ni0A*e^2)/(mi*ep0));
 
    ks = linspace(ks1,ks2,20);

   
    for k = 1:length(ks) 
        
        A = -(nu_in*ks(k)^2*lambda_eA^2)-(2*ks(k)^2*lambda_eA^2*nu_in)-(3*nu_in)-(ks(k)^2*tau*alpha_d)...
            +(nn0*K_i)-(ks(k)^4*tau*alpha_d*lambda_eA^2);
        B = (lambda_eA^2*nu_in^3*ks(k)^2)+(2*w_piA^2*lambda_eA^2*nu_in*ks(k)^2)+(2*ks(k)^2*lambda_eA^2*w_piA^2*nu_in)...
            +(3*ks(k)^4*tau*alpha_d*lambda_eA^2*nu_in^2)+(2*w_piA^2*ks(k)^4*tau*alpha_d*lambda_eA^2)...
            +(w_piA^2*ks(k)^2*lambda_eA^2*nu_in)+(2*nu_in*w_piA^2*ks(k)^2*lambda_eA^2)+(nu_in^3)+(2*w_piA^2*nu_in)...
            +(2*nu_in^2*ks(k)^2*tau*alpha_d)+(2*w_piA^2*nu_in)+(ks(k)^2*tau*alpha_d*nu_in^2)...
            +(2*w_piA^2*ks(k)^2*tau*alpha_d*nu_in)-(w_piA^2*ks(k)^2*nu_in*lambda_eA^2)-(3*nn0*K_i*nu_in^2)...
            -(2*nn0*w_piA^2*K_i);
        C = -(w_piA^4*lambda_eA^2*nu_in*ks(k)^2)-(2*w_piA^2*tau*alpha_d*lambda_eA^2*nu_in^2)-(ks(k)^4*tau*alpha_d*lambda_eA^2*w_piA^4)...
            -(nu_in^3*w_piA^2*ks(k)^2*lambda_eA^2)-(2*w_piA^4*ks(k)^2*lambda_eA^2*nu_in)-(2*w_piA^4*nu_in*lambda_eA^2*ks(k)^2)...
            -(nu_in*w_piA^4)-(2*w_piA^2*nu_in^2*ks(k)^2*tau*alpha_d)-(w_piA^4*ks(k)^2*tau*alpha_d)+(w_piA^2*ks(k)^2*nu_in^3*lambda_eA^2)...
            +(2*w_piA^4*ks(k)^2*nu_in*lambda_eA^2)+(2*w_piA^2*nu_in^2*K_i*nn0)+(w_piA^4*nn0*K_i);
       
        D1 = 1+(ks(k)^2*lambda_eA^2);
        E = -(2*lambda_eA^2*nu_in^2*ks(k)^2)-(ks(k)^2*lambda_eA^2*nu_in^2)-(2*ks(k)^2*lambda_eA^2*w_piA^2)...
            -(ks(k)^4*tau*alpha_d*lambda_eA^2*nu_in)-(2*ks(k)^4*tau*alpha_d*lambda_eA^2*nu_in)-(w_piA^2*lambda_eA^2*ks(k)^2)...
            -(3*nu_in^2)-(nu_in*ks(k)^2*tau*alpha_d)-(2*w_piA^2)-(2*ks(k)^2*tau*alpha_d*nu_in)+(nn0*K_i*nu_in)...
            +(2*nn0*K_i*nu_in);
        F = (2*w_piA^2*lambda_eA^2*nu_in^2*ks(k)^2)+(ks(k)^2*lambda_eA^2*w_piA^4)+(ks(k)^4*tau*alpha_d*lambda_eA^2*nu_in^3)...
            +(2*w_piA^2*tau*alpha_d*lambda_eA^2*nu_in)+(2*w_piA^2*ks(k)^4*tau*alpha_d*lambda_eA^2*nu_in)+(2*nu_in^2*w_piA^2*ks(k)^2*lambda_eA^2)...
            +(nu_in^2*w_piA^2*lambda_eA^2*ks(k)^2)+(2*w_piA^4*ks(k)^2*lambda_eA^2)+(2*w_piA^2*nu_in^2)+(nu_in^3*ks(k)^2*tau*alpha_d)...
            +(4*w_piA^2*ks(k)^2*tau*alpha_d*nu_in)+w_piA^4 -(2*w_piA^2*ks(k)^2*nu_in^2*lambda_eA^2)-(nn0*K_i*nu_in^3)-(4*w_piA^2*nn0*K_i*nu_in);
        G = -(ks(k)^4*tau*alpha_d*lambda_eA^2*nu_in*w_piA^4)-(2*w_piA^4*ks(k)^2*lambda_eA^2*nu_in^2)-(w_piA^6*lambda_eA^2*ks(k)^2)-(nu_in*ks(k)^2*tau*alpha_d*w_piA^4)...
            +(2*w_piA^4*ks(k)^2*nu_in^2*lambda_eA^2)+(nn0*K_i*nu_in*w_piA^4);
        
        P = (lambda_d^2*ks(k)^4*D)+(D*ks(k)^2);
        Q = -(2*w_pd^2*lambda_d^2*ks(k)^4*D)-(2*w_pd^2*D*ks(k)^2)+(lambda_d^2*w_pd^2*ks(k)^2*nu_dn);
        R = (lambda_d^2*ks(k)^4*D*w_pd^4)+(D*ks(k)^2*w_pd^4)-(2*w_pd^4*nu_dn*ks(k)^2*lambda_d^2);
        S = (lambda_d^2*nu_dn*w_pd^6*ks(k)^2);
        T = -(lambda_d^2*ks(k)^2)-1;
        U = (2*w_pd^2*lambda_d^2*ks(k)^2)+(lambda_d^2*w_pd^2*ks(k)^2)+(2*w_pd^2);
        V = (lambda_d^2*w_pd^6*ks(k)^2);
        W = -(lambda_d^2*ks(k)^2*w_pd^4)-(2*w_pd^4*lambda_d^2*ks(k)^2)-w_pd^4;
        
        A1 = (lambda_d^2*D1);
        B1 = (lambda_d^2*E)+(lambda_d^2*D*ks(k)^2*A)-(D1*lambda_d^2*w_pd^2);
        C1 = (F*lambda_d^2)+(lambda_d^2*B*D*ks(k)^2)-(E*lambda_d^2*w_pd^2);
        D1 = (lambda_d^2*C*D*ks(k)^2)-(F*lambda_d^2*w_pd^2)+(lambda_d^2*G);
        E1 = -(G*lambda_d^2*w_pd^2);
        F1 = -(A*lambda_d^2)+(D1*D*ks(k)^2*lambda_d^2);
        G1 = -(lambda_d^2*B)+(E*D*ks(k)^2*lambda_d^2)+(A*lambda_d^2*w_pd^2);
        H1 = -(lambda_d^2*C)+(lambda_d^2*D*ks(k)^2*F)+(B*lambda_d^2*w_pd^2);
        I1 = (G*D*ks(k)^2*lambda_d^2)+(lambda_d^2*w_pd^2*C);
        
        X1 = -(U*lambda_eA^2)-(P*ks(k)^2*tau*alpha_d*lambda_eA^2)+(T*w_piA^2*lambda_eA^2)+(3*nu_in^2*T*lambda_eA^2)-(3*nu_in*lambda_eA^2*P)+(3*nu_in*ks(k)^2*tau*alpha_d*lambda_eA^2*T);
        Y1 = (P*nu_in^3*lambda_eA^2)-(T*ks(k)^2*tau*alpha_d*nu_in^3*lambda_eA^2)-(V*lambda_eA^2)-(Q*ks(k)^2*tau*alpha_d*lambda_eA^2)+(U*w_piA^2*lambda_eA^2)+(3*nu_in^2*U*lambda_eA^2)+(3*nu_in^2*ks(k)^2*tau*alpha_d*P*lambda_eA^2)...
             -(3*nu_in^2*T*w_piA^2*lambda_eA^2)-(3*nu_in*lambda_eA^2*Q)+(3*nu_in*ks(k)^2*tau*alpha_d*lambda_eA^2*U)+(3*P*nu_in*lambda_eA^2*w_piA^2);
        Z1 = (Q*nu_in^3*lambda_eA^2)-(ks(k)^2*tau*alpha_d*nu_in^3*U*lambda_eA^2)-(w_piA^2*nu_in^3*P*lambda_eA^2)-(W*lambda_eA^2)-(R*ks(k)^2*tau*alpha_d*lambda_eA^2)+(V*w_piA^2*lambda_eA^2)...
             +(3*nu_in^2*V*lambda_eA^2)+(3*Q*nu_in^2*ks(k)^2*tau*alpha_d*lambda_eA^2)-(3*nu_in^2*U*w_piA^2*lambda_eA^2)-(3*nu_in*lambda_eA^2*R)+(3*nu_in*lambda_eA^2*w_piA^2*Q)...
             +(3*nu_in*ks(k)^2*tau*alpha_d*lambda_eA^2*V);
        P1 = (R*nu_in^3*lambda_eA^2)-(ks(k)^2*tau*alpha_d*V*nu_in^3*lambda_eA^2)-(Q*w_piA^2*nu_in^3*lambda_eA^2)-(ks(k)^2*tau*alpha_d*S*lambda_eA^2)+(W*w_piA^2*lambda_eA^2)+(3*nu_in^2*W*lambda_eA^2)...
             +(3*nu_in^2*ks(k)^2*tau*alpha_d*R*lambda_eA^2)-(3*V*nu_in^2*w_piA^2*lambda_eA^2)-(3*nu_in*lambda_eA^2*S)+(3*nu_in*ks(k)^2*tau*alpha_d*lambda_eA^2*W)...
             +(3*nu_in*lambda_eA^2*w_piA^2*R);
        Q1 = (nu_in^3*S*lambda_eA^2)-(W*ks(k)^2*tau*alpha_d*nu_in^3*lambda_eA^2)-(w_piA^2*nu_in^3*R*lambda_eA^2)-(3*nu_in^2*w_piA^2*W*lambda_eA^2)+(3*nu_in*lambda_eA^2*w_piA^2*S)+(3*nu_in^2*ks(k)^2*tau*alpha_d*S*lambda_eA^2);
        R1 = -(w_piA^2*nu_in^3*S*lambda_eA^2);
        S1 = (P*lambda_eA^2) - (T*ks(k)^2*tau*alpha_d*lambda_eA^2)-(3*nu_in*lambda_eA^2*T);
        T1 = (T*nu_in^3*lambda_eA^2)+ (Q*lambda_eA^2)-(U*ks(k)^2*tau*alpha_d*lambda_eA^2)-(P*w_piA^2*lambda_eA^2)-(3*nu_in^2*P*lambda_eA^2)+(3*T*nu_in^2*ks(k)^2*tau*alpha_d*lambda_eA^2)-(3*nu_in*lambda_eA^2*U)...
            -(3*nu_in*ks(k)^2*tau*alpha_d*P*lambda_eA^2)+(3*nu_in*lambda_eA^2*w_piA^2*T);
        U1 = (U*nu_in^3*lambda_eA^2)+(P*ks(k)^2*tau*alpha_d*nu_in^3*lambda_eA^2)-(w_piA^2*nu_in^3*T*lambda_eA^2)+(R*lambda_eA^2)-(V*ks(k)^2*tau*alpha_d*lambda_eA^2)-(Q*w_piA^2*lambda_eA^2)-(3*nu_in^2*Q*lambda_eA^2)...
            +(3*U*nu_in^2*ks(k)^2*tau*alpha_d*lambda_eA^2)+(3*nu_in^2*w_piA^2*P*lambda_eA^2)-(3*V*nu_in*lambda_eA^2)-(3*nu_in*ks(k)^2*tau*alpha_d*lambda_eA^2*Q)...
            +(3*U*nu_in*lambda_eA^2*w_piA^2);
        V1 = (V*nu_in^3*lambda_eA^2)+(ks(k)^2*tau*alpha_d*Q*nu_in^3*lambda_eA^2)-(U*w_piA^2*nu_in^3*lambda_eA^2)+(S*lambda_eA^2) - (W*ks(k)^2*tau*alpha_d*lambda_eA^2)-(R*w_piA^2*lambda_eA^2)-(3*nu_in^2*R*lambda_eA^2)...
            +(3*V*ks(k)^2*tau*alpha_d*lambda_eA^2)+(3*nu_in^2*w_piA^2*Q*lambda_eA^2)-(3*W*nu_in*lambda_eA^2)-(3*nu_in*ks(k)^2*tau*alpha_d*lambda_eA^2*R)+(3*V*nu_in*lambda_eA^2*w_piA^2);
        eta1 =(W*nu_in^3*lambda_eA^2)+(R*tau*alpha_d*ks(k)^2*nu_in^3*lambda_eA^2)-(w_piA^2*nu_in^3*V*lambda_eA^2)-(w_piA^2*S*lambda_eA^2)-(3*nu_in^2*S*lambda_eA^2)+(3*nu_in^2*ks(k)^2*tau*alpha_d*W*lambda_eA^2)...
            +(3*nu_in^2*w_piA^2*R*lambda_eA^2)-(3*nu_in*ks(k)^2*tau*alpha_d*lambda_eA^2*S)+(3*nu_in*lambda_eA^2*w_piA^2*W);
        Xi = (ks(k)^2*tau*alpha_d*S*nu_in^3*lambda_eA^2)-(W*w_piA^2*nu_in^3*lambda_eA^2)-(3*nu_in^2*w_piA^2*S*lambda_eA^2);
        
        A2(k) = (A1+(T*lambda_eA^2));
        A3(k) = 0;
        B2(k) = (B1-X1);
        B3(k) = 0;
        C2(k) = (C1-Y1);
        C3(k) =0;
        D2(k) = (D1-Z1);
        D3(k) = 0;
        E2(k) = (E1-P1);
        E3(k) = 0;
        F2(k) = -Q1; 
        F4(k) = 0;
        F3(k) = -R1; 
        G2(k) = (F1-S1);
        H2(k) = (G1-T1);
        I2(k) = (H1-U1);
        J2(k) = (I1-V1);
        K2(k) = eta1;
        L2(k) = Xi;
       
        
        
    end
%         x1 = A2';x2 = B2'; x3 = C2'; x4 = D2'; x5 = E2' ; x6 = F2' ; 
%         x7 = F3' ;
        x1 = A2';x2 = A3'; x3 = B2'; x4 = B3'; x5 = C2' ; x6 = C3' ; x7 = D2' ;
        x8 = D3' ;x9 = E2' ;x10 = D3' ;x11 = F2' ;x12 = F4' ;x13 = F3' ;
        y1 = G2';y2 = H2';y3 = I2';y4 = J2';y5 = K2';y6 = L2';
end