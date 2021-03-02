%%%%%%%%%% Surface wave Dispersion Relation %%%%%%%%%
%%%% Electron-Ion plasma on one side and Dusty Plasma on another side %%%
%% This is the dispersion relation plot of w/w_piB vs k*lambda_e
%%% For Dust Ion Acoustic waves Data has been taken from 'Experimental 
%% quiescent drifting dusty plasmas and temporal dust acoustic wave growth'
%% For ion acoustic waves 'Laboratory observation of the dust-acoustic wave mode'

function [p2,X2] = surface_plot_w_piA_ratio(~)
ev = 1.6e-19;
e =  1.6e-19;
ep0=8.85E-12;
TiA = 0.05*ev;
TeA = 3*ev;
TiB = 0.026*ev;
TeB = 3*ev;
Td = 2*ev;
me = 9.1e-31;
mi = 6.6e-26;
md = 1.15e-12;
ni0A = 5e16;
ne0A = 5e16;
ni0B = 4e16;
ne0B = 1.35e16;
nd0 = 1.325e14;
zd = 200;


for i = 1:3
%%%% Calculated values %%%%%
    v_thiA = sqrt(TiA/mi);
    v_thiB = sqrt(TiB/mi);
    v_td   = sqrt(Td/md);
    v_theA = sqrt(TeA/me);
    v_theB = sqrt(TeB/me);
    lambda_iA = sqrt((ep0*TiA)/(ni0A*e^2));
    lambda_iB = sqrt((ep0*TiB)/(ni0B*e^2));
    lambda_eA = sqrt((ep0*TeA)/(ne0A*e^2));
    lambda_eB = sqrt((ep0*TeB)/(ne0B*e^2));
    w_pd = sqrt((nd0*e^2*zd^2)/(ep0*md));
    w_piA = sqrt((ni0A*e^2)/(mi*ep0));
    w_piB = sqrt((ni0B*e^2)/(mi*ep0));
    w_peA = sqrt((ne0A*e^2)/(me*ep0));
    w_peB = sqrt((ne0B*e^2)/(me*ep0));
    
   
    ks = linspace(0.001,100000,20);
    t = length(ks);
    p = ks.*ks;
    p1 = ks.*ks.*lambda_eB.*lambda_eB;
    p2 = ks.*lambda_eA;

    for k = 1:t 
        a = (lambda_eB^2)-(lambda_eA^2);
        b = (1/(w_piA^2))*((-(2*ks(k)^2*w_piA^2*lambda_eB^2*lambda_eA^2))-(w_piA^2*lambda_eB^2)+(2*ks(k)^2*(w_pd^2+w_piB^2)*lambda_eB^2*lambda_eA^2)+((w_pd^2+w_piB^2)*lambda_eA^2));
        c = (1/(w_piA^4))*((ks(k)^2*lambda_eA^2*lambda_eB^2*w_piA^4)-(lambda_eB^2*lambda_eA^2*ks(k)^2*((w_piB^2+w_pd^2)^2)));
        
        x1(k) = (- b + (sqrt(b^2-(4*a*c))))/(2*a);
        x2(k) = (- b - (sqrt(b^2-(4*a*c))))/(2*a);

        X1(k) = sqrt(x1(k));
        X2(k) = sqrt(x2(k));
        
        t1(k) = x1(k)/(ks(k)^2);
        if (t1 > v_td^2)
%            disp(t)
        end
    end
%     hold all
%     figure(1)
%     plot(p2,X1,'rd-');
% %     ylabel('\omega^2');
% %     xlabel('\kappa^2');
%     ylabel('\omega/\omega_{piA}');
%     xlabel('\kappa\lambda_{eA}');
%     grid on
% %     legend('n_{i0A} - 5\times10^{16}','n_{i0A} - 6\times10^{16}','n_{i0A} - 7\times10^{16}');
%     
% %     hold all
%     figure(2)
%     plot(p2,X2,'*b-');
% %       ylabel('\omega^2');
% %     xlabel('\kappa^2');
%     ylabel('\omega/\omega_{piA}');
%     xlabel('\kappa\lambda_{eA}');
%     grid on
%     legend('n_{i0A} - 1\times10^{16}','n_{i0A} - 2\times10^{16}','n_{i0A} - 3\times10^{16}');
end
end

    