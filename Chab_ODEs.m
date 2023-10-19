function fun=Chab_ODEs(t,y,flag,x,strain_rate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chaboche model (2 back stress) differential equations
%
%James Rouse - 13/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Specify Chaboche material constants
%x(1)=C1;
%x(2)=gamma1;
%x(3)=C2;
%x(4)=gamma2;
%x(5)=A;
%x(6)=B;
%x(7)=D;
%x(8)=b;
%x(9)=Q;
%x(10)=k;
%x(11)=H;
% x(12)=E0;
% x(13)=E1;
% x(14)=eta1;
% x(15)=E2;
% x(16)=eta2;
% x(17)=E3;
% x(18)=eta3;


%Note:
%The following convention for y0 (the etmp state variable array, excluding R) is used
%y0(1)=plastic (inelastic) strain component (abs)
%y0(2)=accumulated plastic (inelastic) strain component (abs)
%y0(3)=Stress (MPa)
%y0(4)=Back Stress Component 1 (MPa)
%y0(5)=Back Stress Component 2 (MPa)
%y0(6)=ev1 (abs)
%y0(7)=ev2 (abs)
%y0(8)=ev3 (abs)
%
%Note R is fully deterministic from accumulated plastic strain
%and is calculated after the ODEs are solved

% if strain_rate>=0
%     x(10)=x(10);
% else
%     x(10)=x(10);
% end


plas_strain_rate=max(((abs(y(3)-(y(4)+y(5))))-((x(9)*(1-exp((-1)*x(8)*y(2))))+(x(11)*y(2)))-x(10)),0);

plas_strain_rate=x(5)*((sinh(plas_strain_rate/x(7)))^(1/x(6)));

plas_strain_rate=plas_strain_rate*sign(y(3)-(y(4)+y(5)));



fun=[plas_strain_rate;
    
    abs(plas_strain_rate);
    
    x(12)*(strain_rate-plas_strain_rate-(((y(3)-(x(13)*y(6)))*(1/x(14))))-(((y(3)-(x(15)*y(7)))*(1/x(16))))-((y(3)-(x(17)*y(8)))*(1/x(18))));
    
    (x(1)*plas_strain_rate)-(x(2)*y(4)*abs(plas_strain_rate));

    (x(3)*plas_strain_rate)-(x(4)*y(5)*abs(plas_strain_rate));
    
    ((y(3)-(x(13)*y(6)))*(1/x(14)));
    
    ((y(3)-(x(15)*y(7)))*(1/x(16)));
    
    ((y(3)-(x(17)*y(8)))*(1/x(18)))] ;





