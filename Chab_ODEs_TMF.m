function fun=Chab_ODEs_TMF(t,y,flag,x,strain_rate,temperature,temperature_rate)

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

x_funcs=x;
clear x

x=zeros(18,1);
x_rate=zeros(18,1);

% if temperature<500
%      temperature=500;
% end

for ii=1:1:18
    poly_temp=x_funcs{ii,1};
    x(ii,1)=polyval(poly_temp,temperature);
    poly_temp=polyder(poly_temp);
    x_rate(ii,1)=polyval(poly_temp,temperature);
    clear poly_temp
end


x(8)=1.89;
x(9)=-64.98;
x(11)=-4.82;


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


plas_strain_rate=max(((abs(y(3)-(y(4)+y(5))))-(y(9))-x(10)),0);

plas_strain_rate=x(5)*((sinh(plas_strain_rate/x(7)))^(1/x(6)));

plas_strain_rate=plas_strain_rate*sign(y(3)-(y(4)+y(5)));


r2=(1/x(8))*(1-exp((-x(8))*y(2)));
r2_rate=abs(plas_strain_rate)*(exp((-x(8))*y(2)));
R_dot=(x(8)*x(9)*r2_rate)+(x(11)*abs(plas_strain_rate));
R_dot=R_dot+(r2*temperature_rate*((x_rate(8)*x(9))+(x_rate(9)*x(8))));
R_dot=R_dot+(x_rate(11)*y(2)*temperature_rate);



fun=[plas_strain_rate;
    
    abs(plas_strain_rate);
    
    x(12)*(strain_rate-plas_strain_rate-(((y(3)-(x(13)*y(6)))*(1/x(14))))-(((y(3)-(x(15)*y(7)))*(1/x(16))))-((y(3)-(x(17)*y(8)))*(1/x(18))));
        
    (x(1)*plas_strain_rate)-(x(2)*y(4)*abs(plas_strain_rate))+((1/x(1))*x_rate(1)*y(4)*temperature_rate);
    
    (x(3)*plas_strain_rate)-(x(4)*y(5)*abs(plas_strain_rate))+((1/x(3))*x_rate(3)*y(5)*temperature_rate);
       
    ((y(3)-(x(13)*y(6)))*(1/x(14)));
    
    ((y(3)-(x(15)*y(7)))*(1/x(16)));
    
    ((y(3)-(x(17)*y(8)))*(1/x(18)));
    
    R_dot] ;





