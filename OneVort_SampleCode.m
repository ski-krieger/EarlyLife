%This is a sample code for the manuscript "Turbulent coherent
%structures and early life below the Kolmogorov scale by Krieger, Sinai,
%Ferrari and Nowak. This code is relevant to Figure 3 and the discussion of
%replicating particles in a single-unsteady-vortex flow

%This code simulates the motion of two particles in an
%unsteady-single-vortex flow

%----------------------

%OUTPUTS:
%
% 'outcell' = the positions of both particles at each time step


%INPUTS:
%
% 'epsi' = the value of epsilon, the amplitude of the unsteady perturbation
% to a fixed one-vortex flow. 

function outcell=OneVort_SampleCode(epsi)

PopArray=[0.4,0.5;0.6,0.5];
popsize=2;
ratecount=0;
rates=0;
A=0.5; %As given in Methods/Supplementary Material

hswitch=0;
outcell=cell(10000,1); %Pre-allocating space
tcount=0;
while hswitch==0
    tcount=tcount+1;
 YYYic=[PopArray(1:popsize,1)',PopArray(1:popsize,2)',ratecount];
  opts = odeset('Events',@(t,YYY) StochEvents(t,YYY), 'RelTol',1e-2, 'InitialStep',0.001);
  [t,YYY,te,ye,ie] = ode45(@(t,YYY) ODEFUN(t,YYY),[0,0.02],YYYic,opts); 
 yout=YYY(end,:);
        PopArray(1:popsize,1)=mod(yout(1:popsize),2);
        PopArray(1:popsize,2)=mod(yout(popsize+1:end-1),1);
outcell{tcount}=PopArray;
  
if te>0 
hswitch=1;
elseif tcount>10000
    hswitch=1;
end

end


 function dydt=ODEFUN(t,y)
        dydt=zeros(2*popsize+1,1); %a column vector  
        ffun=(epsi*sin(2*pi*t))*y(1:popsize).^2+(1-2*epsi*sin(2*pi*t))*y(1:popsize);
        dffun=(2*epsi*sin(2*pi*t))*y(1:popsize)+(1-2*epsi*sin(2*pi*t));
        dydt(1:popsize)=-pi*A*sin(pi*ffun).*cos(pi*y(popsize+1:end-1));
        dydt(popsize+1:end-1)=pi*A*cos(pi*ffun).*sin(pi*y(popsize+1:end-1)).*dffun;
        dydt(end)=sum(rates);
    end


function [position,isterminal,direction] = StochEvents(t,YYY)
position = YYY(1)-1; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end


end %of OneVort_SampleCode