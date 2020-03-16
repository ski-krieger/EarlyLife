%This is a sample code for the manuscript "Turbulent coherent
%structures and early life below the Kolmogorov scale by Krieger, Sinai,
%Ferrari and Nowak. This code is relevant for Figure 4 in the main text,
%where replication occurs according to a Wright-Fisher rule. 

%This code simulates the following process: cooperation is according to the
%replicase R2 rule, the particles are advected by the flow arising from 7 
%point-vortices on the unit flat torus, and replication obeys a 
%Wright-Fisher update rule

%----------------------

%OUTPUTS:
%
% 'outcome' = either 1, indicating "flourishing" or 1000-fold increase from
% the initial inoculum size; or 0, indicating extinction. 
%
% 'tt' = the final time when the simulation has completed, measured in
% expected particle lifetimes


%INPUTS:
%
% 'R' = the radius of interaction/cooperation from the manuscript. R=0.03
% is a reasonable starting value. 
%
% 'popmax' = The population size for the Wright-Fisher (fixed-size) rule. 
%
% 'tflow' = The effective Damkoehler number, since now reproduction occurs
% at discrete intervals. [tflow] gives the length of time the flow code
% will run before a reproduction occurs. 
%


function [outcome,tt]=WrightFisher_SampleCode(R,popmax,tflow)

rng shuffle %So that no two realizations are identical

tt = 0; % current REAL TIME MARKER
popsize=popmax;

%Vortex parameters for 7 point-vortices
vpos=rand(7,2);
nvort=7;
ogamma=zeros(1,10000);
ogamma(1:round(nvort/2))=1;
ogamma(round(nvort/2)+1:nvort)=-1;

%Pre-allocate
PopArray=zeros(popmax,4); %This is the main matrix for the replicator population. Each row represents a (potential)
%biological particle. The number of rows should exceed 1000 times the
%initial size of the inoculum, since this is where some of the realizations
%will stop (success of the population). 

%The columns represent: the x-coordinate, y-coordinate, replicator species,
%whether or not the particle is currently catalyzed, and the particle
%lineage when applicable.

%Initialize population
for typer=1:2
for ij=(1+(popmax/2)*(typer-1)):typer*(popmax/2)
PopArray(ij,:)=[vpos(1,1)+rand()/100,vpos(1,2)+rand()/100,typer,0]; %x,y,type,catalyzed?
end
end
Zpopsize=popmax+7; %This variable tracks the number of vortices + biological particles

%%%%%%%%%%%%%
%-----------
%INITIATE PROCESS
%-----------
%%%%%%%%%%%%%

hitswitch=0;
hittime=0;
while hitswitch==0
hittime=hittime+1;

    %The dynamics proceed in three steps: 1) particles that should be
    %catalyzed are catalyzed; 2) Flow and Wright-Fisher kinetics

%1) ADJACENCY
  [idx,dist]=knnsearch(PopArray(1:popsize,1:2),PopArray(1:popsize,1:2),'k',4);
  propdist=dist(:,2:end);
  propind=idx(:,2:end);
  within=propdist<R;
  typz=reshape(PopArray(propind,3),popsize,3);
  catald=within.*typz;
  catald2=(sum(catald')>0)';
  PopArray(1:popsize,4)=0.8*ones(popsize,1)+0.7*catald2;

 %2) RATES, MOTION
  %ODE structure: vort-x, people-x, vort-y, people-y, rates
  YYYic=[vpos(:,1)',PopArray(1:popsize,1)',vpos(:,2)',PopArray(1:popsize,2)'];
  %opts = odeset('Events',@(t,YYY) StochEvents(t,YYY,targettime), 'RelTol',1e-2, 'AbsTol',1e-2,'InitialStep',0.001/sum(rates));
  [t,YYY] = ode23(@(t,YYY) ODEFUN(t,YYY),[0,tflow],YYYic); 
  
  tt=tt+tflow;
  yout=YYY(end,:);
  PopArray(1:popsize,1)=mod(yout(8:popsize+7),1);
  PopArray(1:popsize,2)=mod(yout(popsize+15:end),1); 
  vpos(:,1)=mod(yout(1:7),1);
  vpos(:,2)=mod(yout(Zpopsize+1:Zpopsize+7),1);
  
  
  %3) REPRODUCTION
  Payoff=PopArray(:,4);
  parents=discretedist(Payoff',popmax);
  NextGen=PopArray;
  for ppl=1:popmax
      NextGen(ppl,:)=[mod(PopArray(parents(ppl),1)+rand()/100,1),mod(PopArray(parents(ppl),2)+rand()/100,1),PopArray(parents(ppl),3),0];
  end
  PopArray=NextGen;
         
  
  
  %End of game check
if isempty(find(PopArray(:,3)==1))
    hitswitch=1;
    outcome=0; %no hypercycle
elseif isempty(find(PopArray(:,3)==2))
    hitswitch=1;
    outcome=1; 
end
  
  
end %of while loop
   
  
  %%%%%%%%%%%%%%%
%AUX FUNCTIONS


function dydt=ODEFUN(t,y) %The set of ODE's that are constantly being integrated. 
        dydt=zeros(2*Zpopsize,1); %a column vector  
        xposes=y(1:Zpopsize);
        yposes=y(Zpopsize+1:end);
        dvx2 = bsxfun(@minus,xposes,xposes');
        dvy2 = bsxfun(@minus,yposes,yposes');
       
        mxm2=cosh(dvx2+4*pi);
        mxm1=cosh(dvx2+2*pi);
        mx0=cosh(dvx2);
        mx1=cosh(dvx2-2*pi);
        mx2=cosh(dvx2-4*pi);
        mym2=cosh(dvy2+4*pi);
        mym1=cosh(dvy2+2*pi);
        my0=cosh(dvy2);
        my1=cosh(dvy2-2*pi);
        my2=cosh(dvy2-4*pi);

         zmatx=-(sin(dvy2)./(mxm2-cos(dvy2)))-(sin(dvy2)./(mxm1-cos(dvy2)))-(sin(dvy2)./(mx0-cos(dvy2)))-(sin(dvy2)./(mx1-cos(dvy2)))-(sin(dvy2)./(mx2-cos(dvy2)));
         zmaty=sin(dvx2)./(mym2-cos(dvx2))+sin(dvx2)./(mym1-cos(dvx2))+sin(dvx2)./(my0-cos(dvx2))+sin(dvx2)./(my1-cos(dvx2))+sin(dvx2)./(my2-cos(dvx2));
         zmatx(isnan(zmatx))=0;
         zmaty(isnan(zmaty))=0;
         dydt(1:Zpopsize)=zmatx*ogamma(1:Zpopsize)';
         dydt(Zpopsize+1:end)=zmaty*ogamma(1:Zpopsize)';
 end %of odes

function BB = discretedist(probs,NNN)
%normalize probs
Pnorm=[0 P]/sum(probs);

%create cumlative distribution
Pcum=cumsum(Pnorm);
R=rand(1,NNN);
V=1:length(P);
[~,inds] = histc(R,Pcum); 
BB = V(inds);
BB=reshape(BB,NNN,1);
end


end %of R2WrightFisher