%This is a universal sample code for the manuscript "Turbulent coherent
%structures and early life below the Kolmogorov scale by Krieger, Sinai,
%Ferrari and Nowak. Provided the code is run enough times, it can replicate
%most all the main results in the paper, and requires minimal changes to
%extend to any relevant result. However, due to the extremely slow run
%times, it will not be feasible to replicate the EXACT results on a
%personal computer; indeed the results shown in the manuscript took many
%months on computer clusters to complete.

%This code simulates the following process: starting from a small inoculum
%(here the size of the inoculum has been slightly increased from the value
%used in the manuscript, to lessen the number of extinction events) of
%replicators obeying a metabolism chosen by the user, the particles are
%advected by the flow arising from 7 point-vortices on the unit flat torus
%and obey a standard time-dependent Gillespie process for birth and death
%events, with rates determined by the metabolism. 

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
% 'nstop' = A parameter that can control how many expected particle lifetimes the code goes through
%without stopping. A useful debugging / time-control argument. 
%
% 'Damkoh' = The *INVERSE* of the Damkohler number as defined in the text
%
% 'metabtype' = The metabolism type; 1 = Replicase R1, 2 = Replicase R2, 3
% = hypercycle


function [outcome,tt]=EarlyLife_SampleCode(R,nstop,Damkoh,metabtype)

rng shuffle %So that no two realizations are identical

ntype=2; %The number of types -- values greater than 2 are only relevant for the hypercycle. For replicases, 
%values greater than 2 simply mean there are different "types" of
%defectors, but their behavior is the same. 

nline=1; %The number of cell "lineages", as though flourescently marked -- this is relevant for Figs. 6 and S2

tt = 0; % current REAL TIME MARKER, measured in expected particle lifetimes
R_counter = 0; % reaction counter for Gillespie process
loopcount=0; %helpful to count how many reactions / flow breaks have occurred
hitswitch=0; %This will tell us when the simulation is over
reacswitch=1; %This controls whether or not to perform a stochastic reaction according to the Gillespie process

%VORTEX PARAMETERS
nvort=7; %Number of vortices; can easily be made dynamic / passed as an argument
vpos=rand(nvort,2); %The vortex positions
ogamma=zeros(1,10000); %Matrix of vortex strengths. Note that there will be some "pseudo-vortices" with zero net circulation. 
ogamma(1:round(nvort/2))=Damkoh; %Vortex strength, here equal to the Damkohler number
ogamma(round(nvort/2)+1:nvort)=-Damkoh; %Half the vortices have the opposite strength

%GILLESPIE PARAMETERS
R_max = 10^7; %The maximum number of reactions considered; this number should just be sufficiently large to
%never be exceeded. It exists because it dramatically speeds up the code to
%pre-allocate Gillespie hitting times
rand_nums = rand(2,R_max); % TWO random numbers for each reaction step - one from the uniform distribution
%to draw the probabilities of particular events, when a hitting time has
%been reached; and one from the exponential distribution to determine the
%hitting times
rand_nums(1,:)=exprnd(1,R_max,1);

%Pre-allocate
PopArray=zeros(10000,5); %This is the main matrix for the replicator population. Each row represents a (potential)
%biological particle. The number of rows should exceed 1000 times the
%initial size of the inoculum, since this is where some of the realizations
%will stop (success of the population). 

%The columns represent: the x-coordinate, y-coordinate, replicator species,
%whether or not the particle is currently catalyzed, and the particle
%lineage when applicable.

%Initialize population
for lineage=1:nline
for typer=1:ntype
for ij=(1+10*(typer-1))+(lineage-1)*ntype*10:10+(10*(typer-1))+(lineage-1)*ntype*20 
PopArray(ij,:)=[rand(),rand(),typer,0,lineage]; %x,y,type,catalyzed?,lineage
end
end
end
popsize=10*ntype*nline;
Zpopsize=popsize+nvort; %This variable tracks the number of vortices + biological particles

%%%%%%%%%%%%%
%-----------
%INITIATE PROCESS
%-----------
%%%%%%%%%%%%%

ratecount=0;

while hitswitch==0

loopcount=loopcount+1;
    if reacswitch==1 %Bookkeeping for stochastic reactions
        reacswitch=0;
        ratecount=0; %reset rate counter
R_counter = R_counter + 1;
targettime=rand_nums(1,R_counter);
    end

    %The dynamics proceed in three steps: 1) particles that should be
    %catalyzed are catalyzed; 2) Flow and Gillespie kinetics; 3) When a stochastic event occurs,
    %choose birth or death, and then do some population bookkeeping
  
  %1a) Find, for each particle, the other particles within a radius R
 [idx,dist]=knnsearch(PopArray(1:popsize,1:2),PopArray(1:popsize,1:2),'k',ntype+2); 
  propdist=dist(:,2:end);
  propind=idx(:,2:end);
  within=propdist<R;
  
  %1b) Cataylze particles according to the metabolism  
  if metabtype==3 %Hypercycle metabolism
    typz=reshape(PopArray(propind,3),popsize,ntype+1); %The following are all bookkeeping commands to find if the correct lineage to catalyze is nearby
    typz2=mod(typz,ntype)+ones(popsize,ntype+1);
     metype=repmat(PopArray(1:popsize,3),1,ntype+1);
     typz3=(typz2==metype);
     catald=within.*typz3;
     catald2=(sum(catald')>0)';
    PopArray(1:popsize,4)=0.8*ones(popsize,1)+0.7*catald2; %Basal birth rate plus catalyzing bonus
  else %Replicase metabolism
  typz=reshape(PopArray(propind,3),popsize,3);
  catald=within.*typz;
  catald2=(sum(catald')>0)';    
  if metabtype==1 %In the Replicase R1 metabolism (which flourishes regardless of spatial structure), replicases can copy themselves
     selfreps=mod(PopArray(1:popsize,3),2);
     catald2=(catald2+selfreps)>0;
  end 
  PopArray(1:popsize,4)=0.8*ones(popsize,1)+0.7*catald2; %Basal birth rate plus catalyzing bonus
  end

 
  

  %2) Fluid flow and updating stochastic kinetics
  rates=[popsize, sum(PopArray(1:popsize,4))]; %Since the death rate is normalized to 1, the relative death
  %rate in the whole population is 1*the population size; the second column
  %is the population birth rate. 
  
  %ODE structure to pass to the ODE solver: vort-x, people-x, vort-y, people-y, rates
  YYYic=[vpos(:,1)',PopArray(1:popsize,1)',vpos(:,2)',PopArray(1:popsize,2)',ratecount];
  opts = odeset('Events',@(t,YYY) StochEvents(t,YYY,targettime), 'RelTol',1e-2, 'AbsTol',1e-2,'InitialStep',0.001/sum(rates));
  [t,YYY,te,ye,ie] = ode23(@(t,YYY) ODEFUN(t,YYY),[0,0.01],YYYic,opts); 

   if ~isempty(te) %If a stochastic reaction will occur, meaning a hitting time has been reached
        reacswitch=1; %Important -- this means that section (3) of the dynamics will take place
        tt=tt+te;
        yez=ye(end,:);
        vpos(:,1)=mod(yez(1:nvort),1);
        vpos(:,2)=mod(yez(Zpopsize+1:Zpopsize+nvort),1);
        PopArray(1:popsize,1)=mod(yez(nvort+1:popsize+nvort),1);
        PopArray(1:popsize,2)=mod(yez(popsize+(2*nvort)+1:end-1),1);     
   else %If a stochastic event did not occur within a reasonable time window (here, 1/100 of an expected particle lifetime)
        tt=tt+0.01;
        yout=YYY(end,:);
        PopArray(1:popsize,1)=mod(yout(nvort+1:popsize+nvort),1);
        PopArray(1:popsize,2)=mod(yout(popsize+(2*nvort)+1:end-1),1); 
        vpos(:,1)=mod(yout(1:nvort),1);
        vpos(:,2)=mod(yout(Zpopsize+1:Zpopsize+nvort),1);
        ratecount=yout(end);
   end
   
   
   %3) Birth-death events
   if reacswitch==1 %meaning a hitting time has been reached
valid_rates = [ones(popsize,1); PopArray(1:popsize,4)]; 
  
selection_intvl = cumsum(valid_rates); %Make a weighted vector of all the weights
selection_intvl = selection_intvl./selection_intvl(end); %normalize
selected_ind = find(selection_intvl>rand_nums(2,R_counter),1,'first'); %pick the reaction to occur 
%uniformly-at-random given the relative propensities

if selected_ind > popsize %birth event
    PopArray(popsize+1,:)=PopArray(selected_ind-popsize,:)+[rand()/30,rand()/30,0,0,0];
else %death event
    PopArray(selected_ind,:)=[0,0,0,0,0]; %Remove particle from the living population
end

 %Re-sort population
PopArray=sortrows(PopArray,-1);
popsize=find(PopArray(:,3)==0, 1, 'first')-1;
Zpopsize=popsize+nvort;


    end % of reaction
    
   
    
  %Check if the process has reached one of the two absorbing states
if popsize<=ntype+1
    hitswitch=1;
    outcome=0; %Population collapse
elseif popsize>5000
    hitswitch=1;
    outcome=1; %Population success
elseif tt>nstop %This can be a helpful switch to keep simulations from running too long, if one wants to
    %change parameters this can affect the runtime dramatically and so it
    %can be nice to set a maximum number of reactions/flow iterations to do
    hitswitch=1;
    outcome='nstop exceeded';
    PopArray=PopArray(1:popsize,:);
end

end %of while loop


%%%%%%%%%%%%%%%
%AUX FUNCTIONS


function dydt=ODEFUN(t,y) %The set of ODE's that are constantly being integrated. 
        dydt=zeros(2*Zpopsize+1,1); %a column vector  
        dydt(end)=sum(rates);  %Except for this term, everything else is solving Eqns. S3-S4 from the manuscript,
        %e.g., the flow arising from a system of 7 point-vortices on the
        %unit flat torus. This term simply adds the rates for the
        %birth-death process at the same time, effectively calculating the
        %next hitting time for an event. 
        xposes=y(1:Zpopsize);
        yposes=y(Zpopsize+1:end-1);
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
         dydt(Zpopsize+1:end-1)=zmaty*ogamma(1:Zpopsize)';
 end %of odes


%This is a standard trigger to exit integration in MATLAB. In this case, it
%halts integration when the stopping time for the Gillespie process has
%been reached. 
function [position,isterminal,direction] = StochEvents(t,YYY,targettime)
position = YYY(end)-targettime; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end



end %of EarlyLife_SampleCode.m