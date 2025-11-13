% This is a one-dimesion wave propagation program with the MPM.
% It was started by Luming Shen on Mar. 16 2001.
% It was modified by Taoran Li on Nov. 12 2023.
% It was designed to calculate a bar with its left end fixed and its right end exerted.
% One material point per cell is used initially.
% For more information about the MPM, please refer to the following papers:
% D. Sulsky, Z. Chen, H.L. Schreyer, 1994, "A particle method for history-dependent materials", Computer methods in applied mechanics and engineering, vol. 118 pp.179-196.
% Z. Chen, W. Hu, L. Shen, X. Xin and R. Brannon, 2002, "An evaluation of the MPM for simulating dynamic failure with damage diffusion", Engineering Fracture Mechanics, vol. 69, pp. 1873-1890.
 
 
%% program starts here
 
close all;
clear;

NE=25; % number of elements

rho=1.; % density 
ys=100.; % Young's modulus 
A=1.;	% cross section area of the bar
L=25.;	% length of the bar

% set up the node number where a fixed boundary condition is applied.

fbc=1;

% input the external force

fbn=0.0; % magnitude of the external force
NB=NE+1; % # of node on which the external force applied

% initial velocity condition
v0=0.75; b1=pi/2/L;

% dt=0.1*0.1*L/NE/sqrt(ys/rho); % time step	
dt=0.01;
tfinal=40;%*L/sqrt(ys/rho); % total running time 
ntime=round(tfinal/dt)+1; % total running steps 

%- input nodal coordinates and connectivities

NN=round(2*NE); %total number of grid nodes 
LOC=zeros(NN,1); 
le=L/NE; % cell ( or element) size
LOC(:)=(0:NN-1)'*le; % grid node coordinates

% set initial material point coordinates

% for i=1:NE
%    xp(i)=0.25*(LOC(i)+LOC(i+1));
% end
xp=0.25:0.5:NE;

np=length(xp);

% Build global material points mass vector

MP=zeros(np,1);
for i=1:np
   MP(i)=rho*A*le/2;
end
 
% Set up initial value of Stress, Strain

SIG=zeros(np,ntime); % matrix to store stress of each point at each time step
xxp=zeros(np,ntime); % matrix to store position of each point at each time step
ssp=zeros(np,1); % stresses of points
snp=zeros(np,1); % strains of points
dssp=zeros(np,1); % stress increments of points
dsnp=zeros(np,1); % strain increments of points
mom=zeros(NN,ntime);
% Set up initial values of velocity, momentum, mass and external force 

vg=zeros(NN,1); % velocities of node
% vp=zeros(NE,1); % velocitie of points
vp = v0*sin(xp*b1);
vol0 = zeros(np,1); % volume of points
vol0(:,1) = A*le/2;
vol = vol0;
Fp = ones(np,1);
dvp=zeros(np,1); % velocity increments of points
fext=zeros(NN,1); % external forces on the nodes

%initialize time and step. 

t=0.0;
n=1;

% iteration starts

while t<=tfinal+0.000000001
   
   MVG=zeros(NN,1); % momentum at nodes
   MG=zeros(NN,1); % mass at nodes
   fint=zeros(NN,1); % internal forces at nodes

   % Find the number of cell in which the material point is located
   for i=1:np
      NC(i)=fix(xp(i)/le)+1;
      xc(i)=(xp(i)-LOC(NC(i)))/le;
   end
   
   % key part of the MPM
   for i=1:np
      MG(NC(i))=MG(NC(i))+MP(i)*(1-xc(i)); % map mass from point to left grid node
      MG(NC(i)+1)=MG(NC(i)+1)+MP(i)*xc(i); % map mass from point to right grid node
      MVG(NC(i))=MVG(NC(i))+MP(i)*vp(i)*(1-xc(i)); % map momentum from point to left grid node
      MVG(NC(i)+1)=MVG(NC(i)+1)+MP(i)*vp(i)*xc(i); % map momentum from point to right grid node
      fint(NC(i))=fint(NC(i))+vol(i)*ssp(i); % map internal force from point to left grid node
      fint(NC(i)+1)=fint(NC(i)+1)-vol(i)*ssp(i); % map internal force from point to right grid node
   end
   % if n==1
   %     fext(NB)=fbn;
   %     else
   %         fext(NB)=0;
   % end
   fext(NB)=fbn; %set force boundary condition 
   f=fint+fext; % calculate total grid node force vector
   f(fbc)=0.0; % apply the fixed boundary condition (total force=0)
   MVG=MVG+f*dt; % update the momenta at the grid nodes
   MVG(fbc)=0.0; % apply the fixed boundary condition (total momemtun=0)
   
   for i=1:np
      % map the nodal velocity increment back to the pariticle
      dvp(i)=f(NC(i))*(1-xc(i))*dt/MG(NC(i))+f(NC(i)+1)*xc(i)*dt/MG(NC(i)+1);
      
      % map the current nodal velocity back to particle
      vpbar(i)=MVG(NC(i))*(1-xc(i))/MG(NC(i))+MVG(NC(i)+1)*xc(i)/MG(NC(i)+1);
   end
   
   vp=vp+dvp'; % compute the current particle velocity
   xp=xp+vpbar*dt; % compute the current particle position
   
   % map the particle momentum back to the cell node
   % MVG=zeros(NN,1); 
   % for i=1:np
   %    MVG(NC(i))=MVG(NC(i))+MP(i)*vp(i)*(1-xc(i));
   %    MVG(NC(i)+1)=MVG(NC(i)+1)+MP(i)*vp(i)*xc(i);
   % end
   
   % find the current nodal velocity 
   for i=1:NN
      vg(NC(i))=MVG(NC(i))/MG(NC(i));
      vg(NC(i)+1)=MVG(NC(i)+1)/MG(NC(i)+1);
   end
   
   vg(fbc)=0.0; % apply the essential boundary conditions (velocity=0)
   
   % get the strain increments from gradient of nodal velocities
   for i=1:np
      dsnp(i)=(-vg(NC(i))+vg(NC(i)+1));
      Fp(i)=(1+dsnp(i)*dt)*Fp(i);
      vol(i)=Fp(i)*vol0(i);
      dsnp(i)=dt*dsnp(i);
      vel_p(n,:)=vp(i);
   end
    
   snp=snp+dsnp; % update the particle strains
   ssp=ssp+dsnp*ys; % undate the particle stresses (linear elasticity is assumed.)
            
   % save the new stress and output data
   
   SIG(:,n)=ssp; % particle stresses at new time step
   xxp(:,n)=xp'; % particle positions at new time step
   mom(:,n)=MVG;
   FPP(:,n)=Fp;
   FVOL(:,n)=vol;
   
   %update n and t
   
   n=n+1;
   t=t+dt;
end

% for output use only, may ignore them.
tt=linspace(0,tfinal,ntime);
xx=linspace(0,L,np);

% end of program   
%% Plot section
xt=0:dt:tfinal;
hold on;
plot(xt,xxp(25,:))
a1=12.5;w1=b1*sqrt(ys/rho);a2=v0/b1/L/w1;
xcm=a1+a2*sin(w1*xt);

plot(xt,xcm)
