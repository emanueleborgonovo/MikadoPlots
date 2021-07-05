function [Sb,Sp,Vy]=mikado(k,n,model,trafo,randomsource,interactions,num_lines)
% MIKADO Sensitivity effects using Liu Owen, with Mikado plot.
%       [SB,SP,V]=LOSI(K,N,MODEL,TRAFO) computes subset and superset 
%		sensitivity for MODEL taking K independent parameters with 
%       distributions given via TRAFO and for N base samples, hence using
%       N*(K+2+K*(K-1)/2) model evaluations, V is output variance
%       [SB,SP]=LOSI(K,N,MODEL,TRAFO,RSOURCE) uses RSOURCE as random
%       number source which is either a function taking 2 arguments
%       or a N x 2K matrix, e.g.
%       [SB,SP]=LOSI(K,N,MODEL,TRAFO,@SOBOLSEQ) gives Sobol' algorithm
%       [SB,SP]=LOSI(K,N,MODEL,TRAFO,[],'plot',L) enables Mikado plots
%               using upto L lines
%
% References: Liu Owen JASA 2006, Fruth et al 2012
%  For graphics: Borgonovo, Rabitti, Plischke 2021
 if (nargin<=5) || isempty(interactions)
  interactions='none'; % filename or 'plot' or 'none'
 else
% Max. Number of Interaction Mikado Lines
  if(nargin<=6) || isempty(num_lines)
   S=30;
  else
   S=num_lines;
  end
 end
% create 2 sets of samples
 if(nargin<=4) || isempty(randomsource)
 % use standard uniform random generator
  ua=rand(n,k);
  ub=rand(n,k);    
 else
 % a function handle?
 % @(n,k)net(scramble(sobolset(k),'MatousekAffineOwen'),n)
 % gives randomized quasi Monte Carlo
  if ~isnumeric(randomsource)
   rsource=randomsource(n,2*k);
   ua=rsource(:,1:k);
   ub=rsource(:,(k+1):end);  
  else
  % a matrix of size n x 2k ?
   [nn,kk]=size(randomsource);
   if(n~=nn) || (k~=kk/2)
    error('LOSI: incompatible size of randomsource');
   end
   ua=randomsource(:,1:k);
   ub=randomsource(:,(k+1):end);
  end
 end
% use the provided transformation or identity mapping
 if(nargin<=3) || isempty(trafo)
  xa=ua;xb=ub;
 else
  xa=trafo(ua);xb=trafo(ub);
 end
%% model evaluation for reference sample runs
 ya=model(xa);
 yb=model(xb);

%% estimating mean and var
 Ey=(mean(ya)+mean(yb))/2;
 za=ya-Ey;zb=yb-Ey;
 Vy=(za'*za+zb'*zb)/(2*n-1);

% subset and superset importance
 Sb=zeros(k,k);Sp=zeros(k,k);
 yi=zeros(n,k); % store mixture outputs
% prepare graphics
 try
 % version 2020
  tiledlayout('float'); % automatic layout
 catch
  pp=1; % counter for subplots    
 end
 
 signx=@(x)sign(x)+(x==0); % sign function with -1 or 1, no zero
 for i = 1:k
  % replace ith parameter from run a with the one from run b
  xi=xa;xi(:,i)=xb(:,i);
  % evaluate model
  yi(:,i)=model(xi);
  % difference
  dy=yi(:,i)-ya;
  % Sobol' fixing of unessential factors
  Sb(i,i)=(yb'*dy)/(n*Vy);
  % for totals, use Jansen's formula
  Sp(i,i)=(dy'*dy)/(2*n*Vy);
 % now, process second order effects
  for j=1:i-1
   xij=xi;xij(:,j)=xb(:,j);
   yij=model(xij);
   d2y=yij-yi(:,j)-dy;
   Sb(i,j)=(yb'*(yij-ya))/(n*Vy);
   Sp(i,j)=(d2y'*d2y)/(4*n*Vy); % Liu Owen superset importance measure
        
   switch(lower(interactions))
   case 'none'
    % ignore
   case 'plot'
    % Perform Interaction Mikado
    z=d2y./signx((xb(:,i)-xa(:,i)).*(xb(:,j)-xa(:,j)));
           [dz,iz]=sort(abs(z),'descend');

    % suppress noise
    ii=iz(dz(1:S)>1e-6*sqrt(Vy));
    % ii=iz(1:S); % take all
    if ~isempty(ii) && abs(z(ii(1)))>0 %0.000001*abs(ya-yb)% if (abs(z(ii(1)))/abs(ya-yb)>0.0001)
     jj=ii(z(ii)>0); % split into positive and negative
     kk=ii(z(ii)<0);
     try
	  nexttile();
     catch
      if(mod(k,2)==0)
       subplot(k/2,k-1,pp);
      else
       subplot((k-1)/2,k,pp);
      end
      pp=pp+1;
     end		
     % swap i and j such that the smaller index is at bottom
     plot([xa(jj,j),xb(jj,j)]',[xa(jj,i),xb(jj,i)]','b-',...
             [xa(kk,j),xb(kk,j)]',[xa(kk,i),xb(kk,i)]','r--');
     xlabel(['x_{' num2str(j) '}']);
     ylabel(['x_{' num2str(i) '}']);
     title(['Magnitude: ' ...
             num2str(abs(z(ii(end))),2) ' to ' num2str(abs(z(ii(1))),2) ]);
     set(gca,'FontSize',14)
    end   
   otherwise
    deltay{i}{j}=dy2./(1-2.*((xb(:,i)-xa(:,i)).*(xb(:,j)-xa(:,j))>=0));
    deltay2{i}{j}=dy2./((xb(:,i)-xa(:,i)).*(xb(:,j)-xa(:,j)));
   end
  end
 end
% save data into a file
 if ~any(strcmpi(interactions,{'none','plot'}))
  save([ interactions '.dat'],'xa','xb','ya','yb','deltay','deltay2','-mat');
 end
end