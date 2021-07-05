function [y] = wingweight(xx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WING WEIGHT FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT AND INPUT:
%
% y  = wing weight
% xx = [Sw, Wfw, A, LamCaps, q, lam, tc, Nz, Wdg, Wp]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin==0)
% variables and input ranges
y={ {'unif',{150,200},'wing area (ft^2)'},...
  {'unif',{220,300},'weight of fuel in the wing (lb)'},...
  {'unif',{6,11},'aspect ratio'},...
  {'unif',{-10,10},'quarter-chord sweep (degrees)'},...
  {'unif',{16,45},'dynamic pressure at cruise (lb/ft^2)'},...
  {'unif',{.5,1},'taper ratio'},...
  {'unif',{.08,.18},'aerofoil thickness to chord ratio'},...
  {'unif',{2.5,6},'ultimate load factor'},...
  {'unif',{1700,2500},'flight design gross weight (lb)'},...
  {'unif',{.025,.08},'paint weight (lb/ft^2)'}};
return
end
[n,k]=size(xx);
if(k==1)
Sw      = xx(1);
Wfw     = xx(2);
A       = xx(3);
LamCaps = xx(4);
q       = xx(5);
lam     = xx(6);
tc      = xx(7);
Nz      = xx(8);
Wdg     = xx(9);
Wp      = xx(10);
else
Sw      = xx(:,1);
Wfw     = xx(:,2);
A       = xx(:,3);
LamCaps = xx(:,4);
q       = xx(:,5);
lam     = xx(:,6);
tc      = xx(:,7);
Nz      = xx(:,8);
Wdg     = xx(:,9);
Wp      = xx(:,10);
end
fact1 = 0.036 * Sw.^0.758 .* Wfw.^0.0035;
fact2 = (A./ ((cos(LamCaps/180*pi)).^2)).^0.6;
fact3 = q.^0.006 .* lam.^0.04;
fact4 = (100.*tc ./ cos(LamCaps/180*pi)).^(-0.3);
fact5 = (Nz.*Wdg).^0.49;

term1 = Sw .* Wp;

y = fact1.*fact2.*fact3.*fact4.*fact5 + term1;

end