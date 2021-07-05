% The desert tortoise Gopherus agassizi
% Leslie population model
% D. J. HODGSON and S. TOWNLEY: Linking management changes to population dynamic 
% responses: the transfer function of a projection matrix perturbation
% Journal of Applied Ecology 2004(41),1155-1161
% x <- Ax
function [x,y,model,trafo]=gopherus(n)
x=[individuals(rand(n,8)),zeros(n,2)]; %rand(n,2)*.05];
y=zeros(n,1);
for(i=1:n), 
    y(i)=extinction(x(i,:)); 
end
if(nargout>2)
    % try to get local functions known
    model=@extinction;
    trafo=@individuals;
end
end

function yys=extinction(x)
[n,k]=size(x);
if(k~=10), error('need 10 subpopulations'); end
yys=zeros(n,1);
for i=1:n
yy=0;
xx=x(i,1:8)';%xxs=xx;
while(sum(xx)>=2)
   
%xx=[ 0 0 0 0 0 1.3 1.98 2.57 ; .716 .567 0 0 0 0 0 0 ; ...
%    0 .194 .576 0 0 0 0 0 ; 0 0 .194 .640 0 0 0 0 ; ...
%    0 0 0 .238 .650 0 0 0 ; 0 0 0 0 .227 .678 0 0 ; ...
%    0 0 0 0 0 .249 .851 0; 0 0 0 0 0 0 .016 .860]*xx;
xx=floor([ 0 0 0 0 0 1.3 1.98 2.57 ; .716 .567 0 0 0 0 0 0 ; ...
    0 .149+x(i,9) .567 0 0 0 0 0 ; 0 0 .149+x(i,9) .604 0 0 0 0 ; ...
    0 0 0 .235+x(i,10) .560 0 0 0 ; 0 0 0 0 .225 .678 0 0 ; ...
    0 0 0 0 0 .249 .851 0; 0 0 0 0 0 0 .016 .860]*xx);
 yy=yy+1; %xxs(:,end+1)=xx;
end
yys(i)=yy;
end
%plot(0:yy,xxs,'+');hold on;set(gca,'ColorOrderIndex',1);
end

function x=individuals(u)
%    x=ceil(max( bsxfun(@plus,[46   76    35    19    13    9    15    2],20*(u-.5)),0));
% x=ceil(max( bsxfun(@plus,2*[89   163    62    27    16    13    29     5],20*(u-.5)),0));
x=floor(max( bsxfun(@plus,[89   163    62    27    16    13    29     5],20*(u)),0));
%x=ceil(max( ones(n,1)*[0   0    0    27    16    13    29     5]+20*(u-.5),0));
% x=max( ones(n,1)*[119    25     2     1     4    56    88   109]+20*(u-.5),0);
end

