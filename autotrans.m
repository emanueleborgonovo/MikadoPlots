function x=autotrans(u,distribparamlist)
%AUTOTRANS Transform unit hypercube using the given parameter list
%dp={ { 'unif', { 3 ,4 }}, {'norm',{ 1}}, {'beta',{1,7}}}
%x=autotrans(rand(n,3),dp);
[n,k]=size(u);
if(length(distribparamlist)~=k), error('Arguments mismatch.\n'); end
x=zeros(n,k);
for i=1:k
    dp=distribparamlist{i};
    distr=dp{1}; %% should check with exist()==2
    if(isa(distr,'function_handle'))
        x(:,i)=distr(u(:,i));
    else
        p=dp{2};
        x(:,i)=feval([ distr 'inv'],u(:,i), p{:});
    end
end
end
