% Mikado Plots
sobolpoints=@(n,k)net(scramble(sobolset(k),'MatousekAffineOwen'),n)
%% Ferrari-Trecate and Muselli (2002) pw. defined
kFTM=2;
trafoFTM=@(u)2*u-1;
modelFTM=@(x)(3+4*x(:,1)+2*x(:,2)).*(0.5*x(:,1)+0.29*x(:,2)>= 0).*(x(:,2)>=0)+...
            (-5-6*x(:,1)+6*x(:,2)).*(0.5*x(:,1)+0.29*x(:,2)<0).*(0.5*x(:,1)-0.29*x(:,2)<0)+...
(-2+4*x(:,1)-2*x(:,2)).* (0.5*x(:,1)-0.29*x(:,2)>= 0).*(x(:,2)<0);
%%

[Sb,Sp]=losi(kFTM,8192*2,modelFTM,trafoFTM,sobolpoints,'plot',30);
hold on;
%plot([-1,1],[0,0],'k',[-1,1],[-.5/.29,.5/.29],'k',[-1,1],[.5/.29,-.5/.29],'k')
patch([-1,0,1,1],[.5/.29,0,0,1],[.9 0.9 .9],'EdgeColor','k');
patch([-1,0,1,1],[-.5/.29,0,0,-1],[.7 .7 .7],'EdgeColor','k')
set(gca,'Children',flipud(get(gca,'Children'))); % patches to background
axis([-1,1,-1,1]);
hold off

%% a sin example (page 7, but not discussed further)
trafo=@(u)2*u-1;
%model=@(x)sin(pi*x*[1.5;1]);
model=@(x)sin(pi*x*[1;1]);
[Sb,Sp]=losi(2,8192*8,model,trafo,sobolpoints,'plot',32);
% I'm missing symmetry, how come?
%% example 6
k=2;
trafo=@(u)2*u-1;
model=@(x)sin(x(:,1))+sin(x(:,2)).*(x(:,1)<0)+sin(2*x(:,2)).*(x(:,1)>=0);
[Sb,Sp]=losi(k,8192*8,model,trafo,sobolpoints,'plot',64);
%% example 8
k=4;
trafo=@(u)u;
model=@(x)x*[1;1;1;0]+x(:,1).*x(:,2).*(x(:,4)>0.5);
[Sb,Sp]=losi(k,8192*8,model,trafo,sobolpoints,'plot',64);
%% example 10
k=2;
trafo=@(u)u+1;
model=@(x)prod(x,2);
[Sb,Sp]=losi(k,8192*8,model,trafo,sobolpoints,'plot',64);
%% log trafo (noise as of range in the title)
[Sb,Sp]=losi(k,8192*8,@(x)log10(model(x)),trafo,sobolpoints,'plot',64);
%% example 13
k=2;
% uniform in [1,2] x [1,5] such that corner point correspond to finite change
%trafo=@(u)1+u*diag([1,4]); 
% discussed with gamma
trafo=@(u)icdf('gam',u,10,1);
model=@(x)x(:,1)./(x(:,1)+x(:,2));
[Sb,Sp]=losi(k,8192*8,model,trafo,sobolpoints,'plot',32);

%% example 14
for k=2:4
    trafo=@(u)u;
    model1=@(x)1.0*any(x>.5,2);
    model2=@(x)1.0*all(x>.5,2);
    [Sb,Sp]=losi(k,8192*4,model1,trafo,sobolpoints,'plot',64)
    sum(diag(Sp))
    drawnow
    pause
    [Sb,Sp]=losi(k,8192*4,model2,trafo,sobolpoints,'plot',64)
    sum(diag(Sp))
    drawnow
    pause
end

%% example 18
 k=2;
 trafo=@(u)norminv(u,[0,1],[1,1]);
 model=@(x)exp(x(:,1)).*abs(sin(x(:,2)));
 [Sb,Sp]=losi(k,8192*32,model,trafo,sobolpoints,'plot',32);

 %% no mikado for spurious examples
 
 %% ishigami
k=3;
n=8192*32;
U=sobolpoints(n+16,2*k);U(1:16,:)=[];
model=@(x)sin(x(:,1)).*(1.0+0.1*x(:,3).^4)+ 7.0*(sin(x(:,2))).^2;
trafo=@(u)(u-.5)*2*pi;
[Sb,Sp]=losi(k,8192*32,model,trafo,U,'plot',32);
 
%%
at=wingweight();
k=length(at);
trafo=@(u)autotrans(u,at);
 [Sb,Sp]=losi(k,8192,@wingweight,trafo,@rand,'plot',32);
x1=trafo(ones(1,k)); y1=wingweight(x1);
x0=trafo(zeros(1,k)); y0=wingweight(x0);
y1-y0

% x10 offers interaction only with x1, all others interaction on the global scale,
% even x4 (but on two orders of magnitude smaller than say x8/x9)
 
%% The turtoise
% hack alert
[~,~,model,trafo]=gopherus(1);
k=8;
[Sb,Sp]=losi(k,8192*2, @(x)model([x,zeros(size(x,1),2)]),trafo,@rand,'plot',20);

