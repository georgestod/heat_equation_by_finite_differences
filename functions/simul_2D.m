function [time,S]=simul_2D(beta, simulation_parameters,data)
%
% perfoms the actual computations based on simulation parameters
%
% remark: the equations are casted in matrix form and matlab's sparsity representation is used for faster computations
%
% by Georges Tod, 2019
% georges.tod@outlook.com
%
%
% physical parameters (betas)
% 1: material conductivity
% -: density*specific heat
% 2: length in z dimension

beta2=1190*1250;
lz=beta(2); %length in z dimesion (in mm) below top layer to reach ambient temperature
tlaser=beta(4);

% unpacking simulation parameters and data
dx=simulation_parameters.dx;
xmax=simulation_parameters.xmax;
ymax=simulation_parameters.ymax;

liquid_temp=10;
M0=liquid_temp*ones(xmax,ymax)+273.15;    % visible layer is overwritten but could be loaded from previous state

for x=1:xmax
    for y=1:ymax
        ic(1).mat(x,y)= M0(x,y);
    end
end

% timestep to guarantee stability of forward Euler method
texp_max=max(data.timeEXP);
dy=dx;
ts=dx^2/(4*beta(1)/beta2);
time=[0:ts:texp_max+tlaser];


% for the part geometry
n=1;

for x=1:xmax
    for y=1:ymax
        if data.model(x,y)==1
            part.coords(n,1)=x;
            part.coords(n,2)=y;

            n=n+1;
        end
    end
end
clear x, clear y

% for the computation zone %%%%%%%%%%%%% here we are implicitely inferring boundary conditions
n=1;
for x=2:xmax-1
    for y=2:ymax-1
        comp.coords(n,1)=x;
        comp.coords(n,2)=y;
        n=n+1;
    end
end
clear x, clear y


%% Into state space form
% Determining A and B matrices + initial conditions
A=-4*eye(length(comp.coords));
B=zeros(length(comp.coords),1);

for i=1:length(comp.coords)
    x=comp.coords(i,1);
    y=comp.coords(i,2);
    
    % west
    [aw]=find(comp.coords(:,1)==x & comp.coords(:,2)==y-1,1);
    % north
    [an]=find(comp.coords(:,1)==x-1 & comp.coords(:,2)==y,1);
    % east
    [ae]=find(comp.coords(:,1)==x & comp.coords(:,2)==y+1,1);
    % south
    [as]=find(comp.coords(:,1)==x+1 & comp.coords(:,2)==y,1);
    
    % if point is a state write in A, else write in B
    if ~isempty(aw)
        A(i,aw)=1;
    else

        B(i,1)=B(i,1)+M0(x,y-1);
    end
    
    if ~isempty(an)
        A(i,an)=1;
    else

        B(i,1)=B(i,1)+M0(x-1,y);
    end
    
    if ~isempty(ae)
        A(i,ae)=1;
    else
        B(i,1)=B(i,1)+M0(x,y+1);
    end
    
    if ~isempty(as)
        A(i,as)=1;
    else
        B(i,1)=B(i,1)+M0(x+1,y);
    end
    
    % intial conditions
    T0(i,1)=M0(x,y);
    
end
% Sparsity attributes to accelerate process
A=sparse(A);
B=sparse(B);

%% Computing solution
alpha=beta(1)/(beta2*dx^2);
Ti=zeros(length(T0),length(time));
Q=zeros(length(comp.coords),length(time));

% fixing Q term only on the part!!!
for i=1:length(comp.coords)
    if ismember(comp.coords(i,1:2),part.coords(:,1:2),'rows')
        Q(i,1)=beta(3);
    end
end

% initial temperature from data
Ti(:,1)=T0;

% next steps using Forward-Euler
Ap=alpha*ts*A;
Bp=alpha*ts*B;
Q=ts/beta2*Q;
for k=1:length(time)-1 % in time
    Ti(:,k+1)=Ap*Ti(:,k)+Bp+Ti(:,k)+Q(:,k);
end

disp('')
disp('coefficient,')
disp(alpha)
disp('')
%% Storing the full story
for i=1:length(time)
    for x=1:xmax
        for y=1:ymax
                S(x,y,i)=ic(1).mat(x,y);
        end
    end
end

for i=1:length(time)
    for j=1:length(comp.coords)
        x=comp.coords(j,1);
        y=comp.coords(j,2);
        S(x,y,i)=Ti(j,i);
        
    end
end
