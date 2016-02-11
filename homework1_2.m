% ECON 5121 HOMEWORK 1
% LINGXIAO ZHAO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 3
clc
clear

% Parameters
global delta alpha A rho sigma mu b kappa beta; 
delta=0.0081; % Separation rate
alpha=0.72; % Elasticity of matching
A=0.158; % Matching efficiency
rho=0.9895; % Autocorrelation of weekly productivity
sigma=0.0034; % Standard Deviation for innovations
mu=0.05; % Bargaining weight for workers
b=0.9; % Unemployment utility 
       % 0.95 is lowered to solve the problem
kappa=0.6; % Posting cost is about 3 days of output
beta=0.999; % Weekly discount rate

% Discretize the AR(1) as a Markov chain
N=50;
p=(1+rho)/2;
q=(1+rho)/2;
PI=[p,1-p;
    1-q,q];
for i=3:N
    % Rouwenhorst Method
    PI=(p.*[PI,zeros(i-1,1);zeros(1,i-1),0]... 
        +(1-p).*[zeros(i-1,1),PI;0,zeros(1,i-1)]...
        +(1-q).*[zeros(1,i-1),0;PI,zeros(i-1,1)]...
        +q.*[0,zeros(1,i-1);zeros(i-1,1),PI]);
    % Normalize to be a stochastic matrix
    Weight=[ones(1,i);0.5*ones(i-2,i);ones(1,i)];
    PI=Weight.*PI;
end;
nu=sigma*sqrt(((N-1)/(1-rho^2)));
% The unconditional mean of the AR(1) process is 1
Z=linspace(1-nu,1+nu,N); 
Z=Z',

% Solve the non-liner equations
global PI Z;
theta0=ones(50,1);
theta=fsolve(@nonlinearsolver,theta0);
p_theta=A.*theta.^(1-alpha);
q_theta=A.*theta.^(-alpha);

% Simulation
% Find invariant distribution pi
[V,d]=eig(PI');
pi=V(:,1)./sum(V(:,1)); % Normalize eigenvector
% Simulation
T=10000;
rng(19920124); % Control random number generation. seed=19920124
x=rand(T,1);
sim=nan(T,1);
sim(1)=find((cumsum(pi)>=x(1)),1); % Initial state
for j=1:T-1
    sim(j+1)=find((cumsum(PI(sim(j),:))>=x(j+1)),1);
end
theta_T=theta(sim);
Z_T=Z(sim);

% Result

subplot(2,1,1)
plot(Z,theta)
xlabel('productivity shock z');
ylabel('theta');

subplot (2,1,2)
plot (Z_T)
xlabel ('Time')
ylabel ('Productivity')

