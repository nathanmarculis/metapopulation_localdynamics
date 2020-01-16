% Set up the computational space
% maximum population size
xmax = 10;
% plot window
ymin = 0;
ymax = 2;
% number of computational points
N = 2^10;
% create the discretization for the population density x
x = linspace(0,2*xmax,N);
%Number of generations
T = 1000;
% right
rightv = find(x>xmax);
rightp = min(rightv-1);

% Growth parameters
y = x;
Kg = zeros(N,N);
growth = zeros(N,N);
phat = zeros(1,N);
% Growth rate
R = 5;
% Carrying capacity
K = 1;
% Growth function
% Beverton-Holt
% g = R*y./(1+(R-1)/K*y);
xmin = 0;
% Ricker
% g = y.*exp(R*(1-y/K));
% xmin = 0;
% Sigmoid Beverton-Holt
delta = 7;
g = R*y.^(delta)./(1+(K^(delta-1)*R-1)/K^delta*y.^(delta));
xK = find(x<K);
xKpoint = max(xK+1);
% simpsons(g(1:xKpoint)-x(1:xKpoint),0,K);
% standard deviation in the growth
sigma_g = 0.1;

% Dispersal parameters
z = x;
% slope of dispersal function
alpha = 0.5;
gamma = 0.1;
lambda = 0.1;
% standard deviation in dispersal
sigma_d = 0.1;
% Function that determines the dispersal
s = alpha*z;
Kd = zeros(N,N);
dispersal = zeros(N,N);

% Initial condition
p = zeros( 1, N);
temp    = find(1<=x & x<=2);
p(temp) = ones(size(p(temp)));

% Growth phase
tic;
P = zeros(T+1,N);
P(1,1:N) = p;
subplot(2,1,1)
plot(x,P(1,1:N),'k','Linewidth',2);
axis([xmin-1 xmax+1 ymin max(P(1,1:N))+0.1]);
title('t=0');
% pause(1);
for j=1:T
    for i=1:N
        Kg(i,1:N) = exp(-(x(i)-g).^2/(2*sigma_g^2))/sqrt(2*pi*sigma_g^2);
        % Kg(i,1:N)=Kg(i,1:N)/simpsons(Kg(i,1:N),xmin,xmax);
        growth(i,1:N) = Kg(i,1:N).*p;
        phat(i) = simpsons(growth(i,1:N),0,2*xmax);
        if phat(i) < 0
            phat(i)=0;
        end
    end
    % Compute the expected dispersal at density z
    Ez = simpsons(s.*phat,0,2*xmax);
    % Dispersal phase
    d = max(y,Ez);
    fract = simpsons(phat,0,2*xmax);
    if fract > 1
        phat = phat/fract;
    end
    fract = simpsons(phat,0,2*xmax);
    for i=1:N
        % Gaussian dispersal
        Kd(i,1:N) = exp(-(x(i)-d).^2/(2*sigma_d^2))/sqrt(2*pi*sigma_d^2);
        dispersal(i,1:N) = Kd(i,1:N).*phat; 
        % Laplce recolonization
        p(i) = simpsons(dispersal(i,1:N),0,2*xmax)...
            + (1-fract)*simpsons(1/(2*lambda)*exp(-(abs(x(i)-y))/lambda).*phat,0,2*xmax);
    end
    total = simpsons(p,0,2*xmax);
    if total>1
        p = p/total;
    end
    P(j+1,1:N) = p;
    % Plot the probability for each generation
    subplot(2,1,2)
    plot(x,P(j+1,1:N),'k','Linewidth',2);
    axis([-0.5 xmax ymin max(P(j+1,1:N))+0.1]);
    title(['t = ',num2str(j)]);
    xlabel('population density ($x$)');
    ylabel('metapopulation density ($p_t(x)$)');
    % pause(1);
    drawnow;
    %simpsons(P(j+1,1:N),0,2*xmax)
    simpsons(p,0,2*xmax)
end
toc;
