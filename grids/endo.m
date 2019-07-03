% endogeneous grid method for consumption problem

% process for earnings
wage = 4.0;
rho = 0.95;
sig = 0.03;
ne = 10;
[prob,eps,z]=tauchen(ne,0.0,rho,sig);
z = exp(z - 0.5*sig^2);

% other parameters
r = 0.03;
R = 1+r;
sigma = 3.0;
beta = 1/R;

% grid for wealth (future wealth)
wmin = 0.0;
wmax = 50.0;
nw = 50;
gapw = (wmax - wmin)/(nw-1);
gridw = wmin:gapw:wmax;
gridw = gridw';

% set up value function
value = zeros(nw,ne);
cons = 0.5*ones(nw,ne);

for e=1:ne
    rhs = beta*R*((cons).^(-sigma))*prob(e,:)';
    optc = rhs.^(-1/sigma);
    optw = (gridw + optc - wage*z(e))./R;
    cap = find(gridw<=optw(1));
    optc(cap) = gridw(cap)*R + wage*z(e) - gridw(1);
    optv = optc.^(1-sigma)./(1-sigma) + beta*(value*prob(e,:)');
    value = interp1(optw,optv,gridw);
end 



