
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Development HW2. Question 2 Part a and b.
%%%% Yongkun Yin.
%%%% CEMFI.
%%%% 2019-02-02.
%%%% Please set Corr_S = 1 if the seasonal parts of consumption and labor 
%%%% are positively correlated; and set Corr_S = -1 if they are negatively
%%%% correlated.
%%%% Please set Corr_NS = 1 if the non-seasonal parts of consumption and labor 
%%%% are positively correlated; and set Corr_NS = -1 if they are negatively
%%%% correlated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;

% set Corr_S = 1 if consumption and labor is positively correlated;
% Set Corr_S = -1 if consumption and labor is negatively correlated.
Corr_S = -1;    % Corr_S means correlation between the seasonal parts.

% set Corr_NS = 1 if the non-seasonal parts of consumption and labor are
% positively correlated;
% Set Corr_NS = -1 if they are negatively correlated.
Corr_NS = -1;    % Corr_S means correlation between the non-seasonal parts.

N = 1000; % number of individuals
Y = 40; % years

beta = 0.99^(1/12);
kappa = -2*0.66/(28.5*30/7)^2;

sigma_e1 = 0.2^(1/2);
sigma_u = 0.2^(1/2);

G = [0.932, 0.863, 0.727;...
    0.845, 0.691, 0.381;...
    1.076, 1.151, 1.303;...
    1.070, 1.140, 1.280;...
    1.047, 1.094, 1.188;...
    1.030, 1.060, 1.119;...
    1.018, 1.037, 1.073;...
    1.018, 1.037, 1.073;...
    1.018, 1.037, 1.073;...
    1.001, 1.002, 1.004;...
    0.984, 0.968, 0.935;...
    0.961, 0.921, 0.843];

if Corr_S==1
    rho_S = 1;
    G_labor = G;
else
    rho_S = -1;
    G_labor = [1.070, 1.140, 1.280;...
    1.076, 1.151, 1.303;...
    0.845, 0.691, 0.381;...
    0.932, 0.863, 0.727;...
    0.961, 0.921, 0.843;...
    0.984, 0.968, 0.935;...
    1.001, 1.002, 1.004;...
    1.018, 1.037, 1.073;...
    1.018, 1.037, 1.073;...
    1.018, 1.037, 1.073;...
    1.030, 1.060, 1.119;...
    1.047, 1.094, 1.188];
end

if Corr_NS==1
    rho_NS= 1;
else
    rho_NS = -1;
end

Sigma_m = [0.043, 0.085, 0.171;...
0.034, 0.068, 0.137;...
0.145, 0.290, 0.580;...
0.142, 0.283, 0.567;...
0.137, 0.273, 0.546;...
0.137, 0.273, 0.546;...
0.119, 0.239, 0.478;...
0.102, 0.205, 0.410;...
0.094, 0.188, 0.376;...
0.094, 0.188, 0.376;...
0.085, 0.171, 0.341;...
0.068, 0.137, 0.273].^(1/2);

c = @(z,g,e2,e1)(z * g * e2 * e1);
h = @(z,g,e2,e1)(z * g * e2 * e1);

rng(20190202);

% generate random numbers for consumption.

U_consum = normrnd(0, sigma_u, N, 1);
U_consum = exp(U_consum);
Z_consum = exp(-sigma_u^2/2) * U_consum;

% generate random numbers for working hours.

U_labor = normrnd(0, sigma_u, N, 1);
U_labor = exp(U_labor);
Z_labor = exp(-sigma_u^2/2) * U_labor * 28.5*30/7;

% generate E2_consum and E2_labor that are jointly distributed.

E1_consum = zeros(N, Y*12);
E1_labor = zeros(N, Y*12);

for i=1:1:N
    for t=1:1:Y*12
        mu = [0,0];
        sigma = [sigma_e1^2, rho_NS*sigma_e1^2; rho_NS*sigma_e1^2, sigma_e1^2];
        r = mvnrnd(mu, sigma, 1);
        E1_consum(i,t) = r(1);
        E1_labor(i,t) = r(2);
    end
end

E1_consum = exp(E1_consum);
E1_consum = exp(-sigma_e1^2/2) * E1_consum;

E1_labor = exp(E1_labor);
E1_labor = exp(-sigma_e1^2/2) * E1_labor;

% generate E2_consum and E2_labor that are jointly distributed.

E2_consum = zeros(N,Y*12,3);
E2_labor = zeros(N,Y*12,3);
for j=1:1:3
    for i=1:1:N
        for y=1:1:Y
            for m=1:1:12
                mu = [0,0];
                sigma = [Sigma_m(m,j)^2, rho_S*Sigma_m(m,j)^2; rho_S*Sigma_m(m,j)^2, Sigma_m(m,j)^2];
                r = mvnrnd(mu, sigma, 1);
                E2_consum(i,(y-1)*12+m,j) = r(1);
                E2_labor(i,(y-1)*12+m,j) = r(2);
                E2_consum(i,(y-1)*12+m,j) = exp(E2_consum(i,(y-1)*12+m,j));
                E2_consum(i,(y-1)*12+m,j) = exp(-Sigma_m(m,j)^2/2) * E2_consum(i,(y-1)*12+m,j);
                E2_labor(i,(y-1)*12+m,j) = exp(E2_labor(i,(y-1)*12+m,j));
                E2_labor(i,(y-1)*12+m,j) = exp(-Sigma_m(m,j)^2/2) * E2_labor(i,(y-1)*12+m,j);     
            end
        end
    end
end


figure1 = figure('color',[1 1 1]);
subplot(2,2,1);
plot(E2_consum(1,:,1),E2_labor(1,:,1),'+');
axis([0 7 0 7]);
xlabel('{exp(-\sigma_m^2/2) \epsilon_m} for consumption');
ylabel('{exp(-\sigma_m^2/2) \epsilon_m} for labor');
title('Low Seasonality');

subplot(2,2,2);
plot(E2_consum(1,:,2),E2_labor(1,:,2),'+');
axis([0 7 0 7]);
xlabel('{exp(-\sigma_m^2/2) \epsilon_m} for consumption');
ylabel('{exp(-\sigma_m^2/2) \epsilon_m} for labor');
title('Middle Seasonality');

subplot(2,2,3);
plot(E2_consum(1,:,3),E2_labor(1,:,3),'+');
axis([0 7 0 7]);
xlabel('{exp(-\sigma_m^2/2) \epsilon_m} for consumption');
ylabel('{exp(-\sigma_m^2/2) \epsilon_m} for labor');
title('High Seasonality');

% initial lifetime utility
% nine cases:
% low deterministic and low stochastic
% low deterministic and middle stochastic
% low deterministic and high stochastic
% middle deterministic and low stochastic
% middle deterministic and middle stochastic
% middle deterministic and high stochastic
% high deterministic and low stochastic
% high deterministic and middle stochastic
% high deterministic and high stochastic


C_o = zeros(N,Y*12,9);
H_o = zeros(N,Y*12,9);
W_o = zeros(N,9);
for j1=1:1:3
    for j2 = 1:1:3
        for i=1:1:N
            for y=1:1:Y
                for m=1:1:12
                    C_o(i,(y-1)*12+m,(j1-1)*3+j2) = c(Z_consum(i,1),G(m,j1),E2_consum(i,(y-1)*12+m,j2),E1_consum(i,(y-1)*12+m));
                    H_o(i,(y-1)*12+m,(j1-1)*3+j2) = h(Z_labor(i,1),G_labor(m,j1),E2_labor(i,(y-1)*12+m,j2),E1_labor(i,(y-1)*12+m));
                    W_o(i,(j1-1)*3+j2) = W_o(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * (log(C_o(i,(y-1)*12+m,(j1-1)*3+j2))+kappa*0.5*H_o(i,(y-1)*12+m,(j1-1)*3+j2)^2);
                end
            end
        end
    end
end

% lifetime utility after removing seasonality in consumption
C_noSC = zeros(N,Y*12,9);
H_noSC = zeros(N,Y*12,9);
W_noSC = zeros(N,9);
for j1=1:1:3
    for j2=1:1:3
        for i=1:1:N
            for y=1:1:Y
                for m=1:1:12
                    C_noSC(i,(y-1)*12+m,(j1-1)*3+j2) = c(Z_consum(i,1),1,1,E1_consum(i,(y-1)*12+m));
                    H_noSC(i,(y-1)*12+m,(j1-1)*3+j2) = h(Z_labor(i,1),G_labor(m,j1),E2_labor(i,(y-1)*12+m,j2),E1_labor(i,(y-1)*12+m));
                    W_noSC(i,(j1-1)*3+j2) = W_noSC(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * (log(C_noSC(i,(y-1)*12+m,(j1-1)*3+j2))+kappa*0.5*H_noSC(i,(y-1)*12+m,(j1-1)*3+j2)^2);
                end
            end
        end
    end
end

% gain (from original to no seasonality in consumption)
Gain_o_noSC = W_noSC-W_o;
Gain_o_noSC = exp(Gain_o_noSC*(1-beta)/(1-beta^(Y*12)))-1;

% lifetime utility after removing seasonality in both consumption and labor
% supply.
C_noSCH = zeros(N,Y*12,9);
H_noSCH = zeros(N,Y*12,9);
W_noSCH = zeros(N,9);
for j1=1:1:3
    for j2=1:1:3
        for i=1:1:N
            for y=1:1:Y
                for m=1:1:12
                    C_noSCH(i,(y-1)*12+m,(j1-1)*3+j2) = c(Z_consum(i,1),1,1,E1_consum(i,(y-1)*12+m));
                    H_noSCH(i,(y-1)*12+m,(j1-1)*3+j2) = h(Z_labor(i,1),1,1,E1_labor(i,(y-1)*12+m));
                    W_noSCH(i,(j1-1)*3+j2) = W_noSCH(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * (log(C_noSCH(i,(y-1)*12+m,(j1-1)*3+j2))+kappa*0.5*H_noSCH(i,(y-1)*12+m,(j1-1)*3+j2)^2);
                end
            end
        end
    end
end

% gain (from no seasonality in consumption to no seasonality in both consumption and labor supply)
Gain_noSC_noSCH = W_noSCH-W_noSC;
Gain_noSC_noSCH = exp(Gain_noSC_noSCH*(1-beta)/(1-beta^(Y*12)))-1;

% gain (from original to no seasonality in both consumption and labor supply)
Gain_o_noSCH = W_noSCH-W_o;
Gain_o_noSCH = exp(Gain_o_noSCH*(1-beta)/(1-beta^(Y*12)))-1;

Results = zeros(9,3);
Results(:,1) = median(Gain_o_noSCH)';
Results(:,2) = median(Gain_o_noSC)';
Results(:,3) = median(Gain_noSC_noSCH)';

Results