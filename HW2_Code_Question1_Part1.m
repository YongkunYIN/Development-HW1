
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Development HW2. Question 1 Part 1.
%%%% Yongkun Yin.
%%%% CEMFI.
%%%% 2019-02-01.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;

N = 1000;                                  % number of individuals
Y = 40;                                    % years

beta = 0.99^(1/12);                        % discount factor
sigma_e1 = 0.2^(1/2);                      % standard deviation of epsilon
sigma_u = 0.2^(1/2);                       % s.d. of u
G = [0.932, 0.863, 0.727;...               % exp(g(m))
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

c = @(z,g,e1)(z * g * e1);            % consumption

rng(20190201);

U = normrnd(0, sigma_u, N, 1);             % draw u randomly
U = exp(U);
Z = exp(-sigma_u^2/2) * U;

E1 = normrnd(0, sigma_e1,N,Y*12);
E1 = exp(E1);
E1 = exp(-sigma_e1^2/2) * E1;

%% Part a.
% original life time utility. 
C_o = zeros(N,Y*12,3);    % Consumption. N agents. Each lives Y years. Three possible levels of uncertainties.
W_o = zeros(N,3);         % life time utility. 
for j=1:1:3
    for i=1:1:N
        for y=1:1:Y
            for m=1:1:12
            C_o(i,(y-1)*12+m,j) = c(Z(i,1),G(m,j),E1(i,(y-1)*12+m));
            W_o(i,j) = W_o(i,j) + beta^((y-1)*12+m-1) * log(C_o(i,(y-1)*12+m,j));
            end
        end
    end
end

% life time utility after removing the seasonal risk.
C_noS = zeros(N,Y*12,3);    % Consumption.
W_noS = zeros(N,3);         % life time utility.
for j=1:1:3
    for i=1:1:N
        for y=1:1:Y
            for m=1:1:12
            C_noS(i,(y-1)*12+m,j) = c(Z(i,1),1,E1(i,(y-1)*12+m)); % replace G(m,j) with 0.
            W_noS(i,j) = W_noS(i,j) + beta^((y-1)*12+m-1) * log(C_noS(i,(y-1)*12+m,j));
            end
        end 
    end
end

Gain_noS = W_noS-W_o;
Gain_noS = exp(Gain_noS*(1-beta)/(1-beta^(Y*12)))-1;

%% Part b.
% life time utility after removing the non-seasonal risk.
% original life time utility.
C_noNonS = zeros(N,Y*12,3);    % Consumption.
W_noNonS = zeros(N,3);         % life time utility.
for j=1:1:3
    for i=1:1:N
        for y=1:1:Y
            for m=1:1:12
            C_noNonS(i,(y-1)*12+m,j) = c(Z(i,1),G(m,j),1); % replace E1(i,(y-1)*12+m,1) with 1.
            W_noNonS(i,j) = W_noNonS(i,j) + beta^((y-1)*12+m-1) * log(C_noNonS(i,(y-1)*12+m,j));
            end
        end 
    end
end

Gain_noNonS = W_noNonS-W_o;
Gain_noNonS = exp(Gain_noNonS*(1-beta)/(1-beta^(Y*12)))-1;

%% Part d
%% ita = 2
ita = 2;
% life time utility original/ after removing seasonal risk/ after removing non-seasonal risk.
W_o_ita2 = zeros(N,3);
W_noS_ita2 = zeros(N,3);
W_noNonS_ita2 = zeros(N,3);

for j=1:1:3
    for i=1:1:N
        for y=1:1:Y
            for m=1:1:12
            W_o_ita2(i,j) = W_o_ita2(i,j) + beta^((y-1)*12+m-1) * (C_o(i,(y-1)*12+m,j))^(1-ita)/(1-ita);
            W_noS_ita2(i,j) = W_noS_ita2(i,j) + beta^((y-1)*12+m-1) * (C_noS(i,(y-1)*12+m,j))^(1-ita)/(1-ita);
            W_noNonS_ita2(i,j) = W_noNonS_ita2(i,j) + beta^((y-1)*12+m-1) * (C_noNonS(i,(y-1)*12+m,j))^(1-ita)/(1-ita);
            end
        end 
    end
end

Gain_noS_ita2 = W_noS_ita2./W_o_ita2;
Gain_noS_ita2 = Gain_noS_ita2.^(1/(1-ita))-1;

Gain_noNonS_ita2 = W_noNonS_ita2./W_o_ita2;
Gain_noNonS_ita2 = Gain_noNonS_ita2.^(1/(1-ita))-1;

%% ita = 4
ita = 4;
% life time utility original/ after removing seasonal risk/ after removing non-seasonal risk.
W_o_ita4 = zeros(N,3);
W_noS_ita4 = zeros(N,3);
W_noNonS_ita4 = zeros(N,3);

for j=1:1:3
    for i=1:1:N
        for y=1:1:Y
            for m=1:1:12
            W_o_ita4(i,j) = W_o_ita4(i,j) + beta^((y-1)*12+m-1) * (C_o(i,(y-1)*12+m,j))^(1-ita)/(1-ita);
            W_noS_ita4(i,j) = W_noS_ita4(i,j) + beta^((y-1)*12+m-1) * (C_noS(i,(y-1)*12+m,j))^(1-ita)/(1-ita);
            W_noNonS_ita4(i,j) = W_noNonS_ita4(i,j) + beta^((y-1)*12+m-1) * (C_noNonS(i,(y-1)*12+m,j))^(1-ita)/(1-ita);
            end
        end 
    end
end

Gain_noS_ita4 = W_noS_ita4./W_o_ita4;
Gain_noS_ita4 = Gain_noS_ita4.^(1/(1-ita))-1;

Gain_noNonS_ita4 = W_noNonS_ita4./W_o_ita4;
Gain_noNonS_ita4 = Gain_noNonS_ita4.^(1/(1-ita))-1;

display ("welfare gain of removing the seasonal component (ita = 1):");
median(Gain_noS)
display ("welfare gain of removing the non-seasonal component (ita = 1):");
median(Gain_noNonS)
display ("welfare gain of removing the seasonal component (ita = 2):");
median(Gain_noS_ita2)
display ("welfare gain of removing the non-seasonal component (ita = 2):");
median(Gain_noNonS_ita2)
display ("welfare gain of removing the seasonal component (ita = 4):");
median(Gain_noS_ita4)
display ("welfare gain of removing the non-seasonal component (ita = 4):");
median(Gain_noNonS_ita4)


% histograms of the gains.
figure(1);
subplot(3,2,1);
h1 = histogram(Gain_noS(:,1));
hold on;
h2 = histogram(Gain_noS(:,2));
hold on;
h3 = histogram(Gain_noS(:,3));
h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';
h1.NumBins = 20;
h2.NumBins = 20;
h3.NumBins = 20;
legend('low seasonality','middle seasonality','high seasonality');
title('\eta=1, removing seasonal risk');

subplot(3,2,2);
h4 = histogram(Gain_noNonS(:,1));
hold on;
h5 = histogram(Gain_noNonS(:,2));
hold on;
h6 = histogram(Gain_noNonS(:,3));
h4.Normalization = 'probability';
h5.Normalization = 'probability';
h6.Normalization = 'probability';
h4.NumBins = 20;
h5.NumBins = 20;
h6.NumBins = 20;
legend('low seasonality','middle seasonality','high seasonality');
title('\eta=1,removing nonseasonal risk');

subplot(3,2,3);
h1 = histogram(Gain_noS_ita2(:,1));
hold on;
h2 = histogram(Gain_noS_ita2(:,2));
hold on;
h3 = histogram(Gain_noS_ita2(:,3));
h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';
h1.NumBins = 20;
h2.NumBins = 20;
h3.NumBins = 20;
legend('low seasonality','middle seasonality','high seasonality');
title('\eta=2, removing seasonal risk');

subplot(3,2,4);
h4 = histogram(Gain_noNonS_ita2(:,1));
hold on;
h5 = histogram(Gain_noNonS_ita2(:,2));
hold on;
h6 = histogram(Gain_noNonS_ita2(:,3));
h4.Normalization = 'probability';
h5.Normalization = 'probability';
h6.Normalization = 'probability';
h4.NumBins = 20;
h5.NumBins = 20;
h6.NumBins = 20;
legend('low seasonality','middle seasonality','high seasonality');
title('\eta=2,removing nonseasonal risk');

subplot(3,2,5);
h1 = histogram(Gain_noS_ita4(:,1));
hold on;
h2 = histogram(Gain_noS_ita4(:,2));
hold on;
h3 = histogram(Gain_noS_ita4(:,3));
h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';
h1.NumBins = 20;
h2.NumBins = 20;
h3.NumBins = 20;
legend('low seasonality','middle seasonality','high seasonality');
title('\eta=4, removing seasonal risk');

subplot(3,2,6);
h4 = histogram(Gain_noNonS_ita4(:,1));
hold on;
h5 = histogram(Gain_noNonS_ita4(:,2));
hold on;
h6 = histogram(Gain_noNonS_ita4(:,3));
h4.Normalization = 'probability';
h5.Normalization = 'probability';
h6.Normalization = 'probability';
h4.NumBins = 20;
h5.NumBins = 20;
h6.NumBins = 20;
legend('low seasonality','middle seasonality','high seasonality');
title('\eta=4,removing nonseasonal risk');
