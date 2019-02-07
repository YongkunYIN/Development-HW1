
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Development HW2. Question 1 Part 2.
%%%% Yongkun Yin.
%%%% CEMFI.
%%%% 2019-02-01.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;

N = 1000; % number of individuals
Y = 40; % years

beta = 0.99^(1/12);
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

Sigma2_m = [0.043, 0.085, 0.171;...
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
0.068, 0.137, 0.273];

c = @(z,g,e2,e1)(z * g * e2 * e1);

rng(20190201);

U = normrnd(0, sigma_u, N, 1);
U = exp(U);
Z = exp(-sigma_u^2/2) * U;

E1 = normrnd(0, sigma_e1,N,Y*12);
E1 = exp(E1);
E1 = exp(-sigma_e1^2/2) * E1;

E2 = zeros(1000,480,3);
for j=1:1:3
    for i=1:1:N
        for y=1:1:Y
            for m=1:1:12
            E2(i,(y-1)*12+m,j) = normrnd(0, Sigma2_m(m,j)^(1/2));
            E2(i,(y-1)*12+m,j) = exp(E2(i,(y-1)*12+m,j));
            E2(i,(y-1)*12+m,j) = exp(-Sigma2_m(m,j)/2) * E2(i,(y-1)*12+m,j); 
            end
        end
    end
end


%% Part a.
% original life time utility.
C_o = zeros(N,Y*12,9);
W_o = zeros(N,9);
for j1=1:1:3
    for j2=1:1:3
        for i=1:1:N
            for y=1:1:Y
                for m=1:1:12
                    C_o(i,(y-1)*12+m,(j1-1)*3+j2) = c(Z(i,1),G(m,j1),E2(i,(y-1)*12+m,j2),E1(i,(y-1)*12+m));
                    W_o(i,(j1-1)*3+j2) = W_o(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * log(C_o(i,(y-1)*12+m,(j1-1)*3+j2));
                end
            end
        end
    end
end

% life time utility after removing the seasonal risk.
C_noS = zeros(N,Y*12,9);
W_noS = zeros(N,9);
for j1=1:1:3
    for j2=1:1:3
        for i=1:1:N
            for y=1:1:Y
                for m=1:1:12
                    C_noS(i,(y-1)*12+m,(j1-1)*3+j2) = c(Z(i,1),1,1,E1(i,(y-1)*12+m)); % replace G(m,j) with 1.
                    W_noS(i,(j1-1)*3+j2) = W_noS(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * log(C_noS(i,(y-1)*12+m,(j1-1)*3+j2));
                end
            end
        end 
    end
end


Gain_noS = W_noS-W_o;
Gain_noS = exp(Gain_noS*(1-beta)/(1-beta^(Y*12)))-1;

%% Part b.
% life time utility after removing the non-seasonal risk.
% original life time utility.
C_noNonS = zeros(N,Y*12,9);
W_noNonS = zeros(N,9);
for j1=1:1:3
    for j2=1:1:3
        for i=1:1:N
            for y=1:1:Y
                for m=1:1:12
                    C_noNonS(i,(y-1)*12+m,(j1-1)*3+j2) = c(Z(i,1),G(m,j1),E2(i,(y-1)*12+m,j2),1); % replace E1(i,(y-1)*12+m,1) with 1.
                    W_noNonS(i,(j1-1)*3+j2) = W_noNonS(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * log(C_noNonS(i,(y-1)*12+m,(j1-1)*3+j2));
                end
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
W_o_ita2 = zeros(N,9);
W_noS_ita2 = zeros(N,9);
W_noNonS_ita2 = zeros(N,9);

for j1=1:1:3
    for j2=1:1:3
        for i=1:1:N
            for y=1:1:Y
                for m=1:1:12
                    W_o_ita2(i,(j1-1)*3+j2) = W_o_ita2(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * (C_o(i,(y-1)*12+m,(j1-1)*3+j2))^(1-ita)/(1-ita);
                    W_noS_ita2(i,(j1-1)*3+j2) = W_noS_ita2(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * (C_noS(i,(y-1)*12+m,(j1-1)*3+j2))^(1-ita)/(1-ita);
                    W_noNonS_ita2(i,(j1-1)*3+j2) = W_noNonS_ita2(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * (C_noNonS(i,(y-1)*12+m,(j1-1)*3+j2))^(1-ita)/(1-ita);
                end
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
W_o_ita4 = zeros(N,9);
W_noS_ita4 = zeros(N,9);
W_noNonS_ita4 = zeros(N,9);

for j1=1:1:3
    for j2=1:1:3
        for i=1:1:N
            for y=1:1:Y
                for m=1:1:12
                    W_o_ita4(i,(j1-1)*3+j2) = W_o_ita4(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * (C_o(i,(y-1)*12+m,(j1-1)*3+j2))^(1-ita)/(1-ita);
                    W_noS_ita4(i,(j1-1)*3+j2) = W_noS_ita4(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * (C_noS(i,(y-1)*12+m,(j1-1)*3+j2))^(1-ita)/(1-ita);
                    W_noNonS_ita4(i,(j1-1)*3+j2) = W_noNonS_ita4(i,(j1-1)*3+j2) + beta^((y-1)*12+m-1) * (C_noNonS(i,(y-1)*12+m,(j1-1)*3+j2))^(1-ita)/(1-ita);
                end
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

