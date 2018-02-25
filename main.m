%*************************************************************************
%  
%  Discrete-time state-feedback switching LPV control with separate Lyapunov functions 
%                        for stability and local performance
%             ----------------------------------------
%
% By Pan Zhao, Control Engineering Lab, UBC, BC, Canada
% Under supervision of Prof. Ryozo Nagamune.
% Creation: Sep 11, 2017.

% Function: Main function
% Reference: P. Zhao and R. Nagamune, \"Discrete-Time State-Feedback Switching LPV Control with Separate Lyapunov
% Functions for Stability and Local Performance"," Accepted by American Control Conference, 2018.

 
clear;
clc;
%% Settings
SeparateLyapunov = 1;               % 1 (0) for using separate (common) Lyapunov functions for
SwLogic = 1;                        % 0: non-switching, 1: arbitrary switching, 2:ADT switching 
BiSearch = 1;                       % do a bisection-search for determining mu value in ADT switching
MaxNum_bs = 10;                     % maximum number of trials for bisection-search 
muSearchRange = [10^0.001 10^2];    % search range for mu in average-dwell-time switching

% for ADT switching
tau = 1000;
mu = 2;
lambda = 1-mu^(-1/tau);
% tau = - log(mu)/log(1-lambda);

% no need for bisearch if not for ADT switching
if SwLogic ~=2
    BiSearch = 0;
end

% test different partitioning
% Theta1 =[-0.7 1];
% Theta1 =[-0.7 0.2;0.1  1];
Theta1 =[-0.7 0;-0.1 0.5;0.4 1];
% Theta1 =[-0.7 0;-0.1  1];

% Theta1 = [-1 1];
% Theta1 = [-1 0.1;0 1];
% Theta1 = [-1 -0.2;-0.3 0.5;0.4 1];

Theta2 =0;
delta = 1;      % maximum parameter deviation at one sampling interval  
thetaMin = Theta1(1,1);
thetaMax = Theta1(end,end);

if SwLogic == 0
    Theta1 = [Theta1(1,1) Theta1(end,end)];
end

%% Plant matrices
syms 
syms theta1 theta2;  
Gasym.A = 0.5*[1-theta1 0 -2+theta1;2-theta1 -1 1-theta1; -1+theta1 1-3*theta1 -theta1] ;
Gasym.B1 = [0;1-theta1;theta1];
Gasym.B2 = [1 0 0]';
Gasym.C1 = [1 1 1];
Gasym.D11 = 0; Gasym.D12 = 0;
Gasym.C2 = eye(3); Gasym.D21 = 0; Gasym.D22 = 0;


F_theta  = @(x) [1 x(1)]; %function for PD matrices
FthetaNum = [1 1]; %affine form

% get the region number information
regnum1 = size(Theta1,1); %num. of subsets for theta1, only partition theta_1
regnum2 = size(Theta2,1); % num. of subsets for theta2;
regnum = regnum1 * regnum2; % Total num. of subsets. S shape from the bottom to order them, an example:
% Determination of regid1 & regid2  
REGID = zeros(regnum,3);
for regid = 1:regnum
   if mod(regid,regnum1) == 0  
        regid1 = regnum1;
    else
        regid1 = mod(regid,regnum1);
   end
   if regid > regnum1
        regid2 = 1+ floor(regid/(regnum1+0.1));
    else
        regid2 = 1;
   end  
   REGID(regid,:) = [regid regid1 regid2];   
end

% get the admissible region of (theta, theta_plus)
for regid1 = 1:regnum1
    [ThetaT{regid1},ThetaPlusT{regid1}] = AdmRegDiscrete(Theta1(regid1,:),thetaMin, thetaMax,delta); 
end
% delta_l = ; delta_u = 0;
if sum(Theta2(1,:)) ~= 0
    for regid2 = 1:regnum2
        [Theta2T{regid2},Theta2PlusT{regid2}] = AdmRegDiscrete(Theta2(regid2,:),thetaMin, thetaMax,delta); 
    end
else       
    Theta2T = {0}; Theta2PlusT = {0};
end

if BiSearch == 0    
    [Gamopt,gamVec,Kopt] = SLPV_SF_discrete(Gasym,Theta1,Theta2,ThetaT,ThetaPlusT, Theta2T, Theta2PlusT,REGID,F_theta,FthetaNum,SwLogic,mu,lambda,SeparateLyapunov);
else    
    Result_bs = ones(MaxNum_bs+2,2)*inf; % 1st column for mu, second column for Gam
    Result_bs([1 2],1) = muSearchRange';
   
    
    for i = 1:MaxNum_bs+2
        if i <=2
            mu = Result_bs(i,1);           
        else
            % sort the mu value in ascending order 
            [~,index] = sort(Result_bs(:,1));
            Result_bs = Result_bs(index,:);
            % sort the gamma value in ascending order 
            [~,index] = sort(Result_bs(:,2));
            % find the mu value corresponding to the smallest two gamma values and take their average as the mu value for next trial
            mu = sum(Result_bs(index(1:2),1))/2;
            Result_bs(i,1) = mu;
        end
        lambda = 1-mu^(-1/tau);
        [Gamopt,gamVec,Kopt] = SLPV_SF_discrete(Gasym,Theta1,Theta2,ThetaT,ThetaPlusT, Theta2T, Theta2PlusT,REGID,F_theta,FthetaNum,SwLogic,mu,lambda,SeparateLyapunov);
        Result_bs(i,2) = Gamopt;
    end
    Result_bs
    [Gamopt, index] = min(Result_bs(:,2));    
    fprintf('The smallest gamma is %.3f obtained at mu= %.3e', Gamopt, Result_bs(index,1));    
end


for regid = 1: regnum
    Ginv(:,:,regid) = inv(Kopt.G(:,:,regid));
end

L = Kopt.L;