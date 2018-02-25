function [Gamopt,gamVec,Kopt] = SLPV_SF_discrete(Gasym,Theta1,Theta2,ThetaT,ThetaPlusT,Theta2T,Theta2PlusT,REGID,F_theta,FthetaNum,SwLogic,mu,lambda,MultiLyapunov)
%% This code is used for state-feedback (switching) LPV control of discrete LPV plants, but can also be used for LTI plant

%% Input parameters :
% Gasym: generalized plant with symbol variables
% Theta is a cell arry: ith element is gridded points for the ith region,
% for instance Theta ={[-1 0.1],[-0.1 0 1]} represents two subsets of [-1
% 0.1] and [-0.1 1] respectively.
% Theta2 is gridded poitns for theta2, Theta2 = [] if only one GS parameter
% d_Theta: d_Theta(1,:) is the bound for d_theta1 (derivative of theta1), d_Theta(2,:) is the bound for d_theta2, .... 


% Output parameters are:
% Gam: is gama value
% xopt: is lmi vector
% Kopt: is Khat vector
% Xopt is Xhat
% Yopt is Xhat
% Fcn_theta: a function handle, Ftheta(theta) will give all the scalar functions for
% the matrix variables, for instance in X = X0+ f1(theta)X1+f2(theta)X2,
% Fcn_theta(theta) = [1 f1(theta) f2(theta)]
%d_Fcn_theta: a function handle for the derivative of the scalar functions
% FthetaNum: a vector, each element show the number of constant matrices and matrices as a function of GS parameters theta1, theta2, ... 
%            for instance, [1 1 1] means one constant matrix, one matrix as
%            a function of theta1, one matrix as a function of theta2
% SwLogic: 0-non-switching, 1- hysteresis switching, 2-average-dwell-time
% switching, not used at this moment

%%
if isempty(Theta2T) || sum(Theta2T{1}) == 0
    GSParaNum = 1;
else
    GSParaNum = 2;
end
% GSParaNum = size(Theta,1); % row number = number of GS parameters
if GSParaNum == 1
    theta = ThetaT{1}(1); %% GS parameter is a column vector
elseif GSParaNum == 2
    theta = [ThetaT{1}(1); Theta2T(1)];
end
Ga = AugPltEval(Gasym,theta);
n =  size(Ga.A(:,:,1),1);
nw = size(Ga.B1(:,:,1),2);
nu = size(Ga.B2(:,:,1),2);
nz = size(Ga.C1(:,:,1),1);
ny = size(Ga.C2(:,:,1),1);

% subset parameters
regnum = max(REGID(:,1));
regnum1 = max(REGID(:,2));
regnum2 = max(REGID(:,3));


% vu_set and vu_global are respectively used for the subset condition and the condition for global stability
vu_subset = 1;
vu_global = 1;
if SwLogic == 1
    mu = 1;     
elseif SwLogic == 2
    if MultiLyapunov
        vu_global = 1-lambda;
    else
        vu_subset = 1-lambda;
        vu_global = 1-lambda;
    end
end

%% Defining LMI variable
setlmis([]);
for nomeainng =1 % just for hiding the following code using a for iteration
    Gam = lmivar(1,[1 1]);
    %% Note the order of matrices: X = X0+ theta1*X1+ theta1^2*X2+ theta2*X3+ theta2^2*X4
    for regid = 1:regnum 
        gam(regid) = lmivar(1,[1 1]);
        for Id_Ftheta=1:sum(FthetaNum)% 
            X(Id_Ftheta,regid)= lmivar(1,[n 1]); 
            if MultiLyapunov == 1
                Xhat(Id_Ftheta,regid)= lmivar(1,[n 1]); 
            else
                Xhat(Id_Ftheta,regid) = X(Id_Ftheta,regid);
            end
            L(Id_Ftheta,regid) = lmivar(2,[nu n]);
        end    
        % use constant G        
        G(regid) = lmivar(2,[n n]);
    end

    lminum = 0;      
    %% LMIs for each subset
    for regid = 1:regnum           
       %% Get the grid points for admissible set   
       regid1 = REGID(regid,2);regid2 = REGID(regid,3);
       thetaT = ThetaT{regid1};theta2T = Theta2T{regid2}; 
       thetaPlusT = ThetaPlusT{regid1};theta2PlusT = Theta2PlusT{regid2}; 
       %% Note that positivity of Lyapunov function is implicitly enforced in the Lyapunov conditions    
        for Id_theta1 = 1:length(thetaT)
            theta1 = thetaT(Id_theta1);
            theta1Plus = thetaPlusT(Id_theta1);  
            for Id_theta2 = 1:length(theta2T)
                theta2 = theta2T(Id_theta2);  
                theta2Plus = theta2PlusT(Id_theta2);                
                
                theta = [theta1;theta2]; 
                thetaPlus = [theta1Plus ;theta2Plus];                

                Ga = AugPltEval(Gasym, theta);
                A = Ga.A;
                B1 = Ga.B1; B2 = Ga.B2;
                C1 = Ga.C1; C2 = Ga.C2;
                D1 = Ga.D11; D2 = Ga.D12;
                
                Ftheta = F_theta(theta); 
                FthetaPlus = F_theta(thetaPlus);
           
                lminum = lminum+1;      
                % subset conditions
                lmiterm([-lminum 1 1 G(regid)],vu_subset,1, 's');
                for Id_Ftheta=1:sum(FthetaNum)   
                    lmiterm([-lminum 1 1 X(Id_Ftheta,regid)],Ftheta(Id_Ftheta)*vu_subset,-1); 
                    lmiterm([-lminum 2 1 L(Id_Ftheta,regid)],B2,Ftheta(Id_Ftheta));
                    lmiterm([-lminum 2 2 X(Id_Ftheta,regid)],FthetaPlus(Id_Ftheta),1);
                    lmiterm([-lminum 3 1 L(Id_Ftheta,regid)],D2,Ftheta(Id_Ftheta));
                end

                lmiterm([-lminum 2 1 G(regid)],A,1);
                lmiterm([-lminum 3 1 G(regid)],C1,1);
                lmiterm([-lminum 3 3 gam(regid)],1,1);


                lmiterm([-lminum 4 2 0],B1');
                lmiterm([-lminum 4 3 0],D1');

                lmiterm([-lminum 4 4 gam(regid)],1,1);

                if MultiLyapunov == 1 && SwLogic ~=0
                    % use Phat for guaranteeing global stability, which is different from P  
                    lminum = lminum+1;   
                    lmiterm([-lminum 1 1 G(regid)],vu_global,1, 's');
                    for Id_Ftheta=1:sum(FthetaNum)   
                        lmiterm([-lminum 1 1 Xhat(Id_Ftheta,regid)],vu_global*Ftheta(Id_Ftheta),-1);  

                        lmiterm([-lminum 2 1 L(Id_Ftheta,regid)],B2,Ftheta(Id_Ftheta));

                        lmiterm([-lminum 2 2 Xhat(Id_Ftheta,regid)],FthetaPlus(Id_Ftheta),1);
                    end
                    lmiterm([-lminum 2 1 G(regid)],A,1);                            
                end
            end % Id_theta2
        end % Id_theta1                
        lminum = lminum + 1;
        lmiterm([lminum 1 1 gam(regid)],1,1);
        lmiterm([lminum 1 1 Gam],-1,1);
    end %regid  
    
    %% LMIs for the monotonic property of Lyapunov function   
    lminum_sub = lminum; 
    %         disp(['The   number of lmis for subsets is ', num2str(lminum_sub)]);
    % switching surface index increases in the theta1 direction first
    % in the way of 
    % 4 5 6
    % 1 2 3
    % and then in theta2 direction in the way of 
    % 2 4
    % 1 3   
    if SwLogic ~= 0        
    %% for switching surfaces in theta1 direction     
    %         Note that to express X == Y in LMI lab, -[varepsi*E=Id_Ftheta X-Y; X-Y Id_Ftheta]<=0 is used
    for regid2 = 1:regnum2
        RegId1 = REGID(find(REGID(:,3)==regid2),2);
        RegId1 = RegId1'; RegId1 = sort(RegId1);     
        RegId = (regid2-1)*regnum1+RegId1;
        for id1 = RegId1(1:end-1)  
            for theta2 = unique(Theta2(regid2,:))%becasue inequalities is affine w.r.t theta1 & theta2, so check of vertices is enough
                    theta1 = Theta1(id1,end);                                                      
                    Ftheta = F_theta([theta1;theta2]);   
                    % mu*Pj > = Pk  Xj< mu*Xk
                    jk = [RegId(id1) RegId(id1+1)];                    
                    lminum = lminum + 1;
                    for Id_Ftheta = 1:length(Ftheta)
                        lmiterm([lminum 1 1 Xhat(Id_Ftheta,jk(1))],Ftheta(Id_Ftheta), 1);
                        lmiterm([lminum 1 1 Xhat(Id_Ftheta,jk(2))],-Ftheta(Id_Ftheta)*mu, 1);
                    end    
                    
                    
                    % the other direction
                    theta1 = Theta1(id1+1,1);
                    Ftheta = F_theta([theta1;theta2]);                                   
                    jk = [RegId(id1+1) RegId(id1)];     
                    lminum = lminum + 1;
                    for Id_Ftheta = 1:length(Ftheta)
                        lmiterm([lminum 1 1 Xhat(Id_Ftheta,jk(1))],Ftheta(Id_Ftheta), 1);
                        lmiterm([lminum 1 1 Xhat(Id_Ftheta,jk(2))],-Ftheta(Id_Ftheta)*mu, 1);
                    end  
            end   
        end %id1
    end % regid2       

    %% for switching surfaces in theta2 direction 
    for regid1 = 1:regnum1
        RegId2 = REGID(find(REGID(:,2)==regid1),3);
        RegId2 = RegId2';RegId2 = sort(RegId2);       
        RegId = (RegId2-1)*regnum1+regid1;
        for id2 = RegId2(1:end-1)
            for theta1 = Theta1(regid1,:) %becasue inequalities is affine w.r.t theta1 & theta2, so check of vertices are enough
                theta2 = Theta2(id2,end);
                Ftheta = F_theta([theta1;theta2]);   
                jk = [RegId(id2) RegId(id2+1)];
                lminum = lminum + 1;
                for Id_Ftheta = 1:length(Ftheta)
                    lmiterm([lminum 1 1 Xhat(Id_Ftheta,jk(1))],Ftheta(Id_Ftheta), 1);
                    lmiterm([lminum 1 1 Xhat(Id_Ftheta,jk(2))],-Ftheta(Id_Ftheta)*mu, 1);
                end 

                % the other direction
                theta2 = Theta2(id2+1,1);
                Ftheta = F_theta([theta1;theta2]);  
                jk = [RegId(id2+1) RegId(id2)];  
                lminum = lminum + 1;
                for Id_Ftheta = 1:length(Ftheta)
                    lmiterm([lminum 1 1 Xhat(Id_Ftheta,jk(1))],Ftheta(Id_Ftheta), 1);
                    lmiterm([lminum 1 1 Xhat(Id_Ftheta,jk(2))],-Ftheta(Id_Ftheta)*mu, 1);
                end
            end    
        end %id2
    end       
    end %if 
    %% Get and solve lmis
    lmisys = getlmis;
    lminum
    nvar = decnbr(lmisys); %number of decision variable.
    c = zeros(nvar,1);
    c(1)=1; 
    options(1)= 1e-4; % relative accuary on the optimal value
    options(2)= 1000; %Number of Iteration
    options(4) =  10; % J, the code terminates when the objective has not decreased by more than the desired relative accuracy during the last J iterations
    options(5)= 0; % 1 for not showing the process
    % initial value 
      
    %slove LMIs
    nvar = decnbr(lmisys); %number of decision variable.
    c = zeros(nvar,1);
    c(1)=1; 

    options(1)= 1e-5;
    options(2)= 500; %Number of Iteration
    options(5)= 0; % 1 not show the process
                    %% initial value
                    % if stpcont~=1
                    %     if XY==1 %Y0 is constatnt
                    %         xinit = mat2dec(lmisys,gamma0,0(:,:,1),Khat0.A(:,:,1),Khat0.B(:,:,1),Khat0.C(:,:,1),Khat0.D(:,:,1),X00(:,:,2),Khat0.A(:,:,2),Khat0.B(:,:,2),Khat0.C(:,:,2),Khat0.D(:,:,2),X00(:,:,3),Khat0.A(:,:,3),Khat0.B(:,:,3),Khat0.C(:,:,3),Khat0.D(:,:,3),Y00);
                    %     else
                    %         xinit = mat2dec(lmisys,gamma0,Y00(:,:,1),Khat0.A(:,:,1),Khat0.B(:,:,1),Khat0.C(:,:,1),Khat0.D(:,:,1),Y00(:,:,2),Khat0.A(:,:,2),Khat0.B(:,:,2),Khat0.C(:,:,2),Khat0.D(:,:,2),Y00(:,:,3),Khat0.A(:,:,3),Khat0.B(:,:,3),Khat0.C(:,:,3),Khat0.D(:,:,3),X00);
                    %     end
                    %     [gammaopt,xopt] = mincx(lmisys,c,options,xinit);
                    % else
                    % end
    [copt,xopt] = mincx(lmisys,c,options);
    % [gammaopt,xopt] = mincx(lmisys,c);
    % [gammaopt,xopt] = mincx(lmisys,c,options,[Khatopt0,Yopt0,Xopt0]);

    %% N and M  
    gamVec = ones(1,regnum)*Inf;
    if ~isempty(xopt)
        Gamopt = dec2mat(lmisys,xopt,Gam); 
        
        for regid = 1:regnum
            for Id_Ftheta = 1:sum(FthetaNum)
                Kopt.L(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,L(Id_Ftheta,regid));            
                Kopt.X(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,X(Id_Ftheta,regid));  
                Kopt.Xhat(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Xhat(Id_Ftheta,regid));  
            end
            Kopt.G(:,:,regid) = dec2mat(lmisys,xopt,G(regid));
            gamVec(regid) = dec2mat(lmisys,xopt,gam(regid));
        end   
    else
        Gamopt = inf;
        Kopt = [];
    end
end           

