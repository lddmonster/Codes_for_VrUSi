


%% Twistzz

% reconflag = reshape(cat(1,ones(160,128),zeros(160,128)),1,320*128);
% R = simData_all(:,reconflag == 1);
R = simData_all;

%TwIST handlers
% Linear operator handlers
hR = @(x) R*x;
hRt = @(x) R'*x; % x in matrix: pixels_numbers * N_elements

% observed data
y = rf_all(:,1:frames_num);

% regularization parameter 
tau = 0.02;

 
% TwIST parameters
lambda1 = 1e-4;  
initarray = zeros(285*128,frames_num);
% stopping theshold
tolA = 1e-2;
Phi = @(x) weighted_L1(x);

% proxOptions = struct;
% proxOptions.Nz = 285;
% proxOptions.Nx = 128;
% proxOptions.RelGapTol = 1e-8;
% proxOptions.AbsGapTol = 1e-12;
% proxOptions.MaxIter = 1000;
% proxOptions.CheckEvery = 5;
% proxOptions.ChunkPixels = 1024;
% proxOptions.FailOnNonconvergence = true;

% Current active implementation: elementwise complex soft thresholding.
Psi = @(x,tau)soft2(x,tau);
% Psi = @(v,gamma) psi_weight_L1(v,gamma,proxOptions);


% stop criterium:  the relative change in the objective function 
% falls below 'ToleranceA'
tic();
[x_twist,x_debias_twist,obj_twist,...
    times_twist,debias_start_twist,mse]= ...
         TwIST(y,hR,tau.* max(max(abs(hRt(y)))), ...
         'Lambda',lambda1,...
         'AT', hRt,...
         'Debias',0,...
         'PSI',Psi,...
         'PHI',Phi,...
         'Monotone',1,...
         'Initialization',initarray,...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'Verbose', 1,...
         'SPARSE',1,...
         'MaxiterA',0.50e2);

h_fig = figure(10);
set(h_fig,'Units','characters',...
        'Position',[30 10 150 45])
subplot(2,2,1)
semilogy(times_twist,obj_twist, 'b', 'LineWidth',2)
legend('TwIST')
st=sprintf('\\lambda_1 = %2.1e',tau);
title(st);
xlabel('CPU time (sec)')
ylabel('Obj. function')
grid

fprintf('twist over, used %d s\n',toc());



