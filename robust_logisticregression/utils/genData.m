% Author: Jakramate Bootkrajang
% Desc: Generate synthetic data generatively/discriminatively
% label from this function is 1...N

function [X Y FD Xt Yt FDt] = genData(CLS, DIM, DPNT, TPNT, CSEP, OPT)

switch OPT
    
    case 'gen'  % generate data in generative way
        
        % mean of distributions
        mu = zeros(CLS, DIM);
        for i = 1:CLS
            % each class sits on each axis
            % pair-wise csep
            %mu(i,i) = mu(i,i) + CSEP*(sqrt(DIM/2));
            
            % classes lie on a line
            % csep between adjacent class is CSEP
            mu(i,:) = mu(i,:) + (i-1)*CSEP;
        end
        
        % covar of distribution, which is now identity matrix
        covar = eye(DIM);%.*randi([1 100],DIM,DIM);

        dMix  = gmm(DIM, CLS, 'full');

        % setup data mixture
        for i = 1:CLS
            dMix.centres(i,:)    = mu(i,:);            
            dMix.covars(:, :, i) = covar ;%+ ((i-1)*4);
        end

        % generating data
        % ensuring there are 2 classes, for small sample size
       % [X  Y] = gmmsamp(dMix, DPNT);
       % while (min(sum(Y==1),sum(Y==2))/DPNT < 0.3)
            [X  Y] = gmmsamp(dMix, DPNT);
       % end
        [Xt Yt] = gmmsamp(dMix, TPNT);        
                
    case 'dis' % favouring discriminative classifier
       
        if(CLS ~= 2)
            error('Expecting CLS=2 for discriminative data generation');
        end
 
        w = ones(DIM,1);
        % 1 relevant features
%        w = [10 zeros(1,DIM-1)]';
        
        % 3 relevant features
       % w = [10/sqrt(3) 10/sqrt(3) 10/sqrt(3) zeros(1,DIM-3)]';       
       
        % Exponentially many relevant features
%          w = zeros(1,DIM)';
%          for i = 1:DIM
%              w(i) = 0.5^(i-1)*sqrt(75);
%          end
        
        dMix                 = gmm(DIM, 1, 'full');        
        dMix.centres(1,:)    = zeros(1,DIM); 
        dMix.covars(:, :, 1) = eye(DIM) * CSEP;
        
        [X  y ]  = gmmsamp(dMix, DPNT);
        [Xt yt]  = gmmsamp(dMix, TPNT);                        
        
        % pass thru sigmoid
        xLogit  = sigmoid(w' *  X');
        xtLogit = sigmoid(w' * Xt');
        
        Y  = xLogit  > 0.5;
        Yt = xtLogit > 0.5;
        
        Y  = Y';
        Yt = Yt';  
                
        Y  = Y  + 1;
        Yt = Yt + 1;
        
        while (min(sum(Y==1),sum(Y==2))/DPNT < 0.3)
            dMix                 = gmm(DIM, 1, 'full');
            dMix.centres(1,:)    = zeros(1,DIM);
            dMix.covars(:, :, 1) = eye(DIM) * CSEP;
            
            [X  y ]  = gmmsamp(dMix, DPNT);
            [Xt yt]  = gmmsamp(dMix, TPNT);
            
            % pass thru sigmoid
            xLogit  = sigmoid(w' *  X');
            xtLogit = sigmoid(w' * Xt');
            
            Y  = xLogit  > 0.5;
            Yt = xtLogit > 0.5;
            
            Y  = Y';
            Yt = Yt';
            
            Y  = Y  + 1;
            Yt = Yt + 1;
        end
        
    case 'uni'
        % generating uniform data
        X   = zeros(DPNT, DIM);
        Xt  = zeros(TPNT, DIM);
        
        LIM = DPNT/2;
        
        X(1:LIM,      :) = unifrnd(-1,3, [LIM DIM]);
        Y(1:LIM,      :) = 0;
        X(LIM+1:DPNT, :) = unifrnd(-3,1, [LIM DIM]);
        Y(LIM+1:DPNT, :) = 1;
        
        LIM = TPNT/2;
        
        Xt(1:LIM,     :) = unifrnd(-1,3, [LIM DIM]);
        Yt(1:LIM,     :) = 0;
        Xt(LIM+1:TPNT,:) = unifrnd(-3,1, [LIM DIM]);
        Yt(LIM+1:TPNT,:) = 1;        

        Y  = Y  + 1;
        Yt = Yt + 1;                
end

FD  = zeros(DPNT,1); % flip descriptor for training data
FDt = zeros(TPNT,1); % flip descriptor for testing  data