function [BestSol,BestFitness,BestFitIter,Pop,obj] = TLBO_OPTIMIZATION(FITNESSFCN,lb,ub,Npop,T,ti,dt,tf,xp,yp,vrt,alphart,fs,kv,kalpha)

% Teaching Learning Based Optimization (TLBO) Algorithm
%
% This function implements the TLBO algorithm to solve optimization problems
% of the form: min F(X) subject to lb <= X <= ub, where X is a solution vector.
%
% Syntax:
% [BestSol, BestFitness, BestFitIter, Pop, obj] = TLBO_OPTIMIZATION(FITNESSFCN, lb, ub, Npop, T, ti, dt, tf, xp, yp, vrt, alphart, fs, kv, kalpha)
%
% Input Arguments:
% FITNESSFCN - Function handle to the fitness function to be minimized.
% lb         - Vector of lower bounds for the solution components.
% ub         - Vector of upper bounds for the solution components.
% Npop       - Size of the population (class size).
% T          - Total number of iterations (teaching periods).
% ti         - Initial time for simulation or analysis.
% dt         - Time step for simulation or analysis.
% tf         - Final time for simulation or analysis.
% xp, yp          - Coordinates representing the position (x, y) of the anemometer in the simulation space.
% vrt, alphart  - Recorded slowly-varying mean wind speed and direction
% fs                - Sampling frequency of the recorded data
% kv, kalpha - Weight coefficients for different components of the fitness function.
%
% Output Arguments:
% BestSol      - Best solution found by the algorithm.
% BestFitness- Fitness value of the best solution.
% BestFitIter- Best fitness value at each iteration.
% Pop        - Final population at the end of the algorithm execution.
% obj        - Objective function values for the final population.
%% Starting of STLBO

% Preallocation to store the best objective function value of every iteration
BestFitIter = NaN(T+1,1);              % Vector to store the fitness function value in every iteration (+1 initialization store)

% Preallocation to store the objective function value of every student
obj = NaN(Npop,1);                      % Vector to store the fitness function value of the population

% Determinig the number of decision variable in the problem (determing the
% size of the problem)
D = length(lb);

% Generation of initial population
Pop = repmat(lb,Npop,1) + repmat((ub-lb),Npop,1).*rand(Npop,D);

% Evaluation of objective function for each member of the population
% NOTE: CAN BE VECTORIZED
for p = 1: Npop
    X = Pop(p,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj(p) = FITNESSFCN(X,ti,dt,tf,xp,yp,vrt,alphart,fs,kv,kalpha);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Storing the best fitness for the zero iteration
BestFitIter(1) = min(obj);

%% Iteration loop

for gen = 1: T
    
    for i = 1: Npop
        %% ------Begining of the Teacher Phase for the ith student--------%
        
        % Determing the mean of the population
        Xmean = mean(Pop);  % Mean Student (NOTA: qui occorre usare le statistiche direzionali per la storm direction e abl direction).
        
        % Determing the location of the teacher
        [~,ind] = min(obj);
        % Copying the solution acting as teacher
        Xbest = Pop(ind,:); % Best Student
        
        % Determination of the teaching factor
        TF = randi([1,2],1,1); % Generating either 1 or 2 randomly for the teaching factor
        
        % Generation of a new solution
        Xnew = Pop(i,:) + rand(1,D).*(Xbest - TF*Xmean);

        % Boundig of the solution Upper Bound (UB)
        for l = 1:D
            if Xnew(l) > ub(l)
                Xnew(l) = lb(l) + rand*(ub(l)-lb(l));
            else
                Xnew(l) = Xnew(l);
            end
        end

        % Boundig of the solution Lower Bound (LB)
        for l = 1:D
            if Xnew(l) < lb(l)
                Xnew(l) = lb(l) + rand*(ub(l)-lb(l));
            else
                Xnew(l) = Xnew(l);
            end
        end
        
        
        
        % Xnew = min(ub,Xnew);    % Bounding the violating variables to their upper bound         
        % Xnew = max(lb,Xnew);    % Bounding the violating variables to their lower bound
        
        % Evaluation of objective function of the newly generated solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        NewObj = FITNESSFCN(Xnew,ti,dt,tf,xp,yp,vrt,alphart,fs,kv,kalpha);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Greedy selection
        if (NewObj < obj(i))
            Pop(i,:) = Xnew;  % Include the new solution in the population
            obj(i) = NewObj;  % Include the fitness function value of the new solution
        end
        %---------Ending of the Teacher Phase for the ith student---------%
        
        
        %% ------Begining of the Learner Phase for the ith student--------%
        
        % Selection of the random partner
        Partner = randi([1 Npop],1,1);
        
        % Ensuring that the current i-th student is not the partner
        while i == Partner
            Partner = randi([1 Npop],1,1);  % Selection of the random partner
        end
        
        % Generation of a new solution
        if (obj(i) < obj(Partner))
            Xnew = Pop(i,:) + rand(1,D).*(Pop(i,:)-Pop(Partner,:));
        else
            Xnew = Pop(i,:) - rand(1,D).*(Pop(i,:)-Pop(Partner,:));
        end
        
        % Boundig of the solution
        for l = 1:D
            if Xnew(l) > ub(l)
                Xnew(l) = lb(l) + rand*(ub(l)-lb(l));
            else
                Xnew(l) = Xnew(l);
            end
        end
        
        % Boundig of the solution LB
        for l = 1:D
            if Xnew(l) < lb(l)
                Xnew(l) = lb(l) + rand*(ub(l)-lb(l));
            else
                Xnew(l) = Xnew(l);
            end
        end
        % Xnew = min(ub,Xnew);    % Bounding the violating variables to their upper bound
        % Xnew = max(lb,Xnew);    % Bounding the violating variables to their lower bound
        
        % Evaluation of objective function of the newly generated solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        NewObj = FITNESSFCN(Xnew,ti,dt,tf,xp,yp,vrt,alphart,fs,kv,kalpha);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Greedy selection
        if (NewObj < obj(i))
            Pop(i,:) = Xnew;  % Include the new solution in the population
            obj(i) = NewObj;  % Include the fitness function value of the new solution
        end
        % --------Ending of the Learner Phase for the ith student---------%
    end
    
    % This is not part of the algorithm but is used to keep track of the
    % best solution determined till the current iteration
    BestFitIter(gen+1) = min(obj);   % Storing the best value of each iteration   

end

%% Extracting the best solution

[BestFitness,ind] = min(obj);
BestSol = Pop(ind,:);





