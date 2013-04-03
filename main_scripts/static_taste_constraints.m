function [H, f, A, b, Aeq, beq, lb, ub] = static_taste_constraints (pin, qin, good, smooth, norm)
%--------------------------------------------------------------------------
% PURPOSE: linear programming constraints and quadratic objective
%          function for the taste adjusted model.
%--------------------------------------------------------------------------

    %----------------------------------------------------------------------
    % BEGIN MAIN PROGRAM
    %----------------------------------------------------------------------

    % U[s] - U[t] + mu[s](q[1,s] - q[1,t]) - lambda[t] * p[t]' * (q[s] - q[t])) <= 0

    %----------------------------------------------------------------------
    % Build A
    T = size(qin, 2);
    for i = 1 : 1 : T;
        % Constraints pertaining to U's
        tmp1 = eye (T);
        tmp1(:, i) = -ones (T, 1);
        tmp1(i, :) = zeros(1, T);

        % Constraints pertaining to lambda's
        tmp2 = zeros (T, T);
        for j = 1 : 1 : T;
            tmp2(j, i) = -pin(:, i)' * (qin(:, j) - qin(:, i));
        end;

        % Constraints pertaining to mu's
        tmp3 = zeros (T, T);
        for j = 1 : 1 : T;
            tmp3(j, i) = (qin(good, j) - qin(good, i));
        end;
        
        % Putting together
        tmp4 = [tmp1 tmp2 tmp3];

        if i == 1;
            A = tmp4;
        else
            A = [A; tmp4]; %#ok<AGROW>
        end;
    end;
    
     
    % Positive shadow price
    Apos = zeros(T, 3*T);
    for t = 1 : T;
        Apos(t, 2*T + t) = 1;
        Apos(t, T + t) = -pin(good, t);
    end;
    A = [A; Apos];
    b = zeros( size( A , 1), 1); 
    
    %----------------------------------------------------------------------
    % Equality constraints to scale problem
    Aeq = zeros(3, 3*T);
    Aeq(1, norm) = 1;
    Aeq(2, 2*T + norm) = 1;
    Aeq(3, T + norm) = 1;
    beq = [0;0;1e-4];
    

    %----------------------------------------------------------------------
    % Build Objective   
    if smooth == 0;
        H = [zeros(3*T, 2*T) [zeros(2*T, T); eye(T) ]];      
    else
        tmpH1 = [zeros(3*T, 2*T) [zeros(2*T, T); eye(T) ]]; 
        
        tmpH2 = eye(T);
        for t = 2:T;
            tmpH2(t, t-1) = -1;
            if t < T;
                tmpH2(t, t+1) = -1;
            end;
        end;
        tmpH2(T/2 + 1, T/2) = 0;
        tmpH2(T/2, T/2 + 1) = 0;
        
        H = [zeros(3*T, 2*T) [zeros(2*T, T); 2*tmpH2]];
        H = 0.1*H + tmpH1;            
    end;
    
     f = zeros(3*T, 1);
     
    %----------------------------------------------------------------------
    % Bounds

    lb = [-Inf(T,1); (1e-10)*ones(T, 1); -Inf(T, 1)];
    ub = Inf*ones(3*T, 1);

    %----------------------------------------------------------------------
    % END MAIN PROGRAM
    %----------------------------------------------------------------------