function [taste, change, exitflag] = static_model(pin, qin, good, smooth, norm)
%--------------------------------------------------------------------------
% PURPOSE: calculate taste change given assumption of static model
%--------------------------------------------------------------------------
    qin = qin(:, find(isnan(mean(qin))==0));
    pin = pin(:, find(isnan(mean(qin))==0));

    T =  size(pin, 2);
%     options = optimset('MaxIter', 5e5);
%     options = optimset('Display', 'Iter');
    if emax(pin, qin) == 1;
        taste = zeros(1, 2*T);
        change = 0;
        exitflag = 1;
    else 
        [H, f, A, b, Aeq, beq, lb, ub] = static_taste_constraints (pin, qin, good, smooth, norm);
        [result, change, exitflag] = quadprog(H, f , A, b, Aeq, beq, lb, ub, []);
        taste = result(T+1:3*T)';
    end;
    
end