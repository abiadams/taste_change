function pout = form_prices_fes(pin, win, xin, K)
%--------------------------------------------------------------------------
% PURPOSE: create consolidated prices using chained Laspayres index
%--------------------------------------------------------------------------
    global g1 g2 g3;

    expnds = win.*(xin*ones(1,size(win,2)));
    expnds = expnds.*(expnds>=0);
    P = size(pin, 1);

    for k = 1 : K;
        eval(['good = g' num2str(k) ' ;']);
        xk = sum(expnds(:, good));
        pp = sum((pin(:,good).*(ones(P,1)*(xk./sum(xk)))),2);
        ww = sum(expnds(:, good), 2)./sum(xk);
        pp = sum(pp.*ww);
        if k == 1;
            pout = pp;
        else
            pout = [pout; pp];
        end;       
    end;
end