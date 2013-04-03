function e = emax(p,q);
% emax(p,q) the highest Afriat Efficiency parameter such
% that the data satisfy GARP
e = 1;
if garp(p,q)==0;
    d = 0.1;
    if garp(p,q)==0;
        while garpe(p,q,e)==0;
                e=e-d;
                garpe(p,q,e);
        end;
    end;
    
    e1 = e;
    e0 = e+d;
    e  = e1+d/2;
    while d>0.00000000001;
        if garpe(p,q,e)==1;
            e1=e;
            d = e0-e1;
            e  = e1+d/2;
        else;
            e0 = e;
            d = e0-e1;
            e  = e1-d/2;
        end;
    end;
    if garpe(p,q,e)==0;
        e=e-d/2;
    end;
end;
