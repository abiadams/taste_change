%--------------------------------------------------------------------------
% PURPOSE:  Multivariate spline for NP Engle Curves
%           Control function approach to control for total exp. endogeneity
%           Censoring correction a la Lewbel
%--------------------------------------------------------------------------
% OUTPUTS:  mhat      fitted budget share
%--------------------------------------------------------------------------

function mhat = engle_curve_spline_censor(p, win, x, inc, z, evaluate)

    global g1 g2 g3 K;
    K = size(p, 1);
    N = length(x);
    lnx = log(x);
    
    % 1. form budget shares
    expnds = win.*(x*ones(1, size(win, 2)));
    expnds = expnds.*(expnds>=0);
    w = NaN(N, K);
    for k = 1 : K;
        eval(['good = g' num2str(k) ' ;']);
        w(:,k) = sum(expnds(:, good), 2)./sum(expnds(:, [g1 g2 g3]), 2);
    end;
    
    % 2. OLS of log(x) on log(m) to get residuals
    instruments = [ones(N,1) log(inc) z];
    pi_hat = (instruments'*instruments)\instruments'*lnx;
    v_tilda = lnx - instruments*pi_hat;
    
    X = [lnx z v_tilda];
    mhat = [];
    for k = 1 : K-1;
        y = w(:, k);
        model_r = aresbuild(X, y);       
        rhat = arespredict(model_r, [lnx z v_tilda]);
        r = arespredict(model_r, [evaluate 0.67 0]);
        
        if k == 1;
            lambda = max(rhat);
            h = 0.9*min(sqrt(sum((rhat - mean(rhat)).^2)/(N - 1)),(prctile(rhat, 75) - prctile(rhat, 25))/1.349)*(N^(-0.2));
            integrate_r = r:(lambda - r)/100:lambda;
            qhat = NaN(length(integrate_r),1);
            for ii = 1 : length(integrate_r);
                 u = rhat - integrate_r(ii);
                 u = u./h;
                 Ku = ((2*pi)^(-0.5))*exp(-0.5*(u.^2));
                 Kuy = Ku.*(w(:,k) > 0);
                 fxy = (1/(N*h))*sum(Kuy);
                 fx = (1/(N*h))*sum(Ku);
                 qhat(ii) = fxy/fx;
            end;
            inv_qhat = ones(length(qhat),1)./qhat;
            integral = trapz(integrate_r, inv_qhat);
            mhat_tmp = lambda - integral;
        else
            mhat_tmp = r;
        end;
        mhat = [mhat; mhat_tmp];
        
    end; % end of good loop
    mhat = [mhat; 1 - sum(mhat)];
end