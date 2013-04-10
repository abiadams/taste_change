%--------------------------------------------------------------------------
% PURPOSE:  Multivariate spline for NP Engle Curves
%           Control function approach to control for total exp. endogeneity
%           Censoring correction a la Lewbel
%--------------------------------------------------------------------------
% OUTPUTS:  mhat      fitted budget share
%--------------------------------------------------------------------------

function mhat = engle_curve_quantile_censor(p, win, x, inc, z, tau, evaluate)

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
    
    % 2. Estimate control variable
    instruments = [ones(N,1) log(inc) z];
    pi_hat = (instruments'*instruments)\instruments'*lnx;
    residual = lnx - instruments*pi_hat;
    [F, xx] = ecdf(residual);
    V_hat = NaN(N,1);
    h = 0.9*min(sqrt(sum((residual - mean(residual)).^2)/(N - 1)),(prctile(residual, 75) - prctile(residual, 25))/1.349)*(N^(-0.2));
    for ii = 1 : N;
         u = xx - residual(ii);
         u = u./h;
         Ku = ((2*pi)^(-0.5))*exp(-0.5*(u.^2));
         Kuy = Ku.*F;
         fxy = (1/(N*h))*sum(Kuy);
         fx = (1/(N*h))*sum(Ku);
         V_hat(ii) = fxy/fx;
    end;
    
    % 3. estimate censoring indicator
    X_ = [lnx lnx.^2 z V_hat];
    X = [ones(N,1) X_];
    y = w(:,1);
    gamma = glmfit(X_, (y == 0), 'binomial');
    delta_hat = glmval(gamma, X_, 'logit');
    delta_hat = 1 - delta_hat;
    
    % initial sample selection
    j0 = find(delta_hat > 1 - tau + 0.1);
    
    % 4. initial inefficient estimator
    rho = @(r) sum(abs(r - (r<0).*r/tau));
    init_ols = X(j0,:)\y(j0,:);
    beta_0 = fminsearch(@(p)rho(y(j0,:) - X(j0,:)*p), init_ols);
    
    % 5. select observations again
    for i = 1 : 3;
        eval(['jtmp = find(X*beta_' num2str(i-1) ' > 1e-4);']);
        init_ols = X(jtmp,:)\y(jtmp,:);
        eval(['beta_' num2str(i) ' = fminsearch(@(p)rho(y(jtmp,:) - X(jtmp,:)*p), init_ols);']);
    end;
    
    % 6. Predicted value at evaluate
    wtmp = [1 evaluate evaluate.^2 1 0.5]*beta_3;
    mhat = [wtmp; 1 - wtmp]; 
end