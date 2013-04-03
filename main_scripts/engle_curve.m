%--------------------------------------------------------------------------
% PURPOSE:  Kernel regression for NP Engle Curves
%           Control function approach to control for total exp. endogeneity
%           Double residual approach to control for demographics
%--------------------------------------------------------------------------
% OUTPUTS:  mhat      fitted budget share
%           ci        confidence interval
%           grid      income grid estimated over 
%           qout      quantitis before final NP regression
%--------------------------------------------------------------------------

function [mhat, ci, grid, qout] = engle_curve(p, win, x, inc, z, evaluate)

global g1 g2 g3 K;
    K = size(p, 1);
    N = length(x);
    lnx = log(x);
    sigma = sqrt(sum((lnx - mean(lnx)).^2)/(N - 1)); % s.deviation of lnx
    iqr = prctile(lnx, 75) - prctile(lnx, 25);
    shat = min(sigma, iqr/1.349);
    
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
    
    mhat = [];
    ci = [];
    q1 = []; q2 = [];
    for k = 1 : K-1;
    
    % 3. regress w, z and v_tilda on lnx 
        ww = w(:, k);
        covs = [ww z v_tilda];
        for j = 1 : 3;
            
            % 3a. kernel 
             y = covs(:,j);
             h = 0.9*shat*(N^(-0.2)); % Silverman density
             ghat = NaN(N,1);
             for ii = 1 : N;
                 u = lnx - lnx(ii);
                 u = u./h;
                 Ku = ((2*pi)^(-0.5))*exp(-0.5*(u.^2));
                 Ku = Ku(isnan(y)==0);
                 Kuy = Ku.*y(isnan(y)==0);
                 fxy = (1/(N*h))*sum(Kuy);
                 fx = (1/(N*h))*sum(Ku);
                 ghat(ii) = fxy/fx;
             end;

            if j == 1;
                w_hat = ghat;
            elseif j==2;
                z_hat = ghat;
            else
                v_hat = ghat;
            end;
            
         end; % end of 3
        
         % 4. coefficients on demo and v
         wresid = ww - w_hat;
         xresid = [(z - z_hat) (v_tilda - v_hat)];
         beta_resid = regress(wresid, xresid);
         
         % 5. kernel to find out demand at the smp income 
         y = w_hat - ones(N,1)*beta_resid(1);
         eval(['q' num2str(k) ' = y;']);
         h = 0.9*shat*(N^(-0.2));
         
         mhat_tmp = [];
         ci_tmp = [];
         if ~isempty(evaluate);
             xx = evaluate;
         else
             xx = min(lnx) : (max(lnx) - min(lnx))/1000 : max(lnx);
         end;
       
         for xxx = 1 : length(xx);
             u = lnx - xx(xxx);
             u = u./h;
             Ku = ((2*pi)^(-0.5))*exp(-0.5*(u.^2));
             Ku = Ku(isnan(y)==0);
             Kuy = Ku.*y(isnan(y)==0);
             fxy = (1/(N*h))*sum(Kuy);
             fx = (1/(N*h))*sum(Ku);
             ghat = fxy/fx;        
             mhat_tmp = [mhat_tmp (ghat + beta_resid(1))];
             s2 = mean((y(isnan(y)==0) - Kuy).^2);
             ci_tmp = [ci_tmp 1.96*sqrt((0.282094792*s2)/(N*h*mean((1/h)*Ku)))];
         end;
         mhat = [mhat; mhat_tmp];
         ci = [ci; ci_tmp];
   end; % end of good loop
   grid = min(lnx) : (max(lnx) - min(lnx))/1000 : max(lnx);
   qout = [q1 q2]; 
   mhat = [mhat; 1 - sum(mhat)];
end