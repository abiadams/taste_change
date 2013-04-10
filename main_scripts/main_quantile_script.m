%--------------------------------------------------------------------------
% PURPOSE: Calculate taste change on tobacco associated with different 
%           age and education cohorts at different quantiles
%--------------------------------------------------------------------------
clear; clc;

% Set paths  
addpath('H:\DPhil\Tastes-Tobacco\Testing\MT12');
addpath('H:\DPhil\Tastes-Tobacco\Testing\HT13');
addpath('H:\DPhil\Tastes-Tobacco\Testing\TT13');
addpath('H:\DPhil\Data and Code\Data');
addpath('H:\DPhil\Data and Code\Functions');

% Load data
load 'H:\DPhil\Tastes-Tobacco\Testing\TT13\quantile_input';

% Problem parameters
good = 1;
norm= 1;
smooth = 0;
K = 2;
T = length(time);
cells = size(all_prices, 1)/K;

p = []; q = []; t=[];
for c = 1: cells;
    t = [t time'];
    p = [p all_prices(c*K - (K-1):c*K, 3:T+2)];
    q = [q all_quantities(c*K - (K-1):c*K, 3:T+2)];
end;
z = find(q(:, 1) == 0);
look_at = setdiff(1:size(p, 2), z);
pp = p(:, look_at);
qq = q(:, look_at);
tt= t(look_at);

disp('Now looking at taste change');
[taste, change, exitflag] = static_model(pp, qq, good, 0, norm);

save quantile_results

