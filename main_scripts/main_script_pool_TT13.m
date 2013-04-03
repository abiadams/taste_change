%--------------------------------------------------------------------------
% PURPOSE: Integrated taste change script
%--------------------------------------------------------------------------
% Programme 1:  Static
%           2:  Fixed Lambda
%--------------------------------------------------------------------------
% Specify   1:  Adaptive kernel?
%           2:  SMP/Non-SMP
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% DATA INPUT
%--------------------------------------------------------------------------
% Set paths 
clear; clc; 
addpath('H:\DPhil\Tastes-Tobacco\Testing\MT12');
addpath('H:\DPhil\Tastes-Tobacco\Testing\TT13');
addpath('H:\DPhil\Data and Code\Data');
addpath('H:\DPhil\Data and Code\Functions');

% Load data
load 'H:\DPhil\Data and Code\Data\boe_quarterly_base';
    interest_rate_data = data;

load 'H:\DPhil\Data and Code\Data\cormacdata';
    age = data(:,107);
    period = data(:,3); 
    month = data(:,78);
    cohort = period-age;
    totx = data(:,4);
    income = data(:,112);
    education = data(:,101);
    hheqsize = data(:,113);
    budget_shares = data(:, 7:77);
    prices = data(:, 114:184)/100;

%--------------------------------------------------------------------------
% GOODS
%--------------------------------------------------------------------------
durables = [31:37 42:45 51:56 59:60 62 66 69:70];
nondurables = setdiff(1:71, durables);
food = 1:26;
tobacco = 29:30;
alcohol = 27:28;

% set good constants
global g1 g2 g3 K;
K = 3;          % number of goods
g1 = tobacco;
g2 = alcohol;
g3 = nondurables;

%--------------------------------------------------------------------------
% COHORT & TIME
%--------------------------------------------------------------------------
% cohort time breaks
breaks = [1955 1960];
global C;
C = length(breaks) - 1;

% smoker 0) all 1) smoker cohort
smoker = 0;

% min and max ages
min_age = 20;
max_age = 60;

% choose time unit: 1) annual 2) quarterly
time_unit = 1;

% set start and end dates for analysis: earliest date: 1978, latest: 2009
start_year = 1980;
end_year = 2000;

if time_unit == 2;
    T = end_year - start_year + 1;
    time = start_year:end_year;
    quarter = ones(size(data, 1), 1);
    for i = 1:size(data, 1);
        if month(i,1) < 4;
            quarter(i,1) = 1;
        elseif month(i, 1) >=4 && month(i,1) <= 6;
            quarter(i, 1) = 2;
        elseif month(i, 1) >=7 && month(i,1) <= 9;
            quarter(i, 1) = 3;
        else
            quarter(i, 1) = 4;
        end;
    end;
    period_quart =[period quarter];
    time = period_quart(period_quart(:, 1) >= start_year, :);
    time = time(time(:, 1) <= end_year, :);
    time = unique(time, 'rows');
else
    time = start_year:end_year;
    time = time';
end;
T = size(time, 1);
    
% education-year cutoff
comp_ed_change = 1958;

%--------------------------------------------------------------------------
% MAIN PROGRAMME
%--------------------------------------------------------------------------
% PROG  1 : STATIC
%       2 : FIXED LAMBDA
%--------------------------------------------------------------------------
PROG = 1;
ADAPT = 1;
SMP=1;
good = 1;
smooth = 0;
norm = 1;


for age_cohort = 1 : C;
    ed_q = [];
    ed_p = [];
    for ed = 1 : 2;
        display(['Now running for age_cohort ' int2str(age_cohort)]);
        display(['                education level ' num2str(ed)]);
        
        % select relevant data
        s = find(   cohort >= breaks(age_cohort) & ...
                    cohort < breaks(age_cohort + 1) & ...
                    income > 0 & ...
                    totx > 0 & ...
                    age >= min_age & ...
                    age <= max_age & ...
                    period >= start_year & ...
                    period <= end_year );
                
        if smoker == 1;
            smokers =  find(sum(budget_shares(s, g1), 2) > 0);
            s = s(smokers);
        end;
                
        if ed == 1 && breaks(age_cohort + 1) > comp_ed_change;
            ee = find(education(s) <= 16 );
        elseif ed == 1 && breaks(age_cohort + 1) <= comp_ed_change;
            ee = find(education(s) <= 15 );
        elseif ed == 2 && breaks(age_cohort + 1) > comp_ed_change;
            ee = find(education(s) > 16 );
        else
            ee = find(education(s) > 15 );
        end; % end of selecting education statement
        
        cohort_people = s(ee);
        N = ones(1, length(time));
           
        % start on time loop
        q = NaN(K,T); p = NaN(K,T); ir = NaN(T,1);
        for t = 1 : T;
            % select data           
            if time_unit == 1;
                time_people = find( period(cohort_people, 1) == time(t));
                time_people = cohort_people(time_people);
                time_i = mean(interest_rate_data(interest_rate_data(:,1) == time(t), 3));
            else
                time_people = find( period_quart(cohort_people, 1) == time(t, 1) & ...
                                period_quart(cohort_people, 2) == time(t, 2));
                time_people = cohort_people(time_people);
                t_unit = find(interest_rate_data(:,1) == time(t, 1) & ...
                              interest_rate_data(:,2) == time(t, 2));
                time_i = mean(interest_rate_data(t_unit, 3));
                time_i = time_i.^(1/4);
            end;
     
            N(1, t) = length(time_people);
            time_p = prices(time_people, :);
            time_w = budget_shares(time_people, :);
            time_x = totx(time_people);
            time_inc = income(time_people);
            time_hheq = hheqsize(time_people);

            % form prices
            ptmp = form_prices_fes(time_p, time_w, time_x, K);
            
            % form quantities
            if SMP == 0;
                evaluate =  median(time_x);     
            else
                if t==1;
                     evaluate =  median(time_x);
                else
                     evaluate =  ptmp'*q(:,t-1);
                end;
            end;
            ln_evaluate =  log(evaluate);
            
            if ADAPT == 1; 
                qtmp = engle_curve_adaptive(ptmp, time_w, time_x, time_inc, time_hheq, ln_evaluate);
            else
               qtmp = engle_curve(ptmp, time_w, time_x, time_inc, time_hheq, ln_evaluate);
            end;  
            qtmp = (evaluate*qtmp)./ptmp;
            
            if sum(isnan(qtmp))>0;
                error('quantities screwy')
            end;
            
            p(:,t) = ptmp;
            q(:,t) = qtmp;
            ir(t) = time_i;
            
        end; % end of time loop
        ed_q = [ed_q q]; 
        ed_p = [ed_p p];
    end; % end of education loop
    
    %---------------------------------------------------------------------
    % TASTE CHANGE WITHIN COHORTS
    %---------------------------------------------------------------------    
    if PROG == 1;
         [taste, change, exitflag] = static_model(ed_p, ed_q, good, 0, norm);
    elseif PROG == 2;
        [taste, change, exitflag] = fixed_lambda_pool_model(ed_p, ed_q, good, smooth, norm, ir, time);
    end;

    lambda = [taste(1:T); taste(T+1: 2*T)];
    alpha = [taste(2*T+1:3*T); taste(3*T+1: 4*T)];
    results = [(p(good, :) - (alpha(1,:)./lambda(1,:)))./p(K, :); ...
                (p(good, :) - (alpha(2,:)./lambda(2,:)))./p(K, :)];
            
    % smoothed results 
    % impose monotonicity/smoothness with cubic hermite spline
    mono = ones(2, T);
    results_hermite = NaN(2, T); 
    for c = 1 : 2;
        for t = 2 : T;
            if any(results(c, 1:t-1) > results(c, t));
                mono(c, t) = 0;
            end;
        end;
        results_hermite(c, :) = pchip(time(mono(c,:) == 1), results(c, mono(c,:) == 1), time');
    end;
    
    figure;
    plot(time, results_hermite(1, :)', '-r' , 'LineWidth', 2);
    hold on 
    plot(time, results_hermite(2, :)', '-b' , 'LineWidth', 2);
    plot(time, (p(1,:)./p(3,:)), ':g', 'LineWidth', 2);
    legend('Taste Adjusted Price: Low Education', 'Taste Adjusted Price: High Education', 'Actual Price');
 
end; % end of cohort age loop


