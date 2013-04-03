function result = garp (p, q)
%--------------------------------------------------------------------------
% PURPOSE: garp (p, q) tests the Generalized Axiom of Revealed Preference
% (GARP) and returns 0 for GARP fails or 1 for GARP passes.
%--------------------------------------------------------------------------
% USAGE: result = garp (p, q)
% where: p = price matrix
%        q = quantity matrix
%--------------------------------------------------------------------------
% RETURNS: result = 0 or 1
%--------------------------------------------------------------------------
% REFERENCE: Varian, H. R. (1982): "The Nonparametric Approach to Demand
%            Analysis," *Econometrica*, 50, 945-973.
%--------------------------------------------------------------------------

    %----------------------------------------------------------------------
    % CHECK USER ERRORS
    %----------------------------------------------------------------------

    % check number of arguments
    if nargin ~= 2;
        error ('Not enough input arguments.');
    end;

    % check price/quantity consistency
    if sum (size (q) == size (p), 2) ~= 2;
        error ('Price and quantity dimensions must agree.');
    end;

    % check minimum number of goods
    if size (q, 1) < 2;
        error ('Not enough goods.');
    end;

    % check minimum number of periods
    if size (q, 2) < 2;
        error ('Not enough periods.');
    end;

    %----------------------------------------------------------------------
    % BEGIN MAIN PROGRAM
    %----------------------------------------------------------------------

    % initialize revealed preference matrix
    M = eye (size (q, 2));

    % summarize directly revealed preferences
    for i = 1 : size (q, 2);
        for j = 1 : size (q, 2);
            M(i, j) = (p(:, i)' * q(:, i) >= p(:, i)' * q(:, j));
        end;
    end;

    % summarize revealed preferences
    for i = 1 : size (q, 2);
        for j = 1 : size (q, 2);
            M(i, :) = max ([M(i, :); M(i, j) * M(j, :)]);
        end;
    end;

    % test GARP
    result = 1;
    for i = 1 : size (q, 2);
        for j = 1 : size (q, 2);
            if M(i, j) == 1;
                if p(:, j)' * q(:, j) > p(:, j)' * q(:, i);
                    result = 0;
                end;
            end;
        end;
    end;

    %----------------------------------------------------------------------
    % END MAIN PROGRAM
    %----------------------------------------------------------------------