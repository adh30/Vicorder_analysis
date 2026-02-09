clc
clearvars -except SBP DBP MAP HR SV_ref dataforSV
% ============================================================
% Compute SV estimates
% ============================================================
SV = unified_SV_models(SBP, DBP, MAP, HR);
SV_mod  = SV.modelflow;
SV_vic  = SV.vicorder;
SV_arc  = SV.arcsolver;
SV_lz  = SV.lz;
SV_picco=SV.picco;
SV_pp=SV.pp;

% ============================================================
% Regression helper function
% ============================================================
regfun = @(x,y) struct( ...
    'slope', polyfit(x,y,1), ...
    'r2',    corrcoef(x,y) );

% ============================================================
% Compute regressions for all pairwise combinations
% ============================================================
pairs = {
    'Reference','Vicorder', SV_ref, SV_vic;
    'Reference','Modelflow',   SV_ref, SV_mod;
    'Reference','ARCSolver', SV_ref, SV_arc;
    'Reference','LZ',  SV_ref, SV_lz;
    'Reference','PiCCO',  SV_ref, SV_picco;
    'Reference','PP',  SV_ref, SV_pp;
 };

fprintf('\n================ Regression Results ================\n');

for k = 1:size(pairs,1)
    name1 = pairs{k,1};
    name2 = pairs{k,2};
    x = pairs{k,3};
    y = pairs{k,4};

    % Fit regression
    p = polyfit(x,y,1);
    yfit = polyval(p,x);
    r = corrcoef(x,y);
    r2 = r(1,2)^2;

    % ICC
    stats = agree(x, y);

    % Print results
    fprintf('\n%s vs %s:\n', name1, name2);
    fprintf('   y = %.3f*x + %.3f\n', p(1), p(2));
    fprintf('   mean difference = %.4f\n', stats.mean_diff);
    fprintf('   SD difference = %.4f\n', stats.sd_diff);
    fprintf('   ICC = %.4f\n', stats.icc);
    fprintf('   r2 = %.4f\n', r2);
end


fprintf('\n===================================================\n');