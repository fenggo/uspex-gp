function RL_export_prior(RL)
% Export learned posterior to local RLPriors/operator_prior.mat (default on).
% Failures are non-fatal: caller is already inside a best-effort flow.
try
   dir_ = 'RLPriors';
   if ~exist(dir_, 'dir'), mkdir(dir_); end
   prior.alpha = RL.alpha;
   prior.beta  = RL.beta;
   prior.inner = RL.inner;
   prior.version = 1;
   safesave([dir_ '/operator_prior.mat'], prior);
catch err
   disp(['RL bandit: export failed (' err.message ')']);
end
end
