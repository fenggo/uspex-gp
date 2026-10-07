function f = RL_agent(popSize)
% RL_agent  Hierarchical contextual bandit over variation operators (L2).
%   f = RL_agent(popSize)
%   Settles the just-finished generation, then draws next-gen arm fractions.
%   State persisted in rl_agent_state.mat; priors read from RLPriors/ (local then global);
%   posterior exported to RLPriors/ when RL_autoExport is on.
%   Any failure is caught by the caller (update_STUFF) and downgraded.
global POP_STRUC
global ORG_STRUC

stateFile = 'rl_agent_state.mat';
PRI = RL_agent_init();          % arm/inner definitions + bounds

% ---------- load or initialise state ----------
RL = [];
if exist(stateFile, 'file')
   try
      L = load(stateFile);  RL = L.RL;
      if ~isfield(RL,'version') || RL.version ~= 1, RL = []; end
   catch
      RL = [];
   end
end
if isempty(RL)
   RL = RL_fresh(PRI);
   RL = RL_fuse_prior(RL, PRI);   % cross-task prior, if any
end

lam   = ORG_STRUC.RL_forget;
w1    = ORG_STRUC.RL_w1;
w2    = ORG_STRUC.RL_w2;
armCap= ORG_STRUC.RL_armCap;
kap   = ORG_STRUC.RL_kappa;

% ---------- settle current generation ----------
N = length(POP_STRUC.POPULATION);
% good-enthalpy baseline for relative gain reward
Hgood = [];
for i = 1:N
   H = POP_STRUC.POPULATION(i).Enthalpies(end);
   if H < 90000, Hgood(end+1,1) = H; end
end
if isempty(Hgood), Hm = 0; Hs = 1; else, Hm = mean(Hgood); Hs = std(Hgood); end
if ~(Hs > 0), Hs = 1; end

% elite set = top of current ranking (size used by update_STUFF)
Npro = max(1, round(N * ORG_STRUC.bestFrac));
elite = POP_STRUC.ranking(1:min(Npro, length(POP_STRUC.ranking)));

trials = zeros(1, PRI.K);
wins   = zeros(1, PRI.K);
uUsed  = RL.lastU;                    % representative params chosen for THIS gen
for i = 1:N
   arm = RL_which_arm(POP_STRUC.POPULATION(i).howCome, PRI);
   if arm == 0, continue; end
   H = POP_STRUC.POPULATION(i).Enthalpies(end);
   r = 0;
   if H < 90000
      inElite = any(elite == i);
      gain    = min(1, max(0, (Hm - H)/Hs));
      r = w1 * double(inElite) + w2 * gain;
   end
   trials(arm) = trials(arm) + 1;
   wins(arm)   = wins(arm) + r;
   % inner posterior update for the parameter that was used by this arm
   if r > 0 && any(arm == [4 5 6])
      for k = 1:size(uUsed,2)
         u = uUsed(arm,k);
         if ~isnan(u)
            nold = lam * RL.inner(arm,k).neff;
            mold = RL.inner(arm,k).m;
            nnew = nold + r;
            RL.inner(arm,k).m    = (mold*nold + u*r) / nnew;
            RL.inner(arm,k).neff = nnew;
         end
      end
   end
end
% outer posterior update
RL.alpha = lam * RL.alpha + wins;
RL.beta  = lam * RL.beta  + (trials - wins);
RL.alpha = max(RL.alpha, 1e-3);
RL.beta  = max(RL.beta,  1e-3);

% ---------- choose next generation ----------
% draw popSize Thompson samples across arms
count = zeros(1, PRI.K);
for s = 1:popSize
   theta = zeros(1, PRI.K);
   for a = 1:PRI.K
      x = randg(RL.alpha(a)); y = randg(RL.beta(a));
      theta(a) = x/(x+y);
   end
   [~, best] = max(theta);
   count(best) = count(best) + 1;
end
f = count / popSize;

% key-arm floors then cap/floor-consistent renormalisation with redistribution
f(1) = max(f(1), 0.10);   % heredity
f(2) = max(f(2), 0.10);   % random
f(6) = max(f(6), 0.10);   % softmutate
f = RL_safe_fractions(f, armCap);

% inner continuous Thompson sampling -> physical parameters for arms 4/5/6
nextU = nan(PRI.K, 2);
for a = [4 5 6]
   for k = 1:PRI.np(a)
      sd = kap / sqrt(RL.inner(a,k).neff) + PRI.innerBaseSD;
      u  = RL.inner(a,k).m + sd * randn();
      nextU(a,k) = min(1, max(0, u));
   end
end
RL.lastU = nextU;

% translate chosen u -> ORG_STRUC parameters
u = nextU(4,1);
ORG_STRUC.rotationAngleMax = PRI.angleLo + u * (PRI.angleHi - PRI.angleLo);
u = nextU(4,2);
ORG_STRUC.translationMax  = PRI.transLo + u * (PRI.transHi - PRI.transLo);
u = nextU(5,1);
ORG_STRUC.mutationRate    = PRI.mutLo + u * (PRI.mutHi - PRI.mutLo);
u = nextU(6,1);
ORG_STRUC.softStepScale   = exp(PRI.logScaleLo + u * (PRI.logScaleHi - PRI.logScaleLo));

% ---------- persist + export ----------
RL.history(end+1).gen    = POP_STRUC.generation;
RL.history(end).trials   = trials;
RL.history(end).wins     = wins;
RL.history(end).f        = f;
RL.history(end).u        = nextU;
RL.version = 1;
safesave(stateFile, RL);
if ORG_STRUC.RL_autoExport
   RL_export_prior(RL);
end

disp(['RL bandit: f = [' sprintf(' %.2f', f) ' ]']);
end

% ================= helpers =================
function PRI = RL_agent_init()
PRI.K = 6;
PRI.keys = {'Heredity','Random','Permutate','Rotate','LatMutate','softmutate'};
PRI.np = [0 0 0 2 1 1];                 % number of continuous params per arm
% physical bounds
PRI.angleLo = pi/12; PRI.angleHi = pi/2;
PRI.transLo = 0.1;   PRI.transHi = 0.8;
PRI.mutLo   = 0.05;  PRI.mutHi   = 0.30;
PRI.logScaleLo = log(0.5); PRI.logScaleHi = log(2.0);
PRI.innerBaseSD = 0.05;
end

function arm = RL_which_arm(howCome, PRI)
arm = 0;
if isempty(howCome), return; end
for a = 1:PRI.K
   if ~isempty(strfind(howCome, PRI.keys{a}))
      arm = a; return;
   end
end
end

function RL = RL_fresh(PRI)
RL.alpha = ones(1, PRI.K);
RL.beta  = ones(1, PRI.K);
RL.inner(PRI.K,2).m = 0.5;
for a = 1:PRI.K
   for k = 1:2
      RL.inner(a,k).m = 0.5;
      RL.inner(a,k).s = 0.25;
      RL.inner(a,k).neff = 1;
   end
end
% defaults aligned with current hardcodes
RL.inner(4,1).m = 1.0;                  % rotationAngleMax = pi/2
RL.inner(4,2).m = (0.5-0.1)/(0.8-0.1);  % translationMax  = 0.5
RL.inner(5,1).m = 0.5;                  % mutationRate midpoint
RL.inner(6,1).m = 0.5;                  % softStepScale = 1
RL.lastU = nan(PRI.K,2);
RL.history = [];
end

function RL = RL_fuse_prior(RL, PRI)
[prior, src] = RL_read_prior_file();
if isempty(prior), return; end
pw = 1.0;   % pseudo-count weight for the prior
try
   RL.alpha = RL.alpha + pw * max(0, prior.alpha - 1);
   RL.beta  = RL.beta  + pw * max(0, prior.beta  - 1);
   for a = 1:PRI.K
      for k = 1:2
         if isfield(prior,'inner') && length(prior.inner) >= a
            p = prior.inner(a,k);
            RL.inner(a,k).m    = p.m;
            RL.inner(a,k).neff = RL.inner(a,k).neff + pw * p.neff;
         end
      end
   end
   disp(['RL bandit: loaded cross-task prior from ' src]);
catch err
   disp(['RL bandit: prior ignored (' err.message ')']);
end
end

function [prior, src] = RL_read_prior_file()
prior = []; src = '';
name = 'operator_prior.mat';
cand = { ['RLPriors/' name], [get_home_dir() '/.config/uspex-gp/RLPriors/' name] };
for i = 1:length(cand)
   if exist(cand{i}, 'file')
      try
         Q = load(cand{i});
         if isfield(Q,'prior') && isfield(Q.prior,'alpha')
            prior = Q.prior; src = cand{i}; return;
         end
      catch
      end
   end
end
end

function h = get_home_dir()
h = getenv('HOME');
if isempty(h), h = char(javaMethod('getProperty','java.lang.System','user.home')); end
end

function f = RL_safe_fractions(f, cap)
% Return fractions that sum to 1, respect key-arm floors already applied, and keep
% every arm at <= cap by repeatedly redistributing the excess over uncapped arms.
f = f / sum(f);
for iter = 1:50
   over = f > cap + 1e-12;
   if ~any(over), break; end
   excess = sum(f(over) - cap);
   f(over) = cap;
   under = ~over;
   room = sum(f(under));
   if room <= 0
      f(1) = f(1) + excess; break;   % degenerate; floors dominate
   end
   f(under) = f(under) + excess * (f(under)/room);
end
f = f / sum(f);
end
