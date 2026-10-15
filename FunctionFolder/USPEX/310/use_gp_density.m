function density = use_gp_density(POP, density)
% Use the GP-predicted density when available, else keep the relaxed (actual) density.
%   density = use_gp_density(POP, density)
% POP     : one POP_STRUC.POPULATION structure
% density : relaxed density computed from the structure (fallback)
% Rule: if POP.density_gp exists and is a finite scalar, use it; otherwise keep density.
try
   if isfield(POP, 'density_gp')
      dg = POP.density_gp;
      if ~isempty(dg) && isnumeric(dg) && isscalar(dg) && ~isnan(dg) && isfinite(dg)
         density = dg;
      end
   end
catch
   % any unexpected field shape -> silently keep the relaxed density
end
end
