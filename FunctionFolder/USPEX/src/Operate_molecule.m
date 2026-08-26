function R = Operate_molecule(R0, Opt, lat, P, PB)
R0_frac = Cart2Frac(R0, lat);
for axis = 1:3
R0_frac_p(:,axis) = R0_frac(:,P(axis));
lat_p = lat(P(axis));
end
for i = 1:size(R0,1)
R_frac_p(i,:)  = R0_frac_p(i,:)*Opt(1:3,:) + Opt(4,:);
end
for axis = 1:3
R_frac(:,axis) = R_frac_p(:,PB(axis));
end
R   = Frac2Cart(R_frac, lat);
