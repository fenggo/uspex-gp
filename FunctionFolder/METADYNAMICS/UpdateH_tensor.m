function [lat_new1, par_new] = UpdateH_tensor(lat_old, pten, pres, width, height, par_old)
global ORG_STRUC
fpath = [ ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');
[S, C] = makeCtensor();
volume = det(lat_old);
av_pten = pten;
av_pten(1,1) = av_pten(1,1) - pres;
av_pten(2,2) = av_pten(2,2) - pres;
av_pten(3,3) = av_pten(3,3) - pres;
av_pten = av_pten + 2*rand(3,3)-1; 
dg_dh = volume * inv(lat_old) * av_pten;
f_c = MattoVec(dg_dh);
f_g = zeros(1,6);
lat_par = MattoVec(lat_old);
Step = size(par_old,1);
if Step > 1
for i = 1 : Step - 1
dp = lat_par - par_old(i,:);
gt = height * exp(-0.5*(norm(dp)^2)/width^2);
f_g = f_g + (dp/width^2)*gt;
end 
end
magn =300; 
if ORG_STRUC.correctionAngle  < max(abs(f_c(4:6)+f_g(4:6)))
ORG_STRUC.correctionAngle  = max(abs(f_c(4:6)+f_g(4:6)));
%  disp(['New correction amplitude (angle): ' num2str(ORG_STRUC.correctionAngle)]);
end
if ORG_STRUC.correctionLength <   max(abs(f_c(1:3)+f_g(1:3)))
ORG_STRUC.correctionLength = 3*max(abs(f_c(1:3)+f_g(1:3)));
%  disp(['New correction amplitude (length): ' num2str(ORG_STRUC.correctionLength)]);
end
badAngle = 0;
CorrectionA = 0; 
lat_tmp = latConverter(lat_old);
[minAngle, minA] = min(lat_tmp(4:6));
[minAngle1, minA1] = min(pi-lat_tmp(4:6));
if minAngle > minAngle1
minAngle = minAngle1;
minA = minA1;
end
if minAngle < 0.5236 
fprintf(fp, 'Worst angle: %6.3f degrees; correction will be performed.\n', minAngle);
badAngle = 1;
CorrectionA = magn*(1 - cos(0.5*pi*(0.5236 - minAngle)/0.5236));
end
CorrectionL1 = 0; 
CorrectionL2 = 0; 
lat_tmp = latConverter(lat_old);
[minLength, minL] = min(lat_tmp(1:3));
if minLength < ORG_STRUC.minVectorLength
fprintf(fp, 'Smallest length: %6.3f; correction will be performed.\n', minLength);
CorrectionL1 = magn*(1 - cos(0.5*pi*(ORG_STRUC.minVectorLength - minLength)/...
ORG_STRUC.minVectorLength));
end
[maxLength, maxL] = max(lat_tmp(1:3));
if maxLength > ORG_STRUC.maxVectorLength
fprintf(fp, ' Largest length: %6.3f; correction will be performed.\n', maxLength);
CorrectionL2 = magn*(1 - cos(0.5*pi*(ORG_STRUC.maxVectorLength - maxLength)/...
ORG_STRUC.maxVectorLength));
end
f_corrA = ones(1,6)*CorrectionA*ORG_STRUC.correctionAngle;
f_corrL = zeros(1,6);
f_corrL(minL) = f_corrL(minL) + CorrectionL1*ORG_STRUC.correctionLength; 
f_corrL(maxL) = f_corrL(maxL) - CorrectionL2*ORG_STRUC.correctionLength; 
f_t = f_g + f_c + f_corrA + f_corrL;
fnorm = norm(f_t);
f_t_m = VectoMat(f_t);
toAdd = MattoVec( tensorProduct(S, f_t_m/fnorm) * (lat_old / volume^(1/3)) )*width;
if min(abs(toAdd)) > width
disp('=========== Attention: changes to vector are too big! ============');
quit;
end
if max(abs(toAdd(1:3))) > width
toAdd = toAdd*width/max(abs(toAdd(1:3)));
end
lat_par1 = lat_par + toAdd;
lat_new1 = VectoMat(lat_par1);
if badAngle == 1 
f_t = f_g + f_c - f_corrA + f_corrL;
fnorm = norm(f_t);
f_t_m = VectoMat(f_t);
toAdd = MattoVec( tensorProduct(S, f_t_m/fnorm) * (lat_old / volume^(1/3)) )*width;
if max(abs(toAdd(1:3))) > width
toAdd = toAdd*width/max(abs(toAdd(1:3)));
end
lat_par2 = lat_par + toAdd;
lat_new2 = VectoMat(lat_par2);
lat_tmp1 = latConverter(lat_new1);
lat_tmp2 = latConverter(lat_new2);
min1 = min(lat_tmp1(3+minA), pi-lat_tmp1(3+minA));
min2 = min(lat_tmp2(3+minA), pi-lat_tmp2(3+minA));
disp(['Lattice1(+): ' num2str(lat_tmp1')]);
disp(['Lattice2(-): ' num2str(lat_tmp2')]);
disp(' ');
if min2 > min1   
lat_par1 = lat_par2;
lat_new1 = lat_new2;
f_corrA = -1*f_corrA;
end
end
par_new = [par_old ; lat_par1];
fprintf(fp, 'f_gaussian: %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', f_g);
fprintf(fp, 'f_tensor  : %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', f_t);
fprintf(fp, 'f_cor_Ang : %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', f_corrA);
fprintf(fp, 'f_cor_len : %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', f_corrL);
fprintf(fp, 'Addlattice: %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', toAdd);
fclose(fp);
