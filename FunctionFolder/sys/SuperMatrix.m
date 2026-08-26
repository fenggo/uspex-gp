function Matrix = SuperMatrix(xmin, xmax, ymin, ymax, zmin, zmax)
Matrix = [];
for i = xmin:1:xmax
for j = ymin:1:ymax
for k = zmin:1:zmax
Matrix = [Matrix; [i,j,k]];
end
end
end
