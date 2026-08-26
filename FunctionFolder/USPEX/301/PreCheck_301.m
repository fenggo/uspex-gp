function PreCheck_301()
global ORG_STRUC
if size(ORG_STRUC.numIons,1) == 1
error('==> USPEX no longer supports single block in 301, please use 300 instead! <==');
end
