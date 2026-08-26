function NEB_Done()
global ORG_STRUC
cd( ORG_STRUC.resFolder );
[nothing, nothing] = unix(['echo -e " " >>  VCNEBReports' ]);
[nothing, nothing] = unix(['echo -e "=====================================  VCNEB  DONE ===================================="   >> VCNEBReports']);
[nothing, nothing] = unix(['echo -e " " >>  VCNEBReports' ]);
