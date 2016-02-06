//basedir = "C:\\Users\\Martin\\Desktop\\cs-experiments-2014-08-18\\!! EXPERIMENTS !!\\astigmatism\\SNR+density";
basedir = "D:\\!CS PAPER\\cs-2014-experiments\\astigmatism\\SNR+density";
methods = newArray("SA-WLSQ","daostorm","CS-MLE");
bkg = 70.0;
tol_export_dist = 200;

run("Camera setup", "isemgain=false pixelsize=80.0 offset=100.0 photons2adu=1.0");

for(i = 0; i < methods.length; i++) {
	method = methods[i];
	
	density = 0.1;
	for(intensity = 2500; intensity <= 2500; intensity += 1000) {
		run("Import ground-truth", "append=false startingframe=1 filepath=["+basedir+"\\density="+density+"\\I="+intensity+"+bkg="+bkg+".csv] fileformat=[CSV (comma separated)]");
		run("Import results", "append=false startingframe=1 rawimagestack= filepath=["+basedir+"\\density="+density+"\\results\\"+method+"\\I="+intensity+"+bkg="+bkg+".csv] livepreview=false fileformat=[CSV (comma separated)]");
		wait(200);	// wait to show both windows (it is necessary to run the evaluation)
		for(tol = 50; tol <= 1000; tol += 50) {
			run("Performance evaluation", "toleranceradius="+tol+"  evaluationspace=xy");
			if(tol == tol_export_dist) {
				wait(200);
				run("Export results", "z=false filepath=["+basedir+"\\density="+density+"\\results\\"+method+"+I="+intensity+"+bkg="+bkg+"_tol="+tol+".csv] fileformat=[CSV (comma separated)] sigma_x=false id=false sigma_y=false frame=true gt_dist_xy=true gt_dist_z=true gt_dist_xyz=true gt_id=true saveprotocol=false offset=false i=false y=false x=false sigma1=false sigma2=false intensity=false bkgstd=false chi2=false uncertainty=false z_rel=false");
			}
		}
		saveAs("Measurements", basedir+"\\density="+density+"\\results\\"+method+"+I="+intensity+"+bkg="+bkg+".xls");
		run("Clear Results");
	}
	
	for(density = 0.5; density <= 20.0; density += 0.5) {
		for(intensity = 2500; intensity <= 2500; intensity += 1000) {
			run("Import ground-truth", "append=false startingframe=1 filepath=["+basedir+"\\density="+density+"\\I="+intensity+"+bkg="+bkg+".csv] fileformat=[CSV (comma separated)]");
			run("Import results", "append=false startingframe=1 rawimagestack= filepath=["+basedir+"\\density="+density+"\\results\\"+method+"\\I="+intensity+"+bkg="+bkg+".csv] livepreview=false fileformat=[CSV (comma separated)]");
			wait(200);	// wait to show both windows (it is necessary to run the evaluation)
			for(tol = 50; tol <= 1000; tol += 50) {
				run("Performance evaluation", "toleranceradius="+tol+"  evaluationspace=xy");
				if(tol == tol_export_dist) {
					wait(200);
					run("Export results", "z=false filepath=["+basedir+"\\density="+density+"\\results\\"+method+"+I="+intensity+"+bkg="+bkg+"_tol="+tol+".csv] fileformat=[CSV (comma separated)] sigma_x=false id=false sigma_y=false frame=true gt_dist_xy=true gt_dist_z=true gt_dist_xyz=true gt_id=true saveprotocol=false offset=false i=false y=false x=false sigma1=false sigma2=false intensity=false bkgstd=false chi2=false uncertainty=false z_rel=false");
				}
			}
			saveAs("Measurements", basedir+"\\density="+density+"\\results\\"+method+"+I="+intensity+"+bkg="+bkg+".xls");
			run("Clear Results");
		}
	}
}
