//basedir = "C:\\Users\\Martin\\Desktop\\cs-experiments-2014-08-18\\!! EXPERIMENTS !!\\astigmatism\\";
basedir = "D:\\!CS PAPER\\cs-2014-experiments\\astigmatism\\";
frames = 100;
bkg = 70.0;

run("Camera setup", "isemgain=false pixelsize=80.0 offset=100.0 photons2adu=1.0");

density = 0.1;
for (intensity = 2500; intensity <= 2500; intensity += 1000) {
  run("Generator of simulated data", "z_range=-400:+400 calibration=["+basedir+"daostorm3d.yaml] psf=[Eliptical Gaussian (3D astigmatism)] maskpath=["+basedir+"\\SNR+density\\square-mask-32x32--20x20.png] height=32 width=32 density="+density+" frames="+frames+" addpoisssonvar="+bkg+" driftangle=0.0 driftdist=0.0 intensityrange="+intensity+" singlefixed=false");
  saveAs("Tiff", basedir + "SNR+density\\density="+density+"\\I="+intensity+"+bkg="+bkg+".tif");
  close();
  run("Export ground-truth", "id=false frame=true sigma1=true filepath=["+basedir+"SNR+density\\density="+density+"\\I="+intensity+"+bkg="+bkg+".csv] sigma2=true bkgstd=true intensity=true offset=true z=true angle=true y=true x=true fileformat=[CSV (comma separated)]");
}

for (density = 0.5; density <= 20.0; density += 0.5) {
  for (intensity = 2500; intensity <= 2500; intensity += 1000) {
    run("Generator of simulated data", "z_range=-400:+400 calibration=["+basedir+"daostorm3d.yaml] psf=[Eliptical Gaussian (3D astigmatism)] maskpath=["+basedir+"SNR+density\\square-mask-32x32--20x20.png] height=32 width=32 density="+density+" frames="+frames+" addpoisssonvar="+bkg+" driftangle=0.0 driftdist=0.0 intensityrange="+intensity+" singlefixed=false");
    saveAs("Tiff", basedir+"SNR+density\\density="+density+"\\I="+intensity+"+bkg="+bkg+".tif");
    close();
    run("Export ground-truth", "id=false frame=true sigma1=true filepath=["+basedir+"SNR+density\\density="+density+"\\I="+intensity+"+bkg="+bkg+".csv] sigma2=true bkgstd=true intensity=true offset=true z=true angle=true y=true x=true fileformat=[CSV (comma separated)]");
  }
}
