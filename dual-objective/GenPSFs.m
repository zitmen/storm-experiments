function I = GenPSFs(cal,gx,gy,x0,y0,z0)

    [sx,sy] = defocusGaussian(cal,z0);
    I1 = genGaussian(gx(:,:,1),gy(:,:,1),x0,y0,sx(1),sy(1),cal.fi);
    I2 = genGaussian(gx(:,:,2),gy(:,:,2),x0,y0,sx(2),sy(2),cal.fi);
    I = [I1,I2];

end

function [wx,wy] = defocusGaussian(cal,z)
    wx = zeros(length(cal.cx),length(z));
    wy = zeros(length(cal.cy),length(z));
    for zi=1:length(z)
        wx(:,zi) = (cal.w0./cal.px./2) .* sqrt(1 + ((z(zi) - cal.cx)./cal.d).^2);
        wy(:,zi) = (cal.w0./cal.px./2) .* sqrt(1 + ((z(zi) - cal.cy)./cal.d).^2);
    end
end

function psf = genGaussian(xgrid,ygrid,x0,y0,xsigma,ysigma,fi)
    sinfi = sin(fi);
    cosfi = cos(fi);
    dx = ((xgrid - x0).*cosfi - (ygrid - y0).*sinfi);
    dy = ((xgrid - x0).*sinfi + (ygrid - y0).*cosfi);
    psf = exp(-0.5 .* ((dx./xsigma).^2 + (dy./ysigma).^2)) ./ (2*pi*xsigma*ysigma);
end