function datagen()

    % Settings
    imsize = 256;     % [px]
    border = 6;       % [px]
    pxsize = 80;      % [nm]
    zrange = 400;     % [nm]
    intensity = 2500; % [photons]
    background = 70;  % [photons]
    AD = 1;           % A/D conversion rate
    offset = 100;     % [digital]
    frames = 100;     % number of frames
    nmols = 10;        % molecules per frame
    
    % Homography: transformation from the first image to the second one
    tx = 0; ty = 0; sx = 1; sy = 1; fx = +1; fy = -1; a = deg2rad(5); shx = 0.3; shy = 0.8; % sxy=scale;txy=translate;fxy=flip;a(lpha)=rot.angle;shxy=shear
    Hctl = [fx*sx*cos(a),shx*sin(a),tx;-shy*sin(a),fy*sy*cos(a),ty;0,0,1];  % manualy approximated "control"!
    
    calibration.width = imsize;
    calibration.height = imsize;
    calibration.homography = Hctl;
    calibration.is_biplane = true;
	calibration.divide_dim = 2;
    calibration.px = pxsize;
    calibration.w0 =    218.4;  % [nm]
    calibration.d  =  400.00;  % [nm]
    calibration.fi =    0.00;  % [rad]
    
    caldata = 'cal1';
    calibration.cx = [-150;+150];  % [nm]
    calibration.cy = [-150;+150];  % [nm]
    
    %caldata = 'cal2';
    %calibration.cx = [-150;+150];  % [nm]
    %calibration.cy = [+150;-150];  % [nm]
    
    %caldata = 'cal3';
    %calibration.cx = [-150;+150];  % [nm]
    %calibration.cy = [+150;+150];  % [nm]
    
    % Prepare grid for rendering
    [xgrid,ygrid] = meshgrid(0.5:imsize,0.5:imsize);
    xgrid = repmat(xgrid,[1,1,2]); ygrid = repmat(ygrid,[1,1,2]);
    xg = reshape(xgrid(:,:,2),[size(xgrid,1)*size(xgrid,2),1]);
    yg = reshape(ygrid(:,:,2),[size(ygrid,1)*size(ygrid,2),1]);
    tg = applyHomography(inv(calibration.homography),[xg,yg],[size(xgrid,1),size(xgrid,2)]);
    xg = tg(:,1);
    yg = tg(:,2);
    xgrid(:,:,2) = reshape(xg,[size(xgrid,1),size(xgrid,2)]);
    ygrid(:,:,2) = reshape(yg,[size(ygrid,1),size(ygrid,2)]);
    
    % Init & run
    data = sprintf('nmols=%d-%s',nmols,caldata);
    file = sprintf('I=%d+bkg=%d',intensity,background);
    if ~exist(data,'dir'), mkdir(data); end;
    header = {'"frame"','"x [nm]"','"y [nm]"','"z [nm]"','"I [photon]"','"background [photon]"'};
    
    table = zeros(frames*nmols,6);
    h_wb = waitbar(0,'Generating ...');
    for ff=1:frames
        
        waitbar(ff/frames,h_wb,sprintf('Generating frame %d / %d ...',ff,frames));

        % Generate positions of molecules
        rawsize = imsize-2*border;
        pos(:,1:2) = rawsize*rand(nmols,2);
        pos(:,3) = 2*rand(nmols,1) - 1;
        
        mols_x = pxsize .* (pos(:,1) + border);
        mols_y = pxsize .* (pos(:,2) + border);
        mols_z = zrange .*  pos(:,3);

        % Render the image (astigmatism)
        img = zeros(imsize,2*imsize);
        for m=1:nmols
            timg = GenPSFs(calibration,xgrid,ygrid,mols_x(m)/pxsize,mols_y(m)/pxsize,mols_z(m));
            img = img + timg;
        end

        I = img .* (intensity/2) + background;
        I = imnoise(uint16(I),'poisson');
        I = I .* AD + offset;
		
		imagesc(I),colorbar,drawnow;
        
        fmode = 'append'; if ff==1, fmode = 'overwrite'; end;
        imwrite(uint16(I),sprintf('%s/%s.tif',data,file),'Compression','none','WriteMode',fmode);
        
        gt_data = [ones(nmols,1)*ff,mols_x(:),mols_y(:),mols_z(:), ...
                   ones(nmols,1)*intensity,ones(nmols,1)*background];
        table((ff-1)*nmols+1:ff*nmols,:) = gt_data;
        
    end
    writeResults(sprintf('%s/%s.csv',data,file),header,table,0);
    close(h_wb);
    
end

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

function writeResults(fpath,header,data,last_writeout)
    if last_writeout == 0
        hstr = header{1};
        for i = 2:length(header), hstr = [hstr,',',header{i}]; end;
        fid = fopen(fpath,'w');
        fprintf(fid,'%s\r\n',hstr);
        fclose(fid);
    end
    dlmwrite(fpath,data(last_writeout+1:end,:),'-append','delimiter',',');
end
