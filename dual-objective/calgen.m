function calgen()

    % Settings
    imsize = 64;      % [px]
    pxsize = 80;      % [nm]
    intensity = 30000; % [photons]
    background = 70;  % [photons]
    AD = 1;           % A/D conversion rate
    offset = 100;     % [digital]
    zrange = -400:10:+400;     % number of frames
    
    nmols = 4;
    pos(1,:) = [10,10,+256];
    pos(2,:) = [50,21,-135];
    pos(3,:) = [19,30,+50];
    pos(4,:) = [30,40,+90];
    
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
    
    data = 'cal1';
    calibration.cx = [-150;+150];  % [nm]
    calibration.cy = [-150;+150];  % [nm]
    
    %data = 'cal2';
    %calibration.cx = [-150;+150];  % [nm]
    %calibration.cy = [+150;-150];  % [nm]
    
    %data = 'cal3';
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
    file = sprintf('I=%d+bkg=%d',intensity,background);
    if ~exist(data,'dir'), mkdir(data); end;
    header = {'"frame"','"x [nm]"','"y [nm]"','"z [nm]"','"I [photon]"','"background [photon]"'};
    
    frames = length(zrange);
    table = zeros(frames*nmols,6);
    h_wb = waitbar(0,'Generating ...');
    for ff=1:frames
        
        waitbar(ff/frames,h_wb,sprintf('Generating frame %d / %d ...',ff,frames));

        % Positions of molecules
        mols_x = pxsize .* pos(:,1);
        mols_y = pxsize .* pos(:,2);
        mols_z = zrange(ff) - pos(:,3);

        % Render the image (astigmatism)
        img = zeros(imsize,2*imsize);
        for m=1:nmols
            timg = GenPSFs(calibration,xgrid,ygrid,mols_x(m)./pxsize,mols_y(m)./pxsize,mols_z(m));
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