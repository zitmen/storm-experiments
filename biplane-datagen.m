function datagen()

    % start dip
    addpath('C:\Program Files\DIPimage 2.6');
    dipstart();
    % use some of the Gu's methods
    addpath('./cs-Gu2014');
    % use some of the 3denseSTORM methods
    addpath('../../astigmatism/SNR+density/cs-sw');

    % Setup
    header = {'"frame"','"x [nm]"','"y [nm]"','"z [nm]"','"I [photon]"','"background [photon]"','"camera_offset [ADU]"'};
    frames = 100;
    img_width = 32; % [px]
    img_a = 0.08;   % pixelsize [um]
    bkgPhotons = 70;
    cameraOffset = 100;
    m_par = MakePar(0.22,-0.15,0.15,0.4,0,0);    % 3D Calibration [um]

    [img_xx,img_yy,img_zz] = meshgrid((0.5:img_width).*img_a,(0.5:img_width).*img_a,0);
    
    densities = [0.1,0.5:0.5:20.0]; % mol/um2
    for di = 1:length(densities)
        density = densities(di);
        fprintf('========== DENSITY %g ===========\n',density);
        data = sprintf('density=%g',density);
        for intensity=2500%500:1000:2500
            fprintf('========== INTENSITY %d ==========\n',intensity);
            file = sprintf('I=%d+bkg=70',intensity);
            molecules = dlmread(['../../astigmatism/SNR+density/',data,'/',file,'.csv'],',',1,0);
            stack = cell(frames,1);
            table = [];
            for ff = 1:frames
                idx = (molecules(:,1) == ff);
                m_num = sum(idx);
                % Gen coordinates
                m_x = molecules(idx, 2) ./ 1e3; % nm -> um
                m_y = molecules(idx, 3) ./ 1e3; % nm -> um
                m_z = molecules(idx,10) ./ 1e3; % nm -> um
                
                % Gen image
                img = zeros(img_width,img_width,2);
                for m=1:m_num
                    timg = GenPSFs(img_xx,img_yy,img_zz,m_x(m),m_y(m),m_z(m),m_par);
                    img = img + timg;
                end

                img = img.*(intensity/2);
                I=img+bkgPhotons;
                I=noise(I,'poisson',1.0);
                I=dip_array(I);
                I=I+cameraOffset;
                stack{ff} = I;
                
                gt_data = [ones(length(m_x),1)*ff,m_x(:)*1000,m_y(:)*1000,m_z(:)*1000, ...
                           ones(length(m_x),1)*intensity,ones(length(m_x),1)*cameraOffset, ...
                           ones(length(m_x),1)*bkgPhotons];
                table = [table; gt_data];
            end
            
            % Save the image and the ground-truth data
            if ~exist(data,'dir'), mkdir(data); end;
            IO.writeResults(sprintf('%s/%s.csv',data,file),header,table,0);
            for ff=1:frames
                im = stack{ff};
                img = [im(:,:,1),im(:,:,2)];
                if ff == 1
                    imwrite(uint16(img),sprintf('%s/%s.tif',data,file),'Compression','none','WriteMode','overwrite');
                else
                    imwrite(uint16(img),sprintf('%s/%s.tif',data,file),'Compression','none','WriteMode','append');
                end
            end
        end
    end

end

function imgs = GenPSFs(xx,yy,zz,x0,y0,z0,par)
    wx1 = GenWidthofPSF(z0(1)-zz(1),par.w0,par.cx1,par.d,par.ax,par.bx);
    wx2 = GenWidthofPSF(z0(1)-zz(1),par.w0,par.cx2,par.d,par.ax,par.bx);
    s=size(xx);
    imgs = zeros([s 2]);
    G1 = GenDoubleGaussianPeak(xx,yy,x0,y0,wx1);
    G2 = GenDoubleGaussianPeak(xx,yy,x0,y0,wx2);
    imgs(:,:,1) = G1 ./ sum(G1(:));
    imgs(:,:,2) = G2 ./ sum(G2(:));
end
