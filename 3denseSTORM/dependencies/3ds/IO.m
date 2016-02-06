classdef IO
    
    methods (Static)
        
        function [info,frames,width,height] = getImageInfo(filepath,cal)
            info = imfinfo(filepath);
            frames = 1:numel(info);
            width  = info(1).Width;
            height = info(1).Height;
            if cal.is_biplane
                if cal.divide_dim == 1
                    height = height / 2;
                else
                    width = width / 2;
                end
            end
        end
        
        function I = readImage(cfg,frame)
            im = imread(cfg.file,frame,'Info',cfg.fileinfo);
            if cfg.calibration.is_biplane
                I(:,:,1) = im(1:cfg.height,1:cfg.width);
                if cfg.calibration.divide_dim == 1
                    I(:,:,2) = im(cfg.height+1:end,1:cfg.width);
                else
                    I(:,:,2) = im(1:cfg.height,cfg.width+1:end);
                end
            else
                I = im;
            end
        end
        
        function table = getResults(frame,cfg,data)
            f = ones(size(data.X)) .* frame;
            x = data.X ./ cfg.zoom .* cfg.camera.pixelsize;   % px --> nm
            y = data.Y ./ cfg.zoom .* cfg.camera.pixelsize;   % px --> nm
            z = data.Z .* 1e3; % um --> nm
            [sx,sy] = cfg.calibration.sigma(data.Z);   % z [nm] --> sigma_{x,y} [px]
            sx = sx .* cfg.camera.pixelsize; sy = sy .* cfg.camera.pixelsize; % px --> nm
            I = data.I; % photons
            off = data.O; % photons
            %
            if cfg.calibration.is_biplane
                table = [f(:),x(:),y(:),z(:),I(:),sx(1,:)',sy(1,:)',sx(2,:)',sy(2,:)',off(:)];
            else
                table = [f(:),x(:),y(:),z(:),I(:),sx(:),sy(:),off(:)];
            end
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
        
    end
    
end
