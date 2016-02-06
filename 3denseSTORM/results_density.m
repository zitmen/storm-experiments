function results_density()
    global bkg pixelsize regionsize nframes frames_per_measurement;

    bkg = 70.0;
    pixelsize = 0.08; % [um]
    regionsize = 20; % [px]
    nframes = 100;
    frames_per_measurement = 20;
    
    max_density = 18.6; % [mols/um2]
    tol = 200; % xy [nm]
    snrI = 2500; % [photons]
    densities = 0.5:0.5:20.0;%[0.1,0.5:0.5:20.0]; % [molecules/um^2]
    intensities = 2500;%500:1000:2500; % [photons]
    
    %% Results from the summary
    
    daostorm = loadResults('astigmatism','daostorm',tol,densities,intensities);
    cs_astigmatism = loadResults('astigmatism','CS-MLE',tol,densities,intensities);
    cs_biplane = loadResults('biplane','CS-MLE',tol,densities,intensities);
    sa = loadResults('astigmatism','SA-WLSQ',tol,densities,intensities);
    
    daostorm.recovered_density = daostorm.TP ./ ((regionsize*pixelsize)^2) ./ nframes;
    cs_astigmatism.recovered_density = cs_astigmatism.TP ./ ((regionsize*pixelsize)^2) ./ nframes;
    cs_biplane.recovered_density = cs_biplane.TP ./ ((regionsize*pixelsize)^2) ./ nframes;
    sa.recovered_density = sa.TP ./ ((regionsize*pixelsize)^2) ./ nframes;
    
    %
    
    figure(1); clf;

    subplot(2,2,1);
    hold on;
    title('F1');
    idx = (sa.I == snrI);
    plot(daostorm.density(idx),daostorm.('F1')(idx),'-*r');
    plot(cs_astigmatism.density(idx),cs_astigmatism.('F1')(idx),'-sb');
    plot(cs_biplane.density(idx),cs_biplane.('F1')(idx),'-sg');
    plot(sa.density(idx),sa.('F1')(idx),'-om');
    xlabel('density [\mum^{-2}]'),ylabel('F1-score [-]');
    legend('3D DAOSTORM','3denseSTORM (astigmatism)','3denseSTORM (biplane)','3D single molecule fitting','Location','Best');
    hold off; 
    
    subplot(2,2,2);
    hold on;
    title('recovered density');
    idx = (sa.I == snrI);
    plot(densities,densities,':k');
    plot(densities,max_density*ones(size(densities)),':r');
    plot(daostorm.density(idx),daostorm.('recovered_density')(idx),'-*r');
    plot(cs_astigmatism.density(idx),cs_astigmatism.('recovered_density')(idx),'-sb');
    plot(cs_biplane.density(idx),cs_biplane.('recovered_density')(idx),'-sg');
    plot(sa.density(idx),sa.('recovered_density')(idx),'-om');
    xlabel('density [\mum^{-2}]'),ylabel('Recovered density [\mum^{-2}]');
	legend('Image density','Theoretical limit','3D DAOSTORM','3denseSTORM (astigmatism)','3denseSTORM (biplane)','3D single molecule fitting','Location','Best');
    hold off; 
    
    subplot(2,2,3);
    hold on;
    title('dxy');
    idx = (sa.I == snrI);
    plot(daostorm.density(idx),daostorm.('dxy')(idx),'-*r');
    plot(cs_astigmatism.density(idx),cs_astigmatism.('dxy')(idx),'-sb');
    plot(cs_biplane.density(idx),cs_biplane.('dxy')(idx),'-sg');
    plot(sa.density(idx),sa.('dxy')(idx),'-om');
    xlabel('density [\mum^{-2}]'),ylabel('Localization error \Delta_{xy} [nm]');
    legend('3D DAOSTORM','3denseSTORM (astigmatism)','3denseSTORM (biplane)','3D single molecule fitting','Location','Best');
    hold off; 
    
    subplot(2,2,4);
    hold on;
    title('dz');
    idx = (sa.I == snrI);
    plot(daostorm.density(idx),daostorm.('dz')(idx),'-*r');
    plot(cs_astigmatism.density(idx),cs_astigmatism.('dz')(idx),'-sb');
    plot(cs_biplane.density(idx),cs_biplane.('dz')(idx),'-sg');
    plot(sa.density(idx),sa.('dz')(idx),'-om');
    xlabel('density [\mum^{-2}]'),ylabel('Localization error \Delta_{z} [nm]');
    legend('3D DAOSTORM','3denseSTORM (astigmatism)','3denseSTORM (biplane)','3D single molecule fitting','Location','Best');
    hold off;
    
    %% Localization errors for each individual molecule
    
    daostorm = loadLocErrs('astigmatism','daostorm',tol,densities,intensities);
    cs_astigmatism = loadLocErrs('astigmatism','CS-MLE',tol,densities,intensities);
    cs_biplane = loadLocErrs('biplane','CS-MLE',tol,densities,intensities);
    sa = loadLocErrs('astigmatism','SA-WLSQ',tol,densities,intensities);
    
    global XCOORD;
    XCOORD = daostorm.density(idx);
    
    figure(2); clf;
    
    figure(2);clf;%subplot(2,2,1);
    title('F1');
    idx = (sa.I == snrI);
    MU(:,1) = daostorm.('F1_mu')(idx);
    MU(:,2) = cs_astigmatism.('F1_mu')(idx);
    MU(:,3) = cs_biplane.('F1_mu')(idx);
    MU(:,4) = sa.('F1_mu')(idx);
    STD(:,1) = daostorm.('F1_std')(idx);
    STD(:,2) = cs_astigmatism.('F1_std')(idx);
    STD(:,3) = cs_biplane.('F1_std')(idx);
    STD(:,4) = sa.('F1_std')(idx);
    myeb(MU,STD);
    xlim([0,max(densities)]); ylim([0,1]);
    xlabel('density [\mum^{-2}]'),ylabel('F1-score [-]');
    legend('3D DAOSTORM','3denseSTORM (astigmatism)','3denseSTORM (biplane)','3D single molecule fitting','Location','Best');
    
    figure(3);clf;%subplot(2,2,2);
    title('recovered density');
    idx = (sa.I == snrI);
    MU(:,1) = daostorm.('recovered_density_mu')(idx);
    MU(:,2) = cs_astigmatism.('recovered_density_mu')(idx);
    MU(:,3) = cs_biplane.('recovered_density_mu')(idx);
    MU(:,4) = sa.('recovered_density_mu')(idx);
    STD(:,1) = daostorm.('recovered_density_std')(idx);
    STD(:,2) = cs_astigmatism.('recovered_density_std')(idx);
    STD(:,3) = cs_biplane.('recovered_density_std')(idx);
    STD(:,4) = sa.('recovered_density_std')(idx);
    myeb(MU,STD);
    hold on;
    plot(densities,densities,':k');
    plot(densities,max_density*ones(size(densities)),':r');
    hold off;
    xlim([0,max(densities)]); ylim([0,max(densities)]);
    xlabel('density [\mum^{-2}]'),ylabel('Recovered density [\mum^{-2}]');
	legend('','','','','3D DAOSTORM','3denseSTORM (astigmatism)','3denseSTORM (biplane)','3D single molecule fitting','Image density','Theoretical limit','Location','Best');
    
    figure(4);clf;%subplot(2,2,3);
    title('dxy');
    idx = (sa.I == snrI);
    MU(:,1) = daostorm.('dxy_mu')(idx);
    MU(:,2) = cs_astigmatism.('dxy_mu')(idx);
    MU(:,3) = cs_biplane.('dxy_mu')(idx);
    MU(:,4) = sa.('dxy_mu')(idx);
    STD(:,1) = daostorm.('dxy_std')(idx);
    STD(:,2) = cs_astigmatism.('dxy_std')(idx);
    STD(:,3) = cs_biplane.('dxy_std')(idx);
    STD(:,4) = sa.('dxy_std')(idx);
    myeb(MU,STD);
    xlim([0,max(densities)]); ylim([0,100]);
    xlabel('density [\mum^{-2}]'),ylabel('Root mean square error \Delta_{xy} [nm]');
    legend('3D DAOSTORM','3denseSTORM (astigmatism)','3denseSTORM (biplane)','3D single molecule fitting','Location','Best');
    
    figure(5);clf;%subplot(2,2,4);
    title('dz');
    idx = (sa.I == snrI);
    MU(:,1) = daostorm.('dz_mu')(idx);
    MU(:,2) = cs_astigmatism.('dz_mu')(idx);
    MU(:,3) = cs_biplane.('dz_mu')(idx);
    MU(:,4) = sa.('dz_mu')(idx);
    STD(:,1) = daostorm.('dz_std')(idx);
    STD(:,2) = cs_astigmatism.('dz_std')(idx);
    STD(:,3) = cs_biplane.('dz_std')(idx);
    STD(:,4) = sa.('dz_std')(idx);
    myeb(MU,STD);
    xlim([0,max(densities)]); ylim([0,300]);
    xlabel('density [\mum^{-2}]'),ylabel('Root mean square error \Delta_{z} [nm]');
    legend('3D DAOSTORM','3denseSTORM (astigmatism)','3denseSTORM (biplane)','3D single molecule fitting','Location','Best');
        
end

%%

function R = loadResults(experiment,method,tol,densities,intensities)
%header = {'Distance radius [nm]','# of TP','# of FP','# of FN','Jaccard index',
%'precision','recall','F1-measure','RMSE x [nm]','RMSE y [nm]','RMSE lateral [nm]',
%'RMSE axial [nm]','RMSE total [nm]'};
    global bkg;
    n = 0;
    for density = densities
        for intensity = intensities
            n = n + 1;
            fpath = sprintf('%s/SNR+density/density=%g/results/%s+I=%d+bkg=%d.xls',...
                experiment,density,method,intensity,bkg);
            M = dlmread(fpath,'\t',1,1);
            idx = (M(:,1) == tol);
            R.I(n) = intensity;
            R.density(n) = density;
            R.TP(n) = M(idx,2);
            R.FP(n) = M(idx,3);
            R.FN(n) = M(idx,4);
            R.Jaccard(n) = M(idx,5);
            R.precision(n) = M(idx,6);
            R.recall(n) = M(idx,7);
            R.F1(n) = M(idx,8);
            R.dx(n) = M(idx,9);
            R.dy(n) = M(idx,10);
            R.dxy(n) = M(idx,11);
            R.dz(n) = M(idx,12);
            R.dxyz(n) = M(idx,13);
        end
    end
end

function R = loadLocErrs(experiment,method,tol,densities,intensities)
    global bkg regionsize pixelsize nframes frames_per_measurement;
    frames = 0:frames_per_measurement:nframes;
    n = 0;
    for density = densities
        for intensity = intensities
            n = n + 1;
            fpath = sprintf('%s/SNR+density/density=%g/results/%s+I=%d+bkg=%d_tol=%d.csv',...
                experiment,density,method,intensity,bkg,tol);
            fpathGT = sprintf('%s/SNR+density/density=%g/I=%d+bkg=%d.csv',experiment,density,intensity,bkg);
            [M,F] = dataread(fpath);
            G = datareadGT(fpathGT);
            R.I(n) = intensity;
            R.density(n) = density;
            dxy = []; dz = []; dxyz = [];
            TP = []; FP = []; FN = [];
            for ff = 1:length(frames)-1
                idx = M(:,1) > frames(ff) & M(:,1) <= frames(ff+1);
                dxy(ff) = mean(M(idx,2));
                dz(ff) = mean(M(idx,3));
                dxyz(ff) = mean(M(idx,4));
                TP(ff) = sum(idx);
                FP(ff) = sum(F > frames(ff) & F <= frames(ff+1));
                FN(ff) = sum(G > frames(ff) & G <= frames(ff+1)) - TP(ff);
            end
            precision = TP ./ (TP + FP);
            recall = TP ./ (TP + FN);
            F1 = 2.*precision.*recall./(precision+recall);
            recovered_density = TP ./ ((regionsize*pixelsize)^2) ./ frames_per_measurement;
            R.recovered_density_mu(n) = mean(recovered_density);
            R.recovered_density_std(n) = std(recovered_density);
            R.F1_mu(n) = mean(F1);
            R.F1_std(n) = std(F1);
            R.dxy_mu(n) = mean(dxy);
            R.dxy_std(n) = std(dxy);
            R.dz_mu(n) = mean(dz);
            R.dz_std(n) = std(dz);
            R.dxyz_mu(n) = mean(dxyz);
            R.dxyz_std(n) = std(dxyz);
        end
    end
end

function [D,FP] = dataread(fpath)
% D = [frame,error_xy,error_z,error_xyz]
    FP = []; D = [];
    items = 0; fpc = 0;
    fid = fopen(fpath,'r');
    fgetl(fid);  % skip over the header
    while ~feof(fid)
        line = fgetl(fid);
        %"frame","gt_id","gt_dist_xy [nm]","gt_dist_z [nm]","gt_dist_xyz [nm]"
        d = textscan(line,'%f%f%s%s%s','delimiter',',','whitespace',' ');
        if d{2} ~= 0    % gt_id != 0
            items = items + 1;
            D(items,:) = [d{1},str2num(d{3}{1}),str2num(d{4}{1}),str2num(d{5}{1})];
        else
            fpc = fpc + 1;
            FP(fpc) = d{1};
        end
    end
    fclose(fid);
end

function D = datareadGT(fpath)
% D = [frame]
    D = [];
    items = 0;
    fid = fopen(fpath,'r');
    fgetl(fid);  % skip over the header
    while ~feof(fid)
        line = fgetl(fid);
        %"frame"
        d = textscan(line,'%f','delimiter',',','whitespace',' ');
        items = items + 1;
        D(items) = d{1}(1);
    end
    fclose(fid);
end


%% 3rd party library for errorbar ploting

function myeb(Y,varargin)
%
% myeb(Y,varargin);
%
% This function makes nice coloured, shaded error bars. Exactly what
% it does depends on Y, and on whether you give it one or two inputs. 
%
% If you only pass it Y, and no other arguments, it assuemd you're
% giving it raw data. 
%
%		myeb(Raw_Data)
%
% 	.) if Y is 2D array, it will then plot mean(Y) with errorbars given
% 	by std(Y). In this case there is only one mean vector with its
% 	errorbars. 
% 
%	.) if Y is 3D array, it will plot size(Y,3) lines with the
%	associated errorbars. Line k will be mean(Y(:,:,k)) with errorbars
%	given by std(Y(:,:,k))
%
% If you pass it 2 arguments, each has to be at most 2D. 
%
%		myeb(mu,std)
%
% 	.) if mu and std are 1D, it just plots one line given by mu with a
% 	shaded region given by std. 
%
%	.) if mu and std are 2D, it will plot size(Y,2) lines in the
%	standard sequence of colours; each line mu(:,k) will have a shaded
%	region in the same colour, but less saturated given by std(:,k)
%
%
% Quentin Huys, 2007
% Center for Theoretical Neuroscience, Columbia University
% Email: qhuys [at] n e u r o theory [dot] columbia.edu
% (just get rid of the spaces, replace [at] with @ and [dot] with .)


global XCOORD;

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;


if length(varargin)==0;

	if length(size(Y))==2 
		m=mean(Y);
		s=std(Y);
		ind1=1:length(m);
		ind2=ind1(end:-1:1);
		hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],.6*ones(1,3));
		set(h,'edgecolor',.6*ones(1,3))
		%plot(ind1,m,'linewidth',2)
        patchline(ind1,m,'linewidth',2,'edgealpha',0.5);
		hold off
	elseif length(size(Y))>2 
		cla; hold on; 
		ind1=1:size(Y,2);
		ind2=ind1(end:-1:1);
		if size(Y,3)>8; col=jet(size(Y,3));ccol=col+.8; ccol(ccol>1)=1;end
		for k=1:size(Y,3)
			m=mean(Y(:,:,k));
			s=std(Y(:,:,k));
			h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],ccol(k,:));
			set(h,'edgecolor',ccol(k,:))
            alpha(h,0.5)
		end
		for k=1:size(Y,3)
			m=mean(Y(:,:,k));
			s=std(Y(:,:,k));
			%plot(ind1,m,'linewidth',2,'color',col(k,:))
            patchline(ind1,m,'edgecolor',col(k,:),'linewidth',2,'edgealpha',0.5);
        end
		hold off 
	end

elseif length(varargin)==1;

	m=Y;
	s=varargin{1};
	if length(size(Y))>2; error;
	elseif min(size(Y))==1;
		if size(m,1)>1; m=m';s=s';end
		ind1=1:length(m);
		ind2=ind1(end:-1:1);
		hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],.6*ones(1,3));
		set(h,'edgecolor',.6*ones(1,3))
		%plot(ind1,m,'linewidth',2)
        patchline(ind1,m,'linewidth',2,'edgealpha',0.5);
		hold off
	else 
		ind1=(1:size(Y,1));
		ind2=ind1(end:-1:1);
		cla; hold on; 
		if size(Y,2)>8; col=jet(size(Y,2));ccol=col+.8; ccol(ccol>1)=1;end
		for k=1:size(Y,2)
			mm=m(:,k)';
			ss=s(:,k)';
            %h=fill([ind1 ind2],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
			h=fill([XCOORD fliplr(XCOORD)],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
            alpha(h,0.5);
			set(h,'edgecolor',ccol(k,:))
        end
		for k=1:size(Y,2);
			mm=m(:,k)';
			ss=s(:,k)';
			%plot(ind1,mm,'linewidth',1.2,'color',col(k,:))
            plot(XCOORD,mm,'linewidth',1.2,'color',col(k,:))
		end
		hold off 
	end
end
end

function p = patchline(xs,ys,varargin)
% Plot lines as patches (efficiently)
%
% SYNTAX:
%     patchline(xs,ys)
%     patchline(xs,ys,zs,...)
%     patchline(xs,ys,zs,'PropertyName',propertyvalue,...)
%     p = patchline(...)
%
% PROPERTIES: 
%     Accepts all parameter-values accepted by PATCH.
% 
% DESCRIPTION:
%     p = patchline(xs,ys,zs,'PropertyName',propertyvalue,...)
%         Takes a vector of x-values (xs) and a same-sized
%         vector of y-values (ys). z-values (zs) are
%         supported, but optional; if specified, zs must
%         occupy the third input position. Takes all P-V
%         pairs supported by PATCH. Returns in p the handle
%         to the resulting patch object.
%         
% NOTES:
%     Note that we are drawing 0-thickness patches here,
%     represented only by their edges. FACE PROPERTIES WILL
%     NOT NOTICEABLY AFFECT THESE OBJECTS! (Modify the
%     properties of the edges instead.)
%
%     LINUX (UNIX) USERS: One test-user found that this code
%     worked well on his Windows machine, but crashed his
%     Linux box. We traced the problem to an openGL issue;
%     the problem can be fixed by calling 'opengl software'
%     in your <http://www.mathworks.com/help/techdoc/ref/startup.html startup.m>.
%     (That command is valid at startup, but not at runtime,
%     on a unix machine.)
%
% EXAMPLES:
%%% Example 1:
%
% n = 10;
% xs = rand(n,1);
% ys = rand(n,1);
% zs = rand(n,1)*3;
% plot3(xs,ys,zs,'r.')
% xlabel('x');ylabel('y');zlabel('z');
% p  = patchline(xs,ys,zs,'linestyle','--','edgecolor','g',...
%     'linewidth',3,'edgealpha',0.2);
%
%%% Example 2: (Note "hold on" not necessary here!)
%
% t = 0:pi/64:4*pi;
% p(1) = patchline(t,sin(t),'edgecolor','b','linewidth',2,'edgealpha',0.5);
% p(2) = patchline(t,cos(t),'edgecolor','r','linewidth',2,'edgealpha',0.5);
% l = legend('sine(t)','cosine(t)');
% tmp = sort(findobj(l,'type','patch'));
% for ii = 1:numel(tmp)
%     set(tmp(ii),'facecolor',get(p(ii),'edgecolor'),'facealpha',get(p(ii),'edgealpha'),'edgecolor','none')
% end
%
%%% Example 3 (requires Image Processing Toolbox):
%%%   (NOTE that this is NOT the same as showing a transparent image on 
%%%         of the existing image. (That functionality is
%%%         available using showMaskAsOverlay or imoverlay).
%%%         Instead, patchline plots transparent lines over
%%%         the image.)
%
% img = imread('rice.png');
% imshow(img)
% img = imtophat(img,strel('disk',15));
% grains = im2bw(img,graythresh(img));
% grains = bwareaopen(grains,10);
% edges = edge(grains,'canny');
% boundaries = bwboundaries(edges,'noholes');
% cmap = jet(numel(boundaries));
% ind = randperm(numel(boundaries));
% for ii = 1:numel(boundaries)
% patchline(boundaries{ii}(:,2),boundaries{ii}(:,1),...
%     'edgealpha',0.2,'edgecolor',cmap(ind(ii),:),'linewidth',3);
% end
%
% Written by Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 5/31/2012
% 
% Revisions:
% 6/26 Improved rice.png example, modified FEX image.
%
% Copyright 2012 MathWorks, Inc.
%
% See also: patch, line, plot

[zs,PVs] = parseInputs(varargin{:});
if rem(numel(PVs),2) ~= 0
    % Odd number of inputs!
    error('patchline: Parameter-Values must be entered in valid pairs')
end

% Facecolor = 'k' is (essentially) ignored here, but syntactically necessary
if isempty(zs)
    p = patch([xs(:);NaN],[ys(:);NaN],'k');
else
    p = patch([xs(:);NaN],[ys(:);NaN],[zs(:);NaN],'k');
end

% Apply PV pairs
for ii = 1:2:numel(PVs)
    set(p,PVs{ii},PVs{ii+1})
end
if nargout == 0
    clear p
end
end

function [zs,PVs] = parseInputs(varargin)
if isnumeric(varargin{1})
    zs = varargin{1};
    PVs = varargin(2:end);
else
    PVs = varargin;
    zs = [];
end
end