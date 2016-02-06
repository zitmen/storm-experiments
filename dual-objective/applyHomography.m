function loc = applyHomography(H, P, imsize)
    P = moveToCenterXY(P, imsize);
    P = applyH(H,P')';
    loc = moveToBoundaryXY(P, imsize);
end

function pout = applyH(H,pin)
    n = size(pin,2);
    pout = H*[pin;ones(1,n)];
    pout = pout(1:2,:)./repmat(pout(3,:),2,1);
end

function P = moveToCenterXY(pts,imsize)
    shift = imsize(1:2) ./ 2;
    P = pts - repmat(shift,[size(pts,1),1]);
end

function P = moveToBoundaryXY(pts,imsize)
    shift = imsize(1:2) ./ 2;
    P = pts + repmat(shift,[size(pts,1),1]);
end
