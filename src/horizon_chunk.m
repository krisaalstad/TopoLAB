function [hout,svfout]=horizon_chunk(xt,yt,x,y,z)
[X,Y]=meshgrid(x,y);
del=x(2)-x(1);
% bins
bins=0:10:350; nbins=numel(bins);

Np=numel(xt);
hout=zeros(Np,nbins);
svfout=zeros(Np,1);
for j=1:Np
    xis=xt(j);
    yis=yt(j);
    
    d=sqrt((X-xis).^2+(Y-yis).^2);
    % Calculate bearings to other grid cell centroids and corners.
    shiftx=[0; -1; -1; 1; 1]; % center, ulc, blc, brc, urc
    shifty=[0; -1; 1; 1; -1]; % -"-
    for k=1:numel(shiftx)
        dX=(X+shiftx(k).*del/2)-xis;
        dY=(Y+shifty(k).*del/2)-yis;
        if k==1
            bearing=zeros(size(dX,1),size(dX,2),numel(shiftx));
        end
        tmp=atan2(dX,dY);
        bearing(:,:,k)=rad2deg(tmp+2.*pi.*(tmp<0));
    end
    
    dz=z-z(d==0); % Change in elevation.
    eang=atand(dz./d); % Elevation angle.
    for nb=1:nbins
        if nb~=nbins
            blim=[bins(nb) bins(nb+1)];
        else
            blim=[bins(nb) 360];
        end
        % Check all staggered grids to also look for grid cell
        % corners that fall inside the bearing range.
        inbin=bearing>=blim(1)&bearing<blim(2);
        inbin=any(inbin,3);
        if any(inbin(:))
            thesea=eang(inbin);
            maxa=inbin&eang==max(thesea);
            thesed=d(maxa);
            here=maxa&d==min(thesed);
            his=uint8(eang(here).*(eang(here)>0)); % 0^deg as lower bound on elevation angle. Relax later.
        else
            his=0;
        end
        hout(j,nb)=mean(his);
        svfout(j)=1-mean(sind(double(hout(j,:))));
    end
  
end

end