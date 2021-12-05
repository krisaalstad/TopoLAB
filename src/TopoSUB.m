function [cn,indsc]=TopoSUB(tp,Nc,maxit,doxy)
%% [cn,indc]=TopoSUB(tp,Nc,maxit,doxy)
% Matlab implementation of the TopoSUB clustering approach.
% outputs:
%       cn : Cluster number \in[1,Nc] for each of the Np input grid cells.
%       indc : Index\in[1,Np] of sample centroids, i.e. the points that are
%               closest to the respective cluster centroids.
Np=numel(tp.svf);
X=[tp.svf -cos(tp.asp) sin(tp.asp) sin(tp.slp) tp.z];
if doxy
    X=[X tp.x tp.y]; % Optional, but important to get the right synoptic weather situation for large areas
end
Zs=zscore(X,1,1); % Does this make sense for bounded variables, or should you transform first?

% Note, may require many iterations (more than standard 1e2) for large Nc
[cn,~,~,D]=kmeans(Zs,Nc,'MaxIter',maxit); % 

% Find points closest to centroids within each cluster, henceforth these
% points are the sample centroids.
indsc=zeros(Nc,1);% Index (possible range: 1 to total number of points) of 
ind=(1:Np)';
for j=1:Nc
    these=(cn==j); % Points in cluster j
    indc=ind(these); % Indices of points in cluster j
    Dc=D(these); % Distances to centroid of points in cluster j.
    here=(Dc==min(Dc)); % Minimum distance to centroid.
    indsc(j)=min(indc(here)); % "min()" is arbitrary and in case of multiple minima.
end

end
