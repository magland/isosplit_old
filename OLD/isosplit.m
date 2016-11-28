function [labels,info]=isosplit(X,opts)
% isosplit - perform clustering based on isotonic regression (jfm, may 2105)
%
% labels = isosplit(X,opts) 
%   X is M x N, M=#dimensions, N=#samples
%   labels: 1xN vector of labels from 1..L
%   opts.split_threshold = a number determining how likely clusters are to
%   split in two. The lower the number, the more likely the split.
%
%   opts.K : The starting number of k-means clusters (to save time).
%            If K is large enough this should not affect the result.
%   opts.minsize
%   opts.verbose, opts.verbose2
%   opts.max_iterations_per_number_clusters
%
% Magland 5/19/2015

if (nargin<1)
	%isosplit_test; %run the test code if there are no parameters
    isosplit_iteration_demo;
	return;
end;

timer_total=tic;
info.T_find_best_pair=0;
info.T_find_centroids=0;
info.T_attempt_redistribution=0;
info.T_isosplit1d=0;
info.T_projection=0;
info.T_sort=0;

%default options
if (nargin<2) opts=struct(); end;
if (~isfield(opts,'split_threshold')) opts.split_threshold=0.9; end;
if (~isfield(opts,'K')) opts.K=25; end;
if (~isfield(opts,'minsize')) opts.minsize=3; end;
if (~isfield(opts,'verbose')) opts.verbose=0; end;
if (~isfield(opts,'verbose2')) opts.verbose2=0; end;
if (~isfield(opts,'max_iterations_per_number_clusters')) opts.max_iterations_per_number_clusters=5000; end;
if (~isfield(opts,'return_iterations')) opts.return_iterations=0; end;

[M,N]=size(X);

%initialize with k-means (just to save time)
timer_initialization=tic;
labels=local_kmeans_sorber(X,opts.K);
info.T_initialization=toc(timer_initialization);

centroids=compute_centroids(X,labels);
distances=compute_distances(centroids);

%Here is a list of the attempted cluster splits/redistributions -- we don't
%ever want to repeat any of these.
attempted_redistributions=zeros(0,1);

if (opts.return_iterations) info.iterations={}; end;

info.num_iterations=0;
num_iterations_with_same_number_of_clusters=0;
while true
	info.num_iterations=info.num_iterations+1;
	old_labels=labels;
	%find the closest two clusters to check for redistribution/merging
	%we want the centroids to be as close as possible, but we want to
	%exclude the pairs we have tried previously
	timer_find_best_pair=tic;
	[label1,label2]=find_best_pair(X,labels,centroids,distances,attempted_redistributions);
	info.T_find_best_pair=info.T_find_best_pair+toc(timer_find_best_pair);
	if (opts.verbose)
		fprintf('Trying (%d,%d) -- %d clusters\n',label1,label2,max(labels));
	end;
	if (label1==0) break; end; %This means that we've tried everything ... we are done!
	timer_find_centroids=tic;
	inds1=find(labels==label1);
	inds2=find(labels==label2);
	%centroid1=mean(X(:,inds1),2); centroid2=mean(X(:,inds2),2);
	centroid1=centroids(:,label1); centroid2=centroids(:,label2);
	%code=get_code(X,inds1,inds2); %This is the code used to add to the attempted_redistributions list
	%code=sum(centroid1.*centroid2);
	info.T_find_centroids=info.T_find_centroids+toc(timer_find_centroids);
	%Here is the core procedure -- look at two clusters, and redistribute
	%the points using the isotonic regression.
	timer_attempt_redistribution=tic;
	[ii1,ii2,redistributed,inf0]=attempt_to_redistribute_two_clusters(X,inds1,inds2,centroid1,centroid2,opts.split_threshold,opts);
	info.T_projection=info.T_projection+inf0.T_projection;
	info.T_isosplit1d=info.T_isosplit1d+inf0.T_isosplit1d;
	info.T_sort=info.T_sort+inf0.T_sort;
	info.T_attempt_redistribution=info.T_attempt_redistribution+toc(timer_attempt_redistribution);
	if (length(ii2)>0)
		num_iterations_with_same_number_of_clusters=num_iterations_with_same_number_of_clusters+1;
	else
		num_iterations_with_same_number_of_clusters=0;
	end;
	%add the code to the list of attempted redistributions
	if (distances(label1,label2)==inf)
		error('Unexpected error.');
	end;
	attempted_redistributions(end+1)=distances(label1,label2);
	if (redistributed)
		%okay, we've changed something. Now let's update the labels
		if (opts.verbose)
			fprintf('  Redistributed (%d,%d)->(%d,%d)\n',length(inds1),length(inds2),length(ii1),length(ii2));
		end;
		labels(ii1)=label1;
		labels(ii2)=label2;
		centroids(:,label1)=mean(X(:,ii1),2);
		tmp0=sum((centroids-repmat(centroids(:,label1),1,size(centroids,2))).^2,1);
		tmp0(ismember(tmp0,attempted_redistributions))=inf;
		distances(:,label1)=tmp0;
		distances(label1,:)=distances(:,label1); distances(label1,label1)=inf;
		if (length(ii2>0))
			centroids(:,label2)=mean(X(:,ii2),2);
			tmp0=sum((centroids-repmat(centroids(:,label2),1,size(centroids,2))).^2,1);
			tmp0(ismember(tmp0,attempted_redistributions))=inf;
			distances(:,label2)=tmp0;
			distances(label2,:)=distances(:,label2); distances(label2,label2)=inf;
		else
			centroids(:,label2)=0;
		end;
		[labels,centroids,distances]=normalize_labels(labels,centroids,distances); % we may have eliminated a label, so let's shift the labelings down
        if (opts.return_iterations) info.iterations{end+1}=labels; end;
	else
		distances(label1,label2)=inf;
		distances(label2,label1)=inf;
	end;

	% The following might be a bad idea.
  	if (num_iterations_with_same_number_of_clusters>opts.max_iterations_per_number_clusters)
  		warning(sprintf('%d iterations with same number of clusters.... stopping',num_iterations_with_same_number_of_clusters));
  		return;
  	end;
	%toc
end;

end

function [ii1,ii2,redistributed,info0]=attempt_to_redistribute_two_clusters(X,inds1,inds2,centroid1,centroid2,split_threshold,opts)

%make sure we have row vectors
timer_projection=tic;
if (~isrow(inds1)) inds1=inds1'; end;
if (~isrow(inds2)) inds2=inds2'; end;
[M,N]=size(X);
inds12=cat(2,inds1,inds2);
X1=X(:,inds1);
X2=X(:,inds2);

V=centroid2-centroid1; %the vector from one centroid to another
%V=find_svm_discriminant_direction(X1,X2);

if (sum(V.^2)==0)
	warning('isosplit: vector V is null.');
end;
V=V/sqrt(sum(V.^2));
XX=V'*X(:,inds12); %Project onto the line connecting the centroids
info0.T_projection=toc(timer_projection);

opts2.verbose=opts.verbose2;
opts2.minsize=opts.minsize;
if (isfield(opts,'m_max')) opts2.m_max=opts.m_max; end;
timer_sort=tic;
XXs=sort(XX); spacings=XXs(2:end)-XXs(1:end-1);
info0.T_sort=toc(timer_sort);

if (length(find(spacings==0))>0)
	warning('spacings=0.');
end;
if (isnan(XX(1)))
	warning('isosplit: isnan');
end;
timer_isosplit1d=tic;
[labels2,score0]=isosplit1d(XX,opts2); %This is the core procedure -- split based on isotonic regression
info0.T_isosplit1d=toc(timer_isosplit1d);
if (score0>split_threshold)
	%It was a statistically significant split -- so let's redistribute!
	ii1=inds12(find(labels2==1));
	ii2=inds12(find(labels2==2));
else
	ii1=inds12;
	ii2=zeros(0,1);
end;

%We have redistributed, unless everything is still the same.
redistributed=1;
if ((length(ii1)==length(inds1))&&(length(ii2)==length(inds2)))
	if (length(find(sort(ii1)==sort(inds1)))==length(ii1))
		if (length(find(sort(ii2)==sort(inds2)))==length(ii2))
			redistributed=0;
		end;
	end;
end;

end


function [label1,label2]=find_best_pair(X,labels,centroids,distances,attempted_redistributions)

%Compute the distances between the centroids...
%Except don't count the pairs we have already attempted (set those to inf)
L=max(labels);
% dists=ones(L,L)*inf;
% 
% for j=1:L
% 	for k=j+1:L
% 		if (j~=k)
% 			%code=get_code(X,inds1,inds2);
% 			code=sum(centroids(:,j).*centroids(:,k));
% 			if (length(find(attempted_redistributions==code))==0)
% 				dist0=(sum((centroids(:,j)-centroids(:,k)).^2));
% 				%dist0=(sum(centroids(:,j)-centroids(:,k)).^2);
% 				dists(j,k)=dist0;
% 			end;
% 		end;
% 	end;
% end;

dists=distances;

if (min(dists(:))==inf)
	%There are none we haven't attempted... we are done!
	label1=0; label2=0;
	return;
end;

%Find the pair that has the minimum distance between centroids (excluding
%those already attempted)
[value, location] = min(dists(:));
[label1,label2] = ind2sub(size(dists),location);

% found=false;
% while ~found
% 	label1=randi(L);
% 	[val,loc]=min(dists(label1,:));
% 	if (~(val==inf))
% 		label2=loc;
% 		found=true;
% 	end;
% end;



end

% function code=get_code(X,inds1,inds2)
% %This is tricky -- it's a code that represents a particular attempt to
% %merge/redistribute two clusters.
% code=mean(mean(X(:,inds1),1),2)*mean(mean(X(:,inds2),1),2); %Not strictly guaranteed to be unique ... but very very likely.
% end

function centroids=compute_centroids(X,labels)

L=max(labels);
M=size(X,1);
centroids=zeros(M,L);
for j=1:L
	inds=find(labels==j);
	centroids(:,j)=mean(X(:,inds),2);
end;

end

function distances=compute_distances(centroids)

L=size(centroids,2);
distances=repmat(sum(centroids.^2,1),L,1)+repmat(sum(centroids.^2,1),L,1)'-2*centroids'*centroids;
for j=1:L
	distances(j,j)=inf;
end;

end

function [labels2,centroids2,distances2]=normalize_labels(labels,centroids,distances)
N=length(labels);
labels2=ones(1,N);
centroids2=zeros(size(centroids,1),0);

maps=zeros(1,max(labels));
map_inds=[];
lll=1;
for k=1:max(labels)
	inds=find(labels==k);
	if (length(inds)>0)
		maps(k)=lll;
		map_inds(end+1)=k;
		lll=lll+1;
	end;
end;

labels2=maps(labels);
centroids2=centroids(:,map_inds);
distances2=distances(map_inds,map_inds);

end

function [L,C]=local_kmeans_sorber(X,k)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.

%   Version: 2013-02-08
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.

L = [];
L1 = 0;

while length(unique(L)) ~= k
    
    % The k-means++ initialization.
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L);
%        D = cumsum(sqrt(dot(D,D,1)));  % orig, seems to be dist (l=1)
        D = cumsum(dot(D,D,1));  % Arthur-Vassilvitskii use dist^2 (l=2)
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1));
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
    end
    
    % The k-means algorithm.
    while any(L ~= L1)
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
    end
    
end

end

function isosplit_iteration_demo

close all;

seed0=randi(10000);
seed0=6239;
rng(seed0);
centers={[0,0],[5,3.8],[-4.5,3]};
pops={1000,800,400};
shapes={[1,1,0],[2,1,0],[1,2,0]};
opts.split_threshold=0.3;
opts.K=25;
opts.return_iterations=1;

fprintf('seed = %d\n',seed0);

samples=zeros(2,0);
true_labels=zeros(1,0);

for j=1:length(centers)
	xx=randn(1,pops{j});
	yy=randn(1,pops{j});
	shape=shapes{j};
	xx2=xx*shape(1)+yy*shape(3);
	yy2=yy*shape(2)-xx*shape(3);
	center=centers{j};
	xx2=xx2+center(1);
	yy2=yy2+center(2);
	tmp=zeros(2,pops{j});
	tmp(1,:)=xx2; tmp(2,:)=yy2;
	samples=[samples,tmp];
	true_labels=[true_labels,ones(1,pops{j})*j];
end;

colors='rgbkymc';
[labels,info]=isosplit(samples,opts);
fprintf('num clusters = %d\n',max(labels));




%# figure




fprintf('seed = %d\n',seed0);
fA=figure; set(fA,'position',[200,200,700,700]);
figure, set(fA, 'Color','white')
for ii=1:length(info.iterations)
    figure(fA);
    labels0=info.iterations{ii};
    for j=1:max(labels0)
        xx=samples(1,labels0==j);
        yy=samples(2,labels0==j);
        col=colors(mod(j-1,length(colors))+1);
        plot(xx,yy,['.',col],'markersize',8);
        if (j==1) hold on; end;
    end
    f = getframe(gca); f0=f;
    if (ii==1)
        set(gca, 'nextplot','replacechildren', 'Visible','off');
        [f,map] = rgb2ind(f.cdata, 256, 'nodither');
        mov = repmat(f, [1 1 1 length(info.iterations)]);
    end;
    % seems to be a problem
    %mov(:,:,1,ii) = rgb2ind(f0.cdata, map, 'nodither');
    
    pause(0.8);
end;

%imwrite(mov, map, 'isosplit_demo.gif', 'DelayTime',0.8, 'LoopCount',0);

end

function isosplit_test
	
close all;

seed0=randi(10000);
seed0=6239;
rng(seed0);
centers={[0,0],[3,3],[-3,5],[9,-6],[0,-6]};
pops={2800,2500,400,90,300};
shapes={[1,1,0],[1,1,0],[1,1,0],[1.3,1.7,0.3],[1.5,1,0]};
opts.split_threshold=0.3;
opts.K=15;

fprintf('seed = %d\n',seed0);

samples=zeros(2,0);
true_labels=zeros(1,0);

for j=1:length(centers)
	xx=randn(1,pops{j});
	yy=randn(1,pops{j});
	shape=shapes{j};
	xx2=xx*shape(1)+yy*shape(3);
	yy2=yy*shape(2)-xx*shape(3);
	center=centers{j};
	xx2=xx2+center(1);
	yy2=yy2+center(2);
	tmp=zeros(2,pops{j});
	tmp(1,:)=xx2; tmp(2,:)=yy2;
	samples=[samples,tmp];
	true_labels=[true_labels,ones(1,pops{j})*j];
end;

colors='rgbkymc';
figure;
for j=1:max(true_labels)
	xx=samples(1,true_labels==j);
	yy=samples(2,true_labels==j);
	col=colors(mod(j-1,length(colors))+1);
	plot(xx,yy,['.',col]);
	if (j==1) hold on; end;
end;
title('Ground truth');
set(gcf,'position',[50,450,400,300]);

% try
% labels=dbscan_daszykowski(samples',4,[])';
% figure;
% for j=1:max(labels)
% 	xx=samples(1,labels==j);
% 	yy=samples(2,labels==j);
% 	col=colors(mod(j-1,length(colors))+1);
% 	plot(xx,yy,['.',col]);
% 	if (j==1) hold on; end;
% end;
% title('dbscan'); set(gcf,'position',[50,0,400,300]);
% catch
% 	warning('Unable to run dbscan_daszykowski');
% end;

labels=local_kmeans_sorber(samples,length(centers));
figure;
for j=1:max(labels)
	xx=samples(1,labels==j);
	yy=samples(2,labels==j);
	col=colors(mod(j-1,length(colors))+1);
	plot(xx,yy,['.',col]);
	if (j==1) hold on; end;
end;
title('kmeans'); set(gcf,'position',[500,0,400,300]);

[labels,info]=isosplit(samples,opts);
figure;
for j=1:max(labels)
	xx=samples(1,labels==j);
	yy=samples(2,labels==j);
	col=colors(mod(j-1,length(colors))+1);
	plot(xx,yy,['.',col]);
	if (j==1) hold on; end;
end;
title('isosplit clustering'); set(gcf,'position',[950,0,400,300]);

fprintf('num clusters = %d\n',max(labels));

fprintf('seed = %d\n',seed0);
	
end

function W=find_svm_discriminant_direction(X,Y)

AA=cat(2,X,Y);
[~,SS]=evalc('svmtrain(transpose(cat(2,ones(1,size(X,2)),ones(1,size(Y,2))*2)),transpose(AA))'); %need to use evalc to suppress verbosity
W=AA(:,SS.sv_indices)*SS.sv_coef;
W=W/sqrt(sum(W.^2));

end