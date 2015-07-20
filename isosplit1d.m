function [labels,pp]=isosplit1d(X,opts)

if (nargin<1) test_isosplit1d; return; end;
if (isstr(X))
	if (strcmp(X,'calibrate'))
		do_calibration(opts);
	else
		error('Invalid string command.');
	end;
	return;
end;

if (nargin<2) opts=struct; end;
opts.minsize=4; %hard coded to be consistent with calibration.

N=length(X);

%read the calibration (speed this up?)
tic1=tic;
[avg,stdev,scores]=read_isosplit_calibration(N);
if (isfield(opts,'m_max'))
	avg=avg(1:opts.m_max);
	stdev=stdev(1:opts.m_max);
end;
curve_len=length(avg);
num_trials=length(scores);
if (length(avg)==0) 
	error(sprintf('isosplit1d: calibration not performed for N=%d\n',N));
end;
%fprintf('tic1: %f\n',toc(tic1));

%compute the curve
tic2=tic;
[curve,cutpoint]=get_isosplit_curve(X,curve_len,opts);
%fprintf('tic2: %f\n',toc(tic2));

tic3=tic;
%compute the score by comparing to calibration avg and stdev
% and then compute the pp
inds=find(stdev~=0);
if (length(inds)>0)
	diff0=(avg(inds)-curve(inds)');
	%we need a change of at least 0.1 to count it
	%this is important in the case where stdev is extremely close to 0
	diff0(abs(diff0)<0.1)=0;
	score0=max(diff0./stdev(inds));
	pp=length(find(scores<score0))/num_trials;
else
	score0=0;
	pp=0;
end;
%fprintf('tic3: %f\n',toc(tic3));

tic4=tic;
labels=zeros(1,length(X));
labels(find(X<cutpoint))=1;
labels(find(X>=cutpoint))=2;
%fprintf('tic4: %f\n',toc(tic4));

end

function fname=get_calibration_fname(N)

dirname=fileparts(mfilename('fullpath'));
fname=sprintf('%s/isosplit_calibration/isosplit_calibration_%d.txt',dirname,N);

end


function [curve,cutpoint]=get_isosplit_curve(X,curve_len,opts)
	N=length(X);
	if (N<opts.minsize*2)
		curve=zeros(1,curve_len);
		cutpoint=0;
		return;
	end;
	X=sort(X);
	spacings=[X(2)-X(1),(X(3:end)-X(1:end-2))/2,X(end)-X(end-1)];
	density=1./spacings;
	logdensity=log(density);
	fit1=jisotonic(logdensity,'updown');
	resid1=logdensity-fit1;
	%fit2=jisotonic(resid1,'downup');
	%fit2=jisotonic(resid1,'downup',spacings);
	fit2=jisotonic(resid1,'downup',sqrt(spacings)); %this weighting seems to work best
	fit2=fit2-mean(fit2); %added 6/15/15, not sure if necessary
	%fit2=jisotonic(resid1,'downup');
	fit2_sorted=sort(fit2(opts.minsize:end-opts.minsize));
	curve=cumsum(fit2_sorted)./(1:length(fit2_sorted));
	if (length(curve)<curve_len)
		curve=[curve,zeros(1,curve_len-length(curve))];
	end;
	curve=curve(1:curve_len);
	[~,ind]=min(fit2(opts.minsize:end-opts.minsize));
	ind=ind+opts.minsize-1;
	if (ind<=1) cutpoint=X(1); return; end;
	if (ind>=length(X)) cutpoint=X(end); return; end;
	if (spacings(ind-1)>spacings(ind+1))
		cutpoint=(X(ind-1)+X(ind))/2;
	else
		cutpoint=(X(ind)+X(ind+1))/2;
	end;
end


function do_calibration(start_kk)

dirname=fileparts(get_calibration_fname(1));
if (~exist(dirname,'dir'))
	mkdir(dirname);
end;

curve_len=1000;
num_trials=1000;
incr=1.02;

kk=start_kk;
while true
	N=floor(incr^kk);
	fname=get_calibration_fname(N);
	if (~(exist(fname,'file')==2))
		fprintf('N=%d, kk=%d\n',N,kk);
		tA=tic;
		[avg,stdev,scores]=do_calibration2(curve_len,num_trials,N);
		F=fopen(fname,'w');
		fprintf(F,'%d\n',curve_len);
		fprintf(F,'%d\n',num_trials);
		for j=1:curve_len
			fprintf(F,'%.4f,%.4f\n',avg(j),stdev(j));
		end;
		for j=1:num_trials
			fprintf(F,'%.4f\n',scores(j));
		end;
		fclose(F);
		elapsed=toc(tA);
		%fprintf('Elapsed: %f\n',elapsed);
		est_time_factor=elapsed/N;
		est_time=0;
		for kk2=kk:700
			est_time=est_time+incr^kk2*est_time_factor;
		end;
		fprintf('Estimated #hours remaining: %f\n',est_time/60/60);
	end;
	kk=kk+1;
end;

end

function [avg,stdev,scores]=do_calibration2(curve_len,num_trials,N)

A=zeros(curve_len,num_trials);
for trial=1:num_trials
	X=rand(1,N);
	opts.minsize=4;
	curve=get_isosplit_curve(X,curve_len,opts);
	A(:,trial)=curve;
end;
Asorted=sort(A,2);
avg=squeeze(mean(A,2)); %curve_len x 1
stdev=squeeze(sqrt(var(A,[],2))); %curve_len x 1
scores=zeros(num_trials,1); %num_trials x 1
inds=find(stdev~=0);
if (length(inds)>0)
	scores(:)=max((repmat(avg(inds,1),1,num_trials)-A(inds,:))./repmat(stdev(inds,1),1,num_trials),[],1);
end;
scores=sort(scores,1);

end

function test_isosplit1d 

close all;
N=100000;
for ct=1:15
	X=[randn(1,N/2),randn(1,N/2)+ct/3];
	opts.verbose=0;
	tA=tic;
	[labels,pp]=isosplit1d(X,opts);
	toc(tA)

	if (pp<0.9) labels=ones(size(labels)); end;
	colors='rgbkymc';
	figure;
	[~,bins0]=hist(X,length(X)/3);
	num_labels=max(labels);
	if (pp<0.9) labels=ones(size(labels)); end;
	for j=1:num_labels
		[v0,b0]=hist(X(labels==j),bins0);
		h=bar(b0,v0);
		col=colors(mod(j-1,length(colors))+1);
		set(h,'FaceColor',col,'EdgeColor',col);
		if (j==1) hold on; end;
	end;
	title(sprintf('pp = %.2f\n',pp));
	pause(0.1);
end;

end
