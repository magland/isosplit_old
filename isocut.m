function cutpoint=isocut(samples1d,opts,varargin)

if (nargin<1) isocut('bootstrap'); return; end;
if (nargin<2) opts=struct; end;
if (nargin<3) varargin={}; end;
if (isstr(samples1d)) 
    param1=samples1d; param2=opts;
    fprintf('Running example: %s\n',param1);
    isocut_example(param1,param2,varargin{:}); return;
end;
if (~isfield(opts,'threshold')) opts.threshold=1.4; end;
if (~isfield(opts,'minsize')) opts.minsize=4; end;
if (~isfield(opts,'single_test')) opts.single_test=0; end;
if (~isfield(opts,'show_histogram')) opts.show_histogram=0; end;
if (~isfield(opts,'verbose')) opts.verbose=0; end;

N=length(samples1d);
samples=sort(samples1d);

N0s=[];
for ii=2:floor(log2(N/2))
    N0s(end+1)=2^ii;
    N0s(end+1)=-2^ii;
end;
N0s(end+1)=N;
if (opts.single_test) N0s=N; end;

timerA=tic;
cutpoint=inf;
for j=1:length(N0s)
    if opts.verbose fprintf('.'); end;
    N0=N0s(j);
    if (N0>0)
        samples0=samples(1:N0);
    else
        samples0=samples(end-abs(N0)+1:end);
    end;
    spacings0=diff(samples0);
    spacings0_fit=jisotonic(spacings0,'downup');
    samples0_fit=cumsum([samples0(1),spacings0_fit]);
    ks0=compute_ks(samples0,samples0_fit)*sqrt(length(samples0));
    if (ks0>=opts.threshold)
        spacings1=spacings0./spacings0_fit;
        spacings1_fit=jisotonic(spacings1,'updown');
        if (N0>=opts.minsize*2)
            [~,ind1]=max(spacings1_fit(opts.minsize:end-opts.minsize+1));
            ind1=ind1+opts.minsize-1;
            cutpoint=(samples0(ind1)+samples0(ind1+1))/2;
            break;
        end;
    end;
end;
if (opts.verbose) fprintf('\n'); end;

if (opts.show_histogram)
    show_histogram(samples,cutpoint);
    title(sprintf('ks = %g',ks0));
end;

end

function ks=compute_ks(S1,S2)
N1=length(S1);
N2=length(S2);
dists=zeros(1,N2);
ii=0;
for j=1:N2
    while ((ii+1<=N1)&&(S1(ii+1)<=S2(j))) ii=ii+1; end;
    dists(j)=abs((j/N2)-(ii/N1));
end;
[ks,ind]=max(dists);
end


function isocut_example(example_name,opts)
if (nargin<2) opts=[]; end;

if (strcmp(example_name,'normal'))
    if (~isfield(opts,'N')) opts.N=5000; end;
    samples=randn(1,opts.N);
elseif (strcmp(example_name,'bimodal'))
    if (~isfield(opts,'N1')) opts.N1=5000; end;
    if (~isfield(opts,'N2')) opts.N2=5000; end;
    if (~isfield(opts,'separation')) opts.separation=3.5; end;
    samples=[randn(1,opts.N1),randn(1,opts.N2)+opts.separation];
elseif (strcmp(example_name,'trimodal'))
    if (~isfield(opts,'N1')) opts.N1=1000; end;
    if (~isfield(opts,'N2')) opts.N2=2000; end;
    if (~isfield(opts,'N3')) opts.N3=3000; end;
    if (~isfield(opts,'separation1')) opts.separation1=3.5; end;
    if (~isfield(opts,'separation2')) opts.separation2=3.5; end;
    samples=[randn(1,opts.N1),randn(1,opts.N2)+opts.separation1,randn(1,opts.N3)+opts.separation1+opts.separation2];
elseif (strcmp(example_name,'small_large'))
    opts.N1=100; opts.N2=10000;
    opts.separation=8;
    isocut_example('bimodal',opts);
    opts.single_test=1;
    isocut_example('bimodal',opts);
    return;
elseif (strcmp(example_name,'bootstrap'))
    bootstrap_test(opts);
    return;
end;

opts.verbose=1;
opts.show_histogram=1;
isocut(samples,opts);

end

function bootstrap_test(opts)
if (~isfield(opts,'N')) opts.N=1000; end;
if (~isfield(opts,'separation')) opts.separation=3.5; end;
if (~isfield(opts,'num_trials')) opts.num_trials=100; end;
if (~isfield(opts,'uniform_baseline')) opts.uniform_baseline=0; end;

fields=fieldnames(opts);
for ii=1:numel(fields)
  fprintf('%s = %s\n',fields{ii},to_string(opts.(fields{ii})));
end

num_cut=0;
for j=1:opts.num_trials
    if opts.uniform_baseline
        samples=rand(1,opts.N);
    else
        samples=randn(1,opts.N);
    end;
    cutpoint=isocut(samples,opts);
    if ~isinf(cutpoint)
        num_cut=num_cut+1;
        fprintf('+');
    else
        fprintf('_');
    end;
    if (j==1) show_histogram(samples,cutpoint); end;
end
fprintf('\n');
title(sprintf('%d/%d (%g%%)',num_cut,opts.num_trials,num_cut/opts.num_trials*100));

num_cut=0;
for j=1:opts.num_trials
    samples=[randn(1,opts.N*0.5),randn(1,opts.N*0.5)+opts.separation];
    cutpoint=isocut(samples,opts);
    if ~isinf(cutpoint)
        num_cut=num_cut+1;
        fprintf('+');
    else
        fprintf('_');
    end;
    if (j==1) show_histogram(samples,cutpoint); end;
end
fprintf('\n');
title(sprintf('%d/%d (%g%%)',num_cut,opts.num_trials,num_cut/opts.num_trials*100));

end

function show_histogram(samples,cutpoint)
figure;
[counts,bins]=hist(samples,ceil(length(samples)/6));
bar(bins,counts,'FaceColor',[0.2,0.2,0.2],'EdgeColor',[0.5,0.5,0.5]);
hold on;
plot([cutpoint,cutpoint],ylim,'Color',[0.8,0.3,0.2],'Linewidth',3);
end

function str=to_string(X)
if (isnumeric(X)) str=sprintf('%g',X);
elseif (isstr(X)) str=X
else str='';
end;
end