function [avg,stdev,scores]=read_isosplit_calibration(N)
global isosplit_data;

if (~isfield(isosplit_data,'calibration_Ns'))
	dirname=fileparts(get_calibration_fname(1));
	list=dir(sprintf('%s/isosplit_calibration_*.txt',dirname));
	isosplit_data.calibration_Ns=zeros(1,length(list));
	for ii=1:length(list)
		str=list(ii).name;
		N0=sscanf(str(22:end-4),'%d');
		isosplit_data.calibration_Ns(ii)=N0;
	end;
end;

N1=0;
N2=0;
found=false;
N1=max(isosplit_data.calibration_Ns(find(isosplit_data.calibration_Ns<=N)));
N2=max(isosplit_data.calibration_Ns(find(isosplit_data.calibration_Ns>=N)));
if (length(N1)==0)||(length(N2)==0)
	avg=[]; stdev=[]; scores=[];
	return;
end;
fname1=get_calibration_fname(N1);
fname2=get_calibration_fname(N2);

[avg1,stdev1,scores1]=read_isosplit_calibration2(fname1);
[avg2,stdev2,scores2]=read_isosplit_calibration2(fname2);
if (N1==N2)
	avg=avg1; stdev=stdev1; scores=scores1;
else
	pct=(N-N1)/(N2-N1);
	avg=avg1+(avg2-avg1)*pct;
	stdev=stdev1+(stdev2-stdev1)*pct;
	scores=scores1+(scores2-scores1)*pct;
end;

end

function [avg,stdev,scores]=read_isosplit_calibration2(fname)

global isosplit_stored_calibrations;
code=sprintf('C%.12d',string2hash(fname));
if (isfield(isosplit_stored_calibrations,code))
	CC=getfield(isosplit_stored_calibrations,code);
	avg=CC.avg;
	stdev=CC.stdev;
	scores=CC.scores;
	return;
end;

lines=strsplit(fileread(fname),'\n');
ii=1;
curve_len=sscanf(lines{ii},'%d'); ii=ii+1;
num_trials=sscanf(lines{ii},'%d'); ii=ii+1;
avg=zeros(curve_len,1);
stdev=zeros(curve_len,1);
scores=zeros(num_trials,1);
for j=1:curve_len
	tmp=sscanf(lines{ii},'%f,%f'); ii=ii+1;
	avg(j)=tmp(1);
	stdev(j)=tmp(2);
end;
for j=1:num_trials
	scores(j)=sscanf(lines{ii},'%f'); ii=ii+1;
end;

CC.avg=avg;
CC.stdev=stdev;
CC.scores=scores;
isosplit_stored_calibrations.(code)=CC;

end

function fname=get_calibration_fname(N)

dirname=fileparts(mfilename('fullpath'));
fname=sprintf('%s/isosplit_calibration/isosplit_calibration_%d.txt',dirname,N);

end

function hash=string2hash(str,type)
% This function generates a hash value from a text string
%
% hash=string2hash(str,type);
%
% inputs,
%   str : The text string, or array with text strings.
% outputs,
%   hash : The hash value, integer value between 0 and 2^32-1
%   type : Type of has 'djb2' (default) or 'sdbm'
%
% From c-code on : http://www.cse.yorku.ca/~oz/hash.html 
%
% djb2
%  this algorithm was first reported by dan bernstein many years ago 
%  in comp.lang.c
%
% sdbm
%  this algorithm was created for sdbm (a public-domain reimplementation of
%  ndbm) database library. it was found to do well in scrambling bits, 
%  causing better distribution of the keys and fewer splits. it also happens
%  to be a good general hashing function with good distribution.
%
% example,
%
%  hash=string2hash('hello world');
%  disp(hash);
%
% Function is written by D.Kroon University of Twente (June 2010)


% From string to double array
str=double(str);
if(nargin<2), type='djb2'; end
switch(type)
    case 'djb2'
        hash = 5381*ones(size(str,1),1); 
        for i=1:size(str,2), 
            hash = mod(hash * 33 + str(:,i), 2^32-1); 
        end
    case 'sdbm'
        hash = zeros(size(str,1),1);
        for i=1:size(str,2), 
            hash = mod(hash * 65599 + str(:,i), 2^32-1);
        end
    otherwise
        error('string_hash:inputs','unknown type');
end

end
