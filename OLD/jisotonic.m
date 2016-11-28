function [B,MSEs]=jisotonic(A,direction,weights)
% jisotonic - isotonic regression (jfm, may 2105)
%
% [B,MSEs] = jisotonic(A,direction,weights)
%   A is the input a vector
%   direction = 'increading', 'decreasing', 'updown', or 'downup'
%   weights is the optional input weight vector (same size as A)
%   B is the output vector (same size as A)
%   MSEs is used internally for 'updown' and 'downup' directions
%
% Magland 5/19/2015

if (nargin<1)
	jisotonic_test; %run the test code
	return;
end;
if (nargin<2)
	direction='increasing';
end;
if (~isrow(A)) A=A'; end;
if (nargin<3)
	weights=ones(size(A));
end;

if (strcmp(direction,'decreasing'))
	[B,MSEs]=jisotonic(-A,'increasing',weights); B=-B;
	return;
elseif (strcmp(direction,'updown'))
	[B1,MSE1]=jisotonic(A,'increasing',weights);
	[B2,MSE2]=jisotonic(A(end:-1:1),'increasing',weights(end:-1:1));
	B2=B2(end:-1:1);
	MSE2=MSE2(end:-1:1);
	MSE0=MSE1+MSE2;
	
	[~,best_ind]=min(MSE0);
	C1=jisotonic(A(1:best_ind),'increasing',weights(1:best_ind));
	C2=jisotonic(A(best_ind:end),'decreasing',weights(best_ind:end));
	B=[C1(1:best_ind),C2(2:end)];
	if (isnan(B(1)))
		warning('jisotonic: NaN');
	end;
	return;
elseif (strcmp(direction,'downup'))
	B=-jisotonic(-A,'updown',weights);
	return;
else
	if (~strcmp(direction,'increasing'))
		error(['invalid direction in jisotonic: ',direction]);
		return;
	end;
end;

try
[B,MSEs]=jisotonic_mex(A,weights);
catch
	error('Unable to run mex file -- You must compile using: mex jisotonic_mex.cpp');
	%mex(sprintf('%s/jisotonic_mex.cpp',fileparts(mfilename('fullpath'))));
	%[B,MSEs]=jisotonic_mex(A,weights);
end;

%assume increasing
% X=cell(1,0);
% N=length(A);
% 
% tmp.unweightedcount=1;
% tmp.count=weights(1);
% tmp.sum=A(1)*weights(1);
% tmp.sumsqr=A(1)^2*weights(1);
% lastind=1;
% X{lastind}=tmp;
% 
% MSEs=zeros(1,N);
% MSEs(1)=0;
% 
% for j=2:N
% 	tmp.unweightedcount=1;
% 	tmp.count=weights(j);
% 	tmp.sum=A(j)*weights(j);
% 	tmp.sumsqr=A(j)^2*weights(j);
% 	X{lastind+1}=tmp;
% 	MSEs(j)=MSEs(j-1);
% 	lastind=lastind+1;
% 	
% 	while true
% 		if (lastind<=1) break; end;
% 		tmp1=X{lastind-1};
% 		tmp2=X{lastind};
% 		prevMSE=tmp1.sumsqr-tmp1.sum^2/tmp1.count + tmp2.sumsqr-tmp2.sum^2/tmp2.count;
% 		if (tmp1.sum/tmp1.count<tmp2.sum/tmp2.count)
% 			break;
% 		else
% 			tmp.unweightedcount=tmp1.unweightedcount+tmp2.unweightedcount;
% 			tmp.count=tmp1.count+tmp2.count;
% 			tmp.sum=tmp1.sum+tmp2.sum;
% 			tmp.sumsqr=tmp1.sumsqr+tmp2.sumsqr;
% 			X{lastind-1}=tmp;
% 			lastind=lastind-1;
% 			newMSE=tmp.sumsqr-tmp.sum^2/tmp.count;
% 			MSEs(j)=MSEs(j)+newMSE-prevMSE;
% 		end;
% 	end;
% end;
% 
% B=zeros(1,N);
% ii=1;
% for k=1:lastind
% 	tmp0=X{k};
% 	for cc=1:tmp0.unweightedcount;
% 		B(ii+cc-1)=tmp0.sum/tmp0.count;
% 	end;
% 	ii=ii+tmp0.unweightedcount;
% end;

end

function jisotonic_test

A=[1.1,3.2,5.3,2.4,4.5,6.6,3.7,5.8,7.9,10.0,9.2,11.2,7.3,3.4,5.5,2.6,3.7,1.8];
[B,BMSE]=jisotonic(A);
[C,CMSE]=jisotonic(A,'decreasing');
D=jisotonic(A,'updown');
figure; plot(1:length(A),A,'b',1:length(B),B,'r',1:length(C),C,'g',1:length(D),D,'k');
figure; plot(1:length(BMSE),BMSE,'k');

end
