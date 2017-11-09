function [ SugiCorr , SugiR , LM , SugiY , SugiX , origY , origX ] = convergent_cross_map( X , Y , tau , E , LMN, ts_length )

% Calculating Sugihara's CMM, L and M causality. 
%
% References:
% - for Sugihara's CCM method > Sugihara, George, et al., Detecting Causality in Complex Ecosystems, Science 26 October 2012, Vol. 338, no. 6106, pp. 496-500.   
% - for L method > Chicharro, Daniel, and Ralph G. Andrzejak, Reliable detection of directional couplings using rank statistics, Physical Review E 80.2 (2009): 026217.
% - for M method > Andrzejak, Ralph G., et al., Bivariate surrogate techniques: necessity, strengths, and caveats, Physical review E 68.6 (2003): 066202.	
% 
% Inputs:
% X,Y - time series with the same length
% tau - time step for the reconstruction
% E   - dimension of the reconstruction  
% LMN - number of neighborhoods for L and M methods 
%       the number of neighborhoods for Sugihara's CCM method is E+1  
%
% Outputs:
% SugiCorr - correlation between the CCM estimation of original data and original data
% SugiR    - sqrt((sum((originaldata-CCMestimaleddata).^2)/numel(originaldata)))/std(origaldata)
% LM - results for L and M methods
% SugiY, SugiX - the CCM estimate of original data 
% origY, origX - original data

switch nargin
    case 5
        ts_length = 21;
    case 4
        ts_length = 21;
        LMN = E+1;
    otherwise
        error('Bad input')
end
Xm = [];
Ym = [];
count = double(idivide(length(X), int32(ts_length)));
L=ts_length;
T=1+(E-1)*tau;
SugiN=E+1;
N = (L-T+1)*count;
SugiX=zeros(N,1);
SugiY=zeros(N,1);
for ts_idx = 1:count
    X_i = X((ts_idx-1)*ts_length+1:ts_idx*ts_length);
    Y_i = Y((ts_idx-1)*ts_length+1:ts_idx*ts_length);
    Xm_i=zeros((L-T+1),E);
    Ym_i=zeros((L-T+1),E);

    %% RECONTRUCTIONS OF ORIGINAL SYSTEMS

    for t=1:(L-T+1)
        Xm_i(t,:)=X_i((T+t-1):-tau:(T+t-1-(E-1)*tau));
        Ym_i(t,:)=Y_i((T+t-1):-tau:(T+t-1-(E-1)*tau));
    end
    Xm = cat(1,Xm,Xm_i);
    Ym = cat(1,Ym,Ym_i);
end
%%
%{
LMj= zeros(2,2,N);


parfor j=1:N
%% neighborhood search 

[n1,d1]=knnsearch(Xm,Xm(j,:),'k',N);
[n2,d2]=knnsearch(Ym,Ym(j,:),'k',N);

%% LM

LMn1=n1(n1~=j);
LMn2=n2(n2~=j);
LMd1=d1(n1~=j);
LMd2=d2(n2~=j);

susXY=arrayfun(@(x) find(LMn1(:) == x,1,'first'), LMn2(1:LMN) );
susYX=arrayfun(@(x) find(LMn2(:) == x,1,'first'), LMn1(1:LMN) );

sum1=sum(LMd1(:))/(N-1);
sum2=sum(LMd2(:))/(N-1);

LMj(:,:,j) = [(N/2-sum(susXY)/LMN)/(N/2-(LMN+1)/2) , (sum1-sum(LMd1(susXY))/LMN)/(sum1-sum(LMd1(1:LMN))/LMN) ; (N/2-sum(susYX)/LMN)/(N/2-(LMN+1)/2) , (sum2-sum(LMd2(susYX))/LMN)/(sum2-sum(LMd2(1:LMN))/LMN)]; 

   % (GN-G(Y|X))/(GN-Gk) => L(Y|X)
   % (GN-G(X|Y))/(GN-Gk) => L(X|Y)

   %(RNY-Rcond(Y|X))/(RNY-RkY) => M(Y|X)
   %(RNX-Rcond(X|Y))/(RNX-RkX) => M(X|Y)
end
%}
    %% CMM

dat=floor((L-T+1)*count/2);

parfor ii=(dat+1):((L-T+1)*count)
    [n1s,d1s]=knnsearch(Xm((ii-dat):(ii-1),:),Xm(ii,:),'k',SugiN);
    [n2s,d2s]=knnsearch(Ym((ii-dat):(ii-1),:),Ym(ii,:),'k',SugiN);
    u1s=exp(-d1s/d1s(1));
    w1s=u1s/sum(u1s);
    SugiY(ii)= w1s*Y(n1s+double(idivide(int32(n1s-1),L-T+1))*(T-1)+ii-(dat+1));

    u2s=exp(-d2s/d2s(1));
    w2s=u2s/sum(u2s);
    SugiX(ii)= w2s*X(n2s+T-1+ii-(dat+1));
end

origY=Y(T:end);
origY=origY((dat+1):(L-T+1)*count);
SugiY=SugiY((dat+1):(L-T+1)*count);
origX=X(T:end);
origX=origX((dat+1):(L-T+1)*count);
SugiX=SugiX((dat+1):(L-T+1)*count);

SugiCorr1=corrcoef(origY,SugiY);
SugiCorr(2,1)=SugiCorr1(1,2);

SugiCorr2=corrcoef(origX,SugiX);
SugiCorr(1,1)=SugiCorr2(1,2);

SugiR(2,1)=sqrt((sum((origY-SugiY).^2)/numel(origY)))/std(origY);
SugiR(1,1)=sqrt((sum((origX-SugiX).^2)/numel(origX)))/std(origX);

%LM = squeeze(mean(LMj,3));
LM = 0;

end