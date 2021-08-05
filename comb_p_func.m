
function [P_chern,P_com,C_sig,P,T,image_data_blur,phi,test_rate,f,ratio,var_c,avg_p,avg_dep]=comb_p_func(image_data,truth,sampsize,sigma,es1,es2,p_val)
%generate datasets
%image_data=randn(10000,sampsize*2);
%smooth first and add true signal
data=zeros(100,100);
image_data_blur=[];
for i=1:sampsize*2
 image_data2d=reshape(image_data(:,i),size(data));
 image_data2d_blur = imgaussfilt(image_data2d, sigma);   
 image_data_blur(:,i)= reshape(image_data2d_blur,numel(image_data2d_blur ),1) ; 
end

%figure;imagesc(reshape(image_data_blur(:,1),100,100));colormap jet;colorbar;%type 2: add true signal first and then smooth

image_data_blur(find(truth==1),1:sampsize)=image_data_blur(find(truth==1),1:sampsize)+es1;
image_data_blur(find(truth==2),1:sampsize)=image_data_blur(find(truth==2),1:sampsize)+es2;

%figure;imagesc(reshape(image_data_blur(:,1),100,100));colormap jet;colorbar;%type 2: add true signal first and then smooth

%two-sample ttest

P=[];T=[];
for j=1: (size(image_data_blur,1))
     x=image_data_blur(j,1:sampsize);
    y=image_data_blur(j,(sampsize+1):sampsize*2);
    T(j)=(mean(x)-mean(y))/sqrt(var(x)/sampsize+var(y)/sampsize);
    P(j)=2*(1-tcdf(abs(T(j)),2*sampsize-2));
end
    
%figure;imagesc(reshape(T,100,100));colormap jet;colorbar;

%calculate the correlation
data=[];corr=[];pval=[];P1d=[];
data=image_data_blur;
%corr=corrcoef(data');
pval=P;

%data1d=data(:);
P1d=-log(pval);
P1d(~isfinite(P1d))=36.0437;
sig1d1=[];datatest=[];CC=[];dd=[];C_sig=[];a=[];
    %signficant index
     sig1d1=find(pval<p_val);    
    datatest=zeros(100,100);
    datatest(sig1d1)=1;
    CC=bwconncomp(datatest);
    [a,~]=cellfun(@size,CC.PixelIdxList(:));
    dd=CC.PixelIdxList(a>1);%dd is C_sig
    
    C_sig=dd;

    
    %explaination of the variables
    %sum_rho: the summation of correlations on off diag elements of the
    %correlation matrix. Refer to corr
    %var_c: variance of the new distribution phi
    %exp_c: expectation of the new distribution phi
    %
    
    
    sum_rho2=[];var_c=[];exp_c=[];f=[];c=[];phi=[];ratio=[];P_com=[];test_rate=[];sum_rho1=[];sum_rho3=[];

for j=1:length(C_sig)
    %sum_rho(j)=sum(sum(triu(corr(C_sig{j},C_sig{j}))))-length(C_sig{j});
    off_corr_m=[];off_corr_m2=[];off_corr_m3=[];
      off_corr_m=corrcoef(data(C_sig{j},:)');
      sum_rho1(j)=sum(sum(off_corr_m))-length(C_sig{j});
      off_corr_m2=off_corr_m.*off_corr_m;
      sum_rho2(j)=sum(sum(off_corr_m2))-length(C_sig{j});
      off_corr_m3=off_corr_m2.*off_corr_m;
      sum_rho3(j)=sum(sum(off_corr_m3))-length(C_sig{j});
%    sum_rho(j)=sum(sum(corrcoef(data(C_sig{j},:)')))-length(C_sig{j});
  
%    sum_rho(j)=sum(sum(triu(corrcoef(data(C_sig{j},:)'))))-length(C_sig{j});
    var_c(j)=4*length(C_sig{j})+3.263*sum_rho1(j)+0.71*sum_rho2(j)+0.027*sum_rho3(j);
    exp_c(j)=2*length(C_sig{j});
f(j)=(2*exp_c(j)^2)/var_c(j);
c(j)=var_c(j)/(2*exp_c(j));
phi(j)=2*sum(P1d(C_sig{j}));
ratio(j)=phi(j)/c(j);

avg_p(j)=mean(P(C_sig{j}));
avg_dep(j)=sum_rho1(j)/(length(C_sig{j})*(length(C_sig{j})-1));

%this is test para
test_rate(j)=phi(j)/length(C_sig{j});
%if ratio(j)>1
%   P_com(j)=((ratio(j)/(2*f(j)))^f(j))*exp(f(j)-ratio(j)/2);
%else
P_com(j)=1- chi2cdf(ratio(j),f(j));    
P_chern(j)=(ratio(j)/(f(j)))^(f(j)/2)*exp(f(j)/2-ratio(j)/2);
%end 
   

end

if isempty(dd)==1
    P_chern=1;
    P_com=1;
    sum_rho=1;
    phi=1;
    test_rate=1;
end
    
    
end