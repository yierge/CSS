function [P_chern,P_com,C_sig,phi,f,var_c,T_max]=CSS_2D(...
    voxel_p,image_data,mask,rho,p_val)
%INPUTS
%voxel_p: voxel-wise p-values in a 1 by V vector
%image_data: input images, V(num of voxels) by N(subjects) matrix
%mask: mask used for the input image, aim for obtaining the dim info
%rho: smoothness level, put 4, 6, or 8
%p_val: primary threshold



%OUTPUTS
%P_chern: CSS for all clusters
%P_com: original combined p-value for all clusters
%C_sig: all clusters
%phi: -sum(log(p)) for each cluster
%f: approximated degree of freedom
%var_c: variance of a*chi_f^2
%cluster_size: size of all clusters


data=image_data;

P1d=-log(voxel_p);
P1d(~isfinite(P1d))=36.0437;

%thresholding
sig1d1=voxel_p<p_val;

datatest=zeros(size(mask));
datatest(sig1d1)=1;
CC=bwconncomp(datatest);
[a,~]=cellfun(@size,CC.PixelIdxList(:));
dd=CC.PixelIdxList(a>1);

C_sig=dd;


%explaination of the variables
%sum_rho: the summation of correlations on off diag elements of the
%correlation matrix. Refer to corr
%var_c: variance of the new distribution phi
%exp_c: expectation of the new distribution phi



P_chern=zeros(length(C_sig),1);
P_com=zeros(length(C_sig),1);
var_c = zeros(length(C_sig),1);
exp_c=zeros(length(C_sig),1);
f=zeros(length(C_sig),1);
c=zeros(length(C_sig),1);
phi=zeros(length(C_sig),1);
ratio=zeros(length(C_sig),1);
test_rate=zeros(length(C_sig),1);
cluster_size=zeros(length(dd),1);

for j=1:length(C_sig)
    
    off_corr_m= corrcoef(data(C_sig{j},:)');
    th_corr_m=variogram_corr_real(C_sig{j},off_corr_m,rho);
    sum_rho1=sum(sum(th_corr_m))-length(C_sig{j});
    th_corr_m2=th_corr_m.*th_corr_m;
    sum_rho2=sum(sum(th_corr_m2))-length(C_sig{j});
    th_corr_m3=th_corr_m2.*th_corr_m;
    sum_rho3=sum(sum(th_corr_m3))-length(C_sig{j});
    
    var_c(j)=4*length(C_sig{j})+3.263*sum_rho1+0.71*sum_rho2+0.027*sum_rho3;
    exp_c(j)=2*length(C_sig{j});
    
    f(j)=(2*exp_c(j)^2)/var_c(j);
    c(j)=var_c(j)/(2*exp_c(j));
    phi(j)=2*sum(P1d(C_sig{j}));
    ratio(j)=phi(j)/c(j);
    
    test_rate(j)=phi(j)/length(C_sig{j});
    
    P_com(j)=1-chi2cdf(ratio(j),f(j));
    P_chern(j)=(ratio(j)/(f(j)))^(f(j)/2)*exp(f(j)/2-ratio(j)/2);
    cluster_size(j)=length(C_sig{j});
    
    
end


if isempty(dd)
    P_chern=1;
    P_com=1;
    phi=1;
    f=1;
    var_c=1;
end


end