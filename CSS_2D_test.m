clear

%%%-----making sample image-----%%%
%create a mask for dim detection
mask=zeros(100,100);

%smoothness parameters
newsigma=8/(2*sqrt(8*log(2)));

%effect size (group difference)
es=0.2;

%generate images
data=zeros(100,100);
data(37:46,37:46)=1; 
data(56:60,40:42)=1;

data1d=data(:);

%# subjects per arm
arm=60;

%primary threshold
p_val=0.001;
%smoothness parameter
smt=4;

image_data=randn(10000,arm*2);
image_data(data1d>0,1:arm)=image_data(data1d>0,1:arm)+es;

image_data_blur=zeros(size(image_data));
for k=1:arm*2
    image_data2d_sub=reshape(image_data(:,k),size(data));
    image_data2d_sub_blur = imgaussfilt(image_data2d_sub, newsigma);
    image_data_blur(:,k)= reshape(image_data2d_sub_blur,numel(image_data2d_sub_blur ),1) ;
end

T=zeros(size(data1d));
P=zeros(size(data1d));
for j=1: (size(image_data_blur,1))
    x=image_data_blur(j,1:arm);
    y=image_data_blur(j,(arm+1):arm*2);
    T(j)=(mean(x)-mean(y))/sqrt(var(x)/arm+var(y)/arm);
    P(j)=2*(1-tcdf(abs(T(j)),2*arm-2));
end

figure;imagesc(reshape(T,100,100));colormap jet;colorbar;

[P_chern,P_com,C_sig,phi,f,var_c,cluster_size]=CSS_2D(P,image_data_blur,mask,smt,p_val);


%permutation test
M=100;
P_min=zeros(M,1);

parfor m=1:M
    group_perm1= randsample(arm*2,arm);
    group_perm0=setdiff(1:arm*2,group_perm1);
    order=[group_perm1' group_perm0];
    P1=zeros(size(data1d));
    T1=zeros(size(data1d));
    for j=1: (size(image_data_blur,1))
        x=image_data_blur(j,group_perm1);
        y=image_data_blur(j,group_perm0);
        T1=(mean(x)-mean(y))/sqrt(var(x)/arm+var(y)/arm);
        P1(j)=2*(1-tcdf(abs(T1),2*arm-2));
    end
    
    [P_chern1,~,~,~,~,~,cluster_size1]=CSS_2D(P1,image_data_blur(:,order),mask,smt,p_val);
    
    P_min(m)=min(P_chern1);
    T_max(m)=max(cluster_size1);
    
    m
end

%output CSS threshold
P_min=sort(P_min);
p_chrf=P_min(.05*M+1)
%output cluster-size threshold
T_max=sort(T_max);
perm_size=T_max(.95*M-1)

%output indices
idx=C_sig{1,P_chern<=p_chrf};

close
