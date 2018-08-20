clear
clc
addpath('function')
fspace     = -90:90;
X          = [10,5];
Y          = [5,10];
Z          = [3,12];
Xs         = 35;
trials     = 100;
iterations = 100;
ind        = 3;
mincp      = 0.1;
maxcp      = 0.6;
lown       = 0.02;
highn      = 0.06;
maxtw      = 10;
decoytype  = {'Inference task - fig.S3A','Perceptual task - fig. S3B'};

for c = 1:length(decoytype)
    clear normcp allcp
    switch c
        case 1
            scaly  = 800;
        case 2
            scaly = 600;
    end
    
    noisemat = fliplr(linspace(lown,highn,3));
    allop = [X(1), Y(1) Z(1)];
    sig   = std(allop);
    sigmat1 = maxtw-(normpdf(fspace,Z(1),sig+Xs)*scaly);
    sigmat1(sigmat1<0.1)=0.1;
    PX1 = normpdf(X(1),fspace,sigmat1);
    PY1 = normpdf(Y(1),fspace,sigmat1);
    PZ1 = normpdf(Z(1),fspace,sigmat1);
    X1  = PX1*fspace';
    Y1  = PY1*fspace';
    Z1  = PZ1*fspace';
    
    allop = [X(2), Y(2) Z(2)];
    sig   = std(allop);
    sigmat2 = maxtw-(normpdf(fspace,Z(2),sig+Xs)*scaly);
    sigmat2(sigmat2<0.1)=0.1;
    PX2 = normpdf(X(2),fspace,sigmat2);
    PY2 = normpdf(Y(2),fspace,sigmat2);
    PZ2 = normpdf(Z(2),fspace,sigmat2);
    X2  = PX2*fspace';
    Y2  = PY2*fspace';
    Z2  = PZ2*fspace';
    allX = X1+X2;
    allY = Y1+Y2;
    allZ = Z1+Z2;
    
    vd_XY = allX-allY;
    vd_YX = allY-allX;
    vd_XZ = allX-allZ;
    vd_ZX = allZ-allX;
    vd_YZ = allY-allZ;
    vd_ZY = allZ-allY;
    
    for n = 1:length(noisemat)
        noise = noisemat(n);
        cp_xy = sigmoidv(vd_XY,0,1,0,noise);
        cp_yx = sigmoidv(vd_YX,0,1,0,noise);
        cp_xz = sigmoidv(vd_XZ,0,1,0,noise);
        cp_zx = sigmoidv(vd_ZX,0,1,0,noise);
        cp_yz = sigmoidv(vd_YZ,0,1,0,noise);
        cp_zy = sigmoidv(vd_ZY,0,1,0,noise);
        
        cp_x = mean([cp_xy,cp_xz]);
        cp_y = mean([cp_yx,cp_yz]);
        cp_z = mean([cp_zx,cp_zy]);
        
        allcp(n,:) = [cp_x,cp_y,cp_z];
    end
    normcp = allcp./sum(allcp,2);
    
    f = figure;
    set(f,'Units','inches','position',[0,0,5,5])
    bar(normcp(:,1:2))
    ylim([mincp maxcp])
    ylabel('choice probabilitys')
    xlabel('Noise level (high to low)')
    title([decoytype{c}])
end
