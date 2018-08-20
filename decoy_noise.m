clear
addpath('function');
fspace     = -90:90;
X          = [15,10];
Y          = [10,15];
sig        = 0;
Xsmat      = 10;
trials     = 100;
iterations = 100;
ind        = 3;
mincp      = 0;
maxcp      = 1;
decoytype  = {'attraction - fig.4a','compromise - fig. 4b'};
%%
for c = 1:length(decoytype)
    switch c
        case 1
            Z = [13,5];
            maxtw = 20;
            scaly = 400;
            noisemat = fliplr(linspace(3,7,4));
            
        case 2
            Z      = [20,5];
            maxtw  = 40;
            scaly  = 140;
            noisemat = fliplr(linspace(0.08,0.3,4));
           
    end
    
    clear allcp
    sigmat1 = maxtw-(normpdf(fspace,Z(1),sig+Xsmat)*scaly);
    sigmat1(sigmat1<0.1) =0.1;
    PX1 = normpdf(X(1),fspace,sigmat1);
    PY1 = normpdf(Y(1),fspace,sigmat1);
    PZ1 = normpdf(Z(1),fspace,sigmat1);
    X1  = PX1*fspace';
    Y1  = PY1*fspace';
    Z1  = PZ1*fspace';
    
    allop = [X(2), Y(2) Z(2)];
    sigmat2 = maxtw-(normpdf(fspace,Z(2),sig+Xsmat)*scaly);
    sigmat(sigmat2<0.1) =0.1;
    PX2 = normpdf(X(2),fspace,sigmat2);
    PY2 = normpdf(Y(2),fspace,sigmat2);
    PZ2 = normpdf(Z(2),fspace,sigmat2);
    X2  = PX2*fspace';
    Y2  = PY2*fspace';
    Z2  = PZ2*fspace';
    
    allX = X1+X2;
    allY = Y1+Y2;
    allZ = Z1+Z2;
    
    for iter = 1:iterations
        for n = 1:length(noisemat)
            nX  = (repmat(allX,[trials,1])+randn([trials,1]).*noisemat(n));
            nY  = (repmat(allY,[trials,1])+randn([trials,1]).*noisemat(n));
            nZ  = (repmat(allZ,[trials,1])+randn([trials,1]).*noisemat(n));
            mchX = (nX>nY & nY>nZ);
            mchY = (nY>nZ & nY>nX);
            mchZ = (nZ>nY & nZ>nX);
            cp_x = mean(mchX);
            cp_y = mean(mchY);
            cp_z = mean(mchZ);
            allcp(iter,:,n) = [cp_x,cp_y,cp_z];
        end
    end
    
    figure;
    barsem(allcp(:,1:ind,:));
    ylim([mincp maxcp])
    ylabel('choice probabilitys')
    xlabel('Noise level (high to low)')
    title([decoytype{c}])
end
