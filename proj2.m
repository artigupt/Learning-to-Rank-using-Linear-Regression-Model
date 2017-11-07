UBitName = 'artigupt';

personNumber = '50170010';

[num,txt,raw]=xlsread('dataset.xlsx');
target=xlsread('dataset2.xlsx');

train=num(1:floor(0.8*size(num,1))+1,:);
validation=num(ceil(0.8*size(num,1))+1:ceil(0.8*size(num,1))+floor(.1*size(num,1)),:);
test=num(ceil(0.8*size(num,1))+ceil(.1*size(num,1)):size(num,1),:);
trainInd1 = (1:floor(0.8*size(num,1))+1)';
validInd1 = (ceil(0.8*size(num,1))+1:ceil(0.8*size(num,1))+floor(.1*size(num,1)))';
target1=target(1:floor(0.8*size(target,1))+1,:);
target2=target(ceil(0.8*size(target,1))+1:ceil(0.8*size(target,1))+floor(.1*size(target,1)),:);
target3=target(ceil(0.8*size(target,1))+ceil(.1*size(target,1)):size(target,1),:);

M1=5;
lambda1 = .1;

Sigma1=diag(var(train)+0.1);
for j = 1:M1,
    Sigma1(:,:,j) = diag(var(train)+0.1);
end

mu1=zeros(size(train,2),M1);
for j=(1:M1)
    mu1(:,j)=train(randi([1,size(train,1)],1),:)';
end;

phi=zeros(size(train,1),M1);
for i=1:size(train,1),
    for j=1:M1,
        if j==1,
            phi(i,j)=1;
        else phi(i,j)= exp((-1/2)*(train(i,:)'-mu1(:,j))'*inv(Sigma1(:,:,j))*(train(i,:)'-mu1(:,j)));
        end
    end
end;

w1=inv(lambda1*eye(size(phi,2)) + phi'*phi)*phi'*target1;

error=((1/2)*sum((target1'-(w1'*phi')).^2)) + lambda1*w1'*w1;

trainPer1=sqrt((2*error)/size(train,1));

%-----------------------------validation Set---------------------%

phi2=zeros(size(validation,1),M1);
for i=1:size(validation,1),
    for j=1:M1,
        if j==1,
            phi2(i,j)=1;
        else phi2(i,j)= exp((-1/2)*(validation(i,:)'-mu1(:,j))'*inv(Sigma1(:,:,j))*(validation(i,:)'-mu1(:,j)));
        end
    end
end


error2=(1/2)*((target2-phi2*w1)'*(target2-phi2*w1));
validPer1=sqrt((2*error2)/size(validation,1));


%------------------------- Synthetic Training Data------------------------%

X=x';

train2=X(1:floor(0.9*size(X,1)),:);
validation2=X(floor(0.9*size(X,1))+1:end,:);
trainInd2 = (1:floor(0.9*size(X,1)))';
validInd2 = (floor(0.9*size(X,1))+1:size(X,1))';
Starget1=t(1:floor(0.9*size(t,1)),:);
Starget2=t(floor(0.9*size(t,1))+1:end,:);

M2=6;
lambda2 = 1;

Sigma2=diag(var(train2)+0.1);
for j = 1:M2,
    Sigma2(:,:,j) = diag(var(train2)+0.1);
end

mu2=zeros(size(train2,2),M2);
for j=(1:M2)
    mu2(:,j)=train2(randi([1,size(train2,1)],1),:)';
end;

syntheticphi=zeros(size(train2,1),M2);
for i=1:size(train2,1),
    for j=1:M2,
        if j==1,
            syntheticphi(i,j)=1;
        else syntheticphi(i,j)= exp((-1/2)*(train2(i,:)'-mu2(:,j))'*inv(Sigma2(:,:,j))*(train2(i,:)'-mu2(:,j)));
        end
    end
end;

w2=inv(lambda2*eye(size(syntheticphi,2)) + syntheticphi'*syntheticphi)*syntheticphi'*Starget1;

Serror=(1/2)*((Starget1-syntheticphi*w2)'*(Starget1-syntheticphi*w2));
trainPer2=sqrt((2*Serror)/size(train2,1));

%-----------------------------Synthetic validation Set---------------------%

Sphi2=zeros(size(validation2,1),M2);
for i=1:size(validation2,1),
    for j=1:M2,
        if j==1,
            Sphi2(i,j)=1;
        else Sphi2(i,j)= exp((-1/2)*(validation2(i,:)'-mu2(:,j))'*inv(Sigma2(:,:,j))*(validation2(i,:)'-mu2(:,j)));
        end
    end
end;

Serror2=(1/2)*((Starget2-Sphi2*w2)'*(Starget2-Sphi2*w2));
validPer2=sqrt((2*Serror2)/size(validation2,1));

%----------------------- Stochastic Closed Form Dataset--------------------------%


size_phi = size(phi,2);
w01 = zeros(size_phi,1);
dw=[];
 eta1=1;
dw1=[];
for i=1:size(phi,1)
    dw=eta1*((target1(i)- w01'*phi(i,:)')*phi(i,:)'-lambda1*w01);
    
    w01 = w01+dw;
    dw1 = [dw1 dw];
end



%----------------------- Stochastic Synthetic Dataset--------------------------%


size_phi2 = size(syntheticphi,2);
w02 = zeros(size_phi2,1);

dws=[];
eta2=1;
dw2=[];
for i=1:size(syntheticphi,1)
    dws=eta2*((Starget1(i)- w02'*syntheticphi(i,:)')*syntheticphi(i,:)'-lambda2*w02)
  
    w02 = w02+dws;
    dw2 = [dw2 dws];
end

save('proj2.mat', 'UBitName', 'personNumber', 'dw', 'dw1', 'dw2', 'dws', 'error', 'error2', 'eta1', 'eta2', 'i', 'j', 'lambda1', 'lambda2', 'M1', 'M2', 'mu1', 'mu2','num', 'phi', 'phi2', 'raw', 'Serror', 'Serror2', 'Sigma1','Sigma2','size_phi', 'size_phi2', 'Sphi2', 'Starget1', 'Starget2', 'syntheticphi', 't', 'target', 'target1', 'target2', 'target3', 'test', 'train','train2', 'trainInd1', 'trainInd2', 'training_size_col', 'trainPer1', 'trainPer2', 'validation','validation2','validInd1', 'validInd2', 'validPer1', 'validPer2','w01','w02', 'w1', 'w2', 'x','X' );
 










