
% 2a
load frey.mat
ff = ff(:,1:100);
x = im2double(ff);
[W H] = naiveNMF(x);
for i = 1:10
  subplot(2,5,i); cla;
  imagesc(reshape(W(:,i),20,28)'); grid off; 
  set(gca,'XTick',[]); set(gca,'YTick',[]); title(strcat('Basis ', int2str(i)));
end;

% 2b
load swimmer.mat
X=reshape(Y,1024,256);
[W H] = naiveNMF(X);
for i = 1:10
  subplot(2,5,i); cla;
  imagesc(reshape(W(:,i),32,32)); grid off; 
  set(gca,'XTick',[]); set(gca,'YTick',[]); title(strcat('Basis ', int2str(i)));
end;

% 2c
load frey.mat
x = im2double(ff);
[p q] = r1d_sparseA(x,[],[],20,0.75,0.75);
for i = 1:10
  subplot(2,5,i); cla;
  imagesc(reshape(p(:,i),20,28)'); grid off; 
  set(gca,'XTick',[]); set(gca,'YTick',[]); title(strcat('Basis ', int2str(i)));
end;

for i = 11:20
  subplot(2,5,i-10); cla;
  imagesc(reshape(p(:,i),20,28)'); grid off; 
  set(gca,'XTick',[]); set(gca,'YTick',[]); title(strcat('Basis ', int2str(i)));
end;

% 2d
load swimmer.mat
X=reshape(Y,1024,256);
[p q] = r1d_sparseA(X,[],[],20,0.75,0.75);
for i = 1:10
  subplot(2,5,i); cla;
  imagesc(reshape(p(:,i),32,32)'); grid off; 
  set(gca,'XTick',[]); set(gca,'YTick',[]); title(strcat('Basis ', int2str(i)));
end;

for i = 11:20
  subplot(2,5,i-10); cla;
  imagesc(reshape(p(:,i),32,32)'); grid off; 
  set(gca,'XTick',[]); set(gca,'YTick',[]); title(strcat('Basis ', int2str(i)));
end;