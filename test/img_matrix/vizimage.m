%
% CS 472/672 Assignment #3 - vizimage.m
%
% It reads the output from main.c and displays the denoised image.
%

load -ascii img.m

N=length(img);
n=sqrt(N);
A=zeros(n,n);
X=zeros(n,n);
for i=1:n,
  X(i,:)=img((i-1)*n+1:i*n)';
end;

figure(1);
colormap(gray);
imagesc(X);
axis('square');

