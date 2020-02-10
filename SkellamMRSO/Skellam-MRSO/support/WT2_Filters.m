function [Fd,Fs] = udwt_filters(Htype,M)

% generate filter coefficients for the wavelet transform

[L0,H0,L1,H1]=wfilters(Htype);
hd = sum(2.^[0:M-1])*(length(L0))+1; % size of the largest filter
Hd = zeros(hd,hd,3*M+1);

% initialize
LL = [L0 0]'*[L0 0];
ll = [L0 0]'*[L0 0];
lh = [L0 0]'*[H0 0];
hl = [H0 0]'*[L0 0];
hh = [H0 0]'*[H0 0];
Hd((hd-size(lh,1))/2+1:end-(hd-size(lh,1))/2,(hd-size(lh,1))/2+1:end-(hd-size(lh,1))/2,2  )=lh;
Hd((hd-size(hl,1))/2+1:end-(hd-size(hl,1))/2,(hd-size(hl,1))/2+1:end-(hd-size(hl,1))/2,3  )=hl;
Hd((hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,(hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,4  )=hh;
k  = 5;
for m=1:M-1
    ll = conv2(up(ll),LL); 
    lh = conv2(up(lh),LL); Hd((hd-size(lh,1))/2+1:end-(hd-size(lh,1))/2,(hd-size(lh,1))/2+1:end-(hd-size(lh,1))/2,k  )=lh;
    hl = conv2(up(hl),LL); Hd((hd-size(hl,1))/2+1:end-(hd-size(hl,1))/2,(hd-size(hl,1))/2+1:end-(hd-size(hl,1))/2,k+1)=hl;
    hh = conv2(up(hh),LL); Hd((hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,(hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,k+2)=hh;
    k = k+3;
end
Hd(:,:,1)=ll;

hd = sum(2.^[0:M-1])*(length(L1))+1; % size of the largest filter
Hs = zeros(hd,hd,3*M+1);
% initialize
LL = [0 L1]'*[0 L1];
ll = [0 L1]'*[0 L1];
lh = [0 L1]'*[0 H1];
hl = [0 H1]'*[0 L1];
hh = [0 H1]'*[0 H1];
Hs((hd-size(lh,1))/2+1:end-(hd-size(lh,1))/2,(hd-size(lh,1))/2+1:end-(hd-size(lh,1))/2,2  )=lh;
Hs((hd-size(hl,1))/2+1:end-(hd-size(hl,1))/2,(hd-size(hl,1))/2+1:end-(hd-size(hl,1))/2,3  )=hl;
Hs((hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,(hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,4  )=hh;
k  = 5;
for m=1:M-1
    ll = conv2(up(ll),LL); 
    lh = conv2(up(lh),LL); Hs((hd-size(lh,1))/2+1:end-(hd-size(lh,1))/2,(hd-size(lh,1))/2+1:end-(hd-size(lh,1))/2,k  )=lh;
    hl = conv2(up(hl),LL); Hs((hd-size(hl,1))/2+1:end-(hd-size(hl,1))/2,(hd-size(hl,1))/2+1:end-(hd-size(hl,1))/2,k+1)=hl;
    hh = conv2(up(hh),LL); Hs((hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,(hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,k+2)=hh;
    k = k+3;
end
Hs(:,:,1)=ll;

max_d = 0;
% align wavelet filters
for m=3*M-1:3*M+1
    hd=abs(Hd(:,:,m));

    id = (sum(sum(hd.*repmat([1:size(Hd,1)]',[1,size(Hd,2)])))/sum(hd(:)));
    jd = (sum(sum(hd.*repmat(1:size(Hd,2),[size(Hd,1),1])))/sum(hd(:)));
    max_d = max(max_d,round(abs((size(Hd,1)+1)/2-id)));
    max_d = max(max_d,round(abs((size(Hd,1)+1)/2-jd)));
   
    for n=m-3:-3:2
        hd=abs(Hd(:,:,n));
        Id = (sum(sum(hd.*repmat([1:size(Hd,1)]',[1,size(Hd,2)])))/sum(hd(:)));
        Jd = (sum(sum(hd.*repmat(1:size(Hd,2),[size(Hd,1),1])))/sum(hd(:)));

        kd = round([id;jd]-[Id;Jd]);    
        I = zeros(size(Hd(:,:,1)));
        I((end+1)/2+kd(1),(end+1)/2+kd(2))=1;
        Hd(:,:,n) = conv2(Hd(:,:,n),I,'same');
        Hs(:,:,n) = conv2(Hs(:,:,n),I(end:-1:1,end:-1:1),'same');
    end    
end
max_d = 2*round(max_d);

% center wavelet filters
Fd = zeros([size(Hd,1)+max_d,size(Hd,2)+max_d,size(Hd,3)]);
Fs = zeros([size(Hd,1)+max_d,size(Hd,2)+max_d,size(Hd,3)]);

Id = (size(Hd,1)+1)/2;
Jd = (size(Hd,2)+1)/2;
for m=3*M-1:3*M+1
    hd=abs(Hd(:,:,m));

    id = (sum(sum(hd.*repmat([1:size(Hd,1)]',[1,size(Hd,2)])))/sum(hd(:)));
    jd = (sum(sum(hd.*repmat(1:size(Hd,2),[size(Hd,1),1])))/sum(hd(:)));
   
    kd = round([Id;Jd]-[id;jd]);    
    Fd( kd(1)+(end+1)/2+[-(Id-1):(Id-1)], kd(2)+(end+1)/2+[-(Jd-1):(Jd-1)],m:-3:2) = Hd(:,:,m:-3:2);
    Fs(-kd(1)+(end+1)/2+[-(Id-1):(Id-1)],-kd(2)+(end+1)/2+[-(Jd-1):(Jd-1)],m:-3:2) = Hs(:,:,m:-3:2);
end

m=1;
hd=abs(Hd(:,:,m));

id = (sum(sum(hd.*repmat([1:size(Hd,1)]',[1,size(Hd,2)])))/sum(hd(:)));
jd = (sum(sum(hd.*repmat(1:size(Hd,2),[size(Hd,1),1])))/sum(hd(:)));

kd = round([Id;Jd]-[id;jd]);    
Fd( kd(1)+(end+1)/2+[-(Id-1):(Id-1)], kd(2)+(end+1)/2+[-(Jd-1):(Jd-1)],1) = Hd(:,:,1);
Fs(-kd(1)+(end+1)/2+[-(Id-1):(Id-1)],-kd(2)+(end+1)/2+[-(Jd-1):(Jd-1)],1) = Hs(:,:,1);

%%%%%%%%%%%%%%%%%%

% compensate for overcomplete wavelet redundancies
Fs(:,:,1)=Fs(:,:,1)/(2^(2*M));

for m=1:M
    Fs(:,:,3*m-1)=Fs(:,:,3*m-1)/(2^(2*m));
    Fs(:,:,3*m  )=Fs(:,:,3*m  )/(2^(2*m));
    Fs(:,:,3*m+1)=Fs(:,:,3*m+1)/(2^(2*m));
end

% trim zeros in the border
i  = min(sum(abs(Fd+Fs),3)==0,[],2);
j  = min(sum(abs(Fd+Fs),3)==0,[],1);

bi=0;
for a=1:length(i)
    if i(a)==1
        bi=bi+1;
    else
        break;
    end    
end

bj=0;
for a=1:length(j)
    if j(a)==1
        bj=bj+1;
    else
        break;
    end    
end

Bi=0;
for a=length(i):-1:1
    if i(a)==1
        Bi=Bi+1;
    else
        break;
    end    
end

Bj=0;
for a=length(j):-1:1
    if j(a)==1
        Bj=Bj+1;
    else
        break;
    end    
end

ki = min(bi,Bi);
kj = min(bj,Bj);

Fd = Fd(1+ki:end-ki,1+kj:end-kj,:);
Fs = Fs(1+ki:end-ki,1+kj:end-kj,:);

% reorder to go from coarse to fine
Fd(:,:,2:end)=Fd(:,:,end:-1:2);
Fs(:,:,2:end)=Fs(:,:,end:-1:2);
    
function H = up(h)
H=zeros(size(h)*2-[1 1]);
H(1:2:end,1:2:end)=h;

