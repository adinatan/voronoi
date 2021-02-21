%% make synthetic data
clear
N=1e5; % # of events
x=rand(1,N);
y=x+randn(1,N)/10 ;
d=1/min([min(abs(diff(y(:)))) min(abs(diff(x(:))))]);

h=hist3([x(:) y(:)],[256,256]);

%%
S = qtdecomp(h,15); % sparse
blocks = zeros(size(S));
pow = 2.^(8:-1:0);
bt=[];
for n=1:numel(pow)
    numblocks = numel(find(S==pow(n)));
    if (numblocks > 0)
        values = ones(pow(n),pow(n),numblocks);%repmat(uint8(1),[dim dim numblocks]);
        values(2:pow(n),2:pow(n),:) = 0;
        blocks = qtsetblk(blocks,S,pow(n),values);
        
        [binvals,r,c] = qtgetblk(h,S,pow(n));
        btree{n}=squeeze(sum(sum(binvals)));
        bt=[bt ; squeeze(sum(sum(binvals)))];
    end
end

blocks(end,1:end) = 1;
blocks(1:end,end) = 1;
imagesc(blocks+h./max(h(:)))

%%
% pow= 2.^(8:-1:0)
% for n=1:numel(pow)
%
%      [vals,r,c] = qtgetblk(h,S,pow(n));
%      for n2=1:size(vals,3)
%          b{n}=sum(sum(vals(:,:,n2)))
%
%

%%


