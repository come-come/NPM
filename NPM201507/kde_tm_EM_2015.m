function [Idx] = kde_tm_EM_2015(D,K,Par)
%[Idx,Pf,Pf_d,logPf_d,Li,Kr,alpha] = kde_tm_EM_2015(D,K,kk,Par)

if ~isfield(Par,'anchor')
    Par.anchor  = 100;
end
if ~isfield(Par,'maxit')
    Par.maxit  = 400;
end
if ~isfield(Par,'Leps')
    Par.Leps   = 1;
end

kk = Par.anchor;

% if nargin<4
%    Par.maxit  = 400;
%    Par.Leps   = 1;
% end
   Par.landmarker = 1;
   Par.bandwidth = 1;

dim = size(D{1},1);
allx = [];
for i = 1:length(D)
    allx = [allx,D{i}];
end
allx = allx';

n_allx = size(allx,1);
n_d = length(D); % # of documents


if Par.landmarker == 1
    % ---------- set anchor points - sampling
    N = kk;
    anchors = allx(randsamp(n_allx,N),:);
else
    % ---------- set anchor points - grid
    mdist=std(allx);
    xmax=max(allx)+mdist;
    xmin=min(allx)-mdist;
    if dim == 2
        x1=linspace(xmin(1),xmax(1),kk);
        x2=linspace(xmin(2),xmax(2),kk);
        [px,py]=meshgrid(x1,x2);
        anchors=[px(:) py(:)];
    elseif dim == 3
        x1=linspace(xmin(1),xmax(1),kk);
        x2=linspace(xmin(2),xmax(2),kk);
        x3=linspace(xmin(3),xmax(3),kk);
        [px,py,pz]=meshgrid(x1,x2,x3);
        anchors=[px(:) py(:) pz(:)];
    else
        disp('dimension error');
        return;
    end
end

% figure;
% plot(allx(:,1),allx(:,2),'r.');
% hold on;
% plot(anchors(:,1),anchors(:,2),'b.');

M = size(anchors,1);

% set bandwidth
% rule of thumb bandwidth suggested by Bowman and Azzalini (1997) p.31
%h = norm(median(abs(anchors-repmat(median(anchors),M,1)))/0.6745*(1/M)^(1/6));

h = Par.bandwidth*(1.06*mean(pdist(anchors))/M^(1/5));

% calculate kernel function
kerf=@(t)2*exp(-t/2)/(sqrt(2*pi)*h);
% kerf=@(t)2*exp(-t.*t/2)/(sqrt(2*pi)*h);
for t = 1:n_d
    xx = D{t}';
    mt = size(xx,1);
    Kr{t} = zeros(mt,M);
    for k = 1:M
    z = sum((xx-anchors(k+zeros(mt,1),:)).*(xx-anchors(k+zeros(mt,1),:)),2)/(h*h);
    Kr{t}(:,k) = kerf(z)+eps;
    end
end

% initialize
Li = [];
maxit = Par.maxit;
[Pf,alpha] = kde_tm_init(M,K);
logPf_d = zeros(n_d,K);
for t = 1:n_d
    for j = 1:K
        logPdt_fj = sum( log( Kr{t}*alpha(:,j) ) );
        logPf_d(t,j) = log(Pf(j))+logPdt_fj;
    end
end

% EM algorithm
for it = 1:maxit
    fprintf('Iteration %d ',it);

    % E-step
    [Pf_d] = kde_tm_Estep(D,Pf,logPf_d);

    % M-step
    [Pf,alpha,logPd_f,logPf_d] = kde_tm_Mstep(D,Pf_d,alpha,Kr);
    %logPf_d

    % Evaluate data log-likelihood
    Li(it) = kde_tm_logL(D,Pf,logPf_d);

    dLi = 0;
    if it > 1
         dLi    = Li(it) - Li(it-1);
         if (dLi<0)
             disp('error!!!!!!');
         end
         if abs(dLi) < Par.Leps
             %figure;plot(Li);
             disp('over');
             break
         end
    end
    fprintf('dLi=%f \n',dLi);
end

[~,Idx] = max(Pf_d,[],2);

function [Pf,alpha] = kde_tm_init(M,K)

Pf = ones(1,K)/K;           % uniform prior on topics

alpha = rand(M,K);
C    = 1./sum(alpha,1);     % normalize to sum to 1
alpha = alpha * diag(C);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E step compute posterior on z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pf_d] = kde_tm_Estep(D,Pf,logPf_d)

K = length(Pf);
n_d = length(D);
Pf_d = zeros(n_d,K);

% for t = 1:n_d
%     for j = 1:K
%         for i = 1:size(D{t},2)
%             logPxi_f(i) = log(Kr{t}(i,:)*alpha(:,j));
%         end
%         logPf_dt(j) = log(Pf(j))+sum(logPxi_f);
%     end
%     for j = 1:K
%         Pf_d(t,j) = 1/sum(exp(logPf_dt-logPf_dt(j)));
%     end
% end

for t = 1:n_d
    for j = 1:K
        Pf_d(t,j) = 1/(sum(  exp( logPf_d(t,:)-logPf_d(t,j) ) ) ) +eps;
    end
end

C = sum(Pf_d,2);
Pf_d = diag(1./C) * Pf_d;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M step, maximazize log-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pf,alpha,logPd_f,logPf_d] = kde_tm_Mstep(D,Pf_d,old_alpha,Kr)

n_d = length(D);
K = size(Pf_d,2);
M = size(Kr{1},2);

Pf = sum(Pf_d,1)/n_d;

% for j = 1:K
%     for k = 1:M
%         for t = 1:n_d
%             for i = 1:length(D{t})
%                 nu = Kr{t}(i,k);
%                 de = Kr{t}(i,:)*old_alpha(:,j);
%                 frac(i) = nu/de;
%             end
%             ss(t) = Pf_d(t,j)*sum(frac);
%         end
%         alpha(k,j) = sum(ss);
%     end
% end

for j = 1:K
    ss = [];
    for t = 1:n_d
        mt = size(Kr{t},1);
        nu = repmat(old_alpha(:,j)',mt,1).*Kr{t};
        de = repmat((Kr{t}(:,:)*old_alpha(:,j)+eps),1,M);
        ss = [ss,sum(nu./de)'];
    end
    alpha(:,j) = ss*Pf_d(:,j);
end

% for j = 1:K
%     for k = 1:M
%         ss = [];
%         for t = 1:n_d
%             nu = old_alpha(k,j)*Kr{t}(:,k);
%             de = Kr{t}(:,:)*old_alpha(:,j)+eps;
%             %frac = nu./de+dir-1;
%             frac = nu./de;
%             ss(t) = Pf_d(t,j)*sum(frac);
%         end
%         alpha(k,j) = sum(ss);
%     end
% end

C = 1./sum(alpha,1);     % normalize to sum to 1
alpha = alpha * diag(C);

logPf_d = zeros(n_d,K);
logPd_f = zeros(n_d,K);
for t = 1:n_d
    for j = 1:K
        logPd_f(t,j) = sum( log( Kr{t}*alpha(:,j) +eps));
        %logPd_f(t,j) = sum( log( Kr{t}*alpha(:,j) +eps) + (dir-1)*sum(log(alpha(:,j)+eps) ));
        logPf_d(t,j) = log(Pf(j))+logPd_f(t,j);
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data log-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = kde_tm_logL(D,Pf,logPf_d)

% n_d = length(D);
% maxlog = max(max(logPf_d));
% 
% for t = 1:n_d
%     L(t) = sum(exp( logPf_d(t,:)-repmat(maxlog,size(logPf_d(t,:)))));
% end
% 
% L = sum(log(L)+maxlog);

maxlog = max(logPf_d,[],2);
L = sum(  log( sum(exp(logPf_d-repmat(maxlog,1,size(logPf_d,2))),2) ) +maxlog  );

% n_d = length(D);
% 
% for t = 1:n_d
%     L(t) = Pf*logPd_f(t,:)';
% end
% 
% L = sum(L);
% 
%  L = sum(log(sum(exp(logPf_d),2)));

return;
