function  [Convergence_curve FEs_counts BestR Fbest fitness_history position_history  Trajectories it]=ICO(F_index,nPop,CostFunction,VarMin, VarMax , nVar,Parameter,DisplayResults,fpt)%
global    maxFEs bestStar

%Constants
Smin = 0;       % Minimum Number of colones
beta=100;       %Constant in (Eq. 13)
Exponent = 2;   %Variance Reduction Exponent (Eq. 18) 
mu=4;           %constant in (eq. 1)
gamma=1e-19;    %Constant in (eq. 15)
%Parameters
Smax=Parameter.Smax;
try
    sigma_initial2 =Parameter.sigma_initial2;
catch
    sigma_initial2=0.5;      % Initial Value of Standard Deviation
end
try
    sigma_final2 =Parameter.sigma_final2;
catch
    sigma_final2=0.1;      % Final Value of Standard Deviation
end

%%Start

mK=0.5*mean([(maxFEs/(Parameter.Smax*nPop)),(maxFEs/(1*nPop))]); %(eq. 9)

NumIt=round(maxFEs/Parameter.nPop);
fitness_history=zeros(nPop,NumIt);
position_history=zeros(nPop,NumIt,nVar);
Trajectories=zeros(nPop,NumIt);
Convergence_curve=zeros(1,NumIt);



mFES=0;
M=(VarMax(1)-VarMin(1))/2; % (Eq. 14)

dim=nVar;
VarSize = [1 nVar]; 


try
    fpt; 
catch
    fpt=0;
end

%% Initialization


%% Initialization
X(1,:) = unifrnd(0, 1, VarSize);
for i = 1:nPop
    % Initialize Solutions
    X(i+1,:) = mu* X(i,:).*(1- X(i,:)); %(Eq. 1)
end 
X(1,:)=[];
for i = 1:nPop
    X(i,:)=mapping(X(i,:),dim,VarMax,VarMin);
    % Evaluation
    Costs(i,1) = CostFunction(X(i,:),F_index,dim);
    mFES=mFES+1;
end

% Initialize Best Cost History
BestCostL=0;
[BestCost index]= min(Costs);
it=1;
Fbest=BestCost;
[Convergence_curve(it)] = BestCost;
fitness_history(:,it)=Costs;
if fpt==1
    position_history(:,it,:)=X;
    Trajectories(:,it)=X(:,1);
end
FEs_counts(it)=mFES;


%%ICO  Main Loop
acounter=0;
V=zeros(nPop,nVar);
it=2;
while mFES<=maxFEs
    
    Z=exp(-beta*it/mK);% (Eq. 13)
    if Z<=gamma
        beta=-log(10*gamma)*mK/it; %(Eq. 15)
        Z=exp(-beta*it/mK); %(Eq. 13)
        acounter=acounter+1;
    end
    
    % Update Standard Deviation
    sigma_iter= ((mK - it)/(mK - 1))^Exponent * (sigma_initial2 - sigma_final2) + sigma_final2; %(Eq. 18)
    alpha_iter=(10*(log(M)))*Z; % (Eq. 12)
    
    BestCost=min(Costs);
    WorstCost = max(Costs);
     
    
    if BestCost==WorstCost
        ratio=ones(1,nPop);
    else
        ratio=((Costs-WorstCost)./(BestCost-WorstCost))'; % (Eq. 4),        
    end
    
    A=TTEq(ratio,X,alpha_iter,it,mK);

    
    
    XL=[];VL=[];AL=[];Costs3=[];CostsG=[];
    XB=[];VB=[];AB=[];
    % Clone
    for i = 1:nPop
        
        S = floor(Smin + (Smax - Smin)*ratio(i)); %(Eq. 3)
        mFES=mFES+S; %Update the number of fitness evaluations
        for j = 1:S
            
            pm=rand;
            
            if pm<sigma_iter
                
                X2= X(i,:) +alpha_iter*randn(VarSize);% Second part of (Eq. 19)
                
                X2=space_bound(X2,VarMax,VarMin,1,nVar);
                V2=rand(1,dim).*V(i,:)+A(i,:);
                mCost=CostFunction(X2,F_index,dim);
                XL=[XL;X2];
                VL=[VL;V2];
                AL=[AL;A(i,:)];
                Costs3=[Costs3;mCost];
                
                
            else
                [X2 V2]=cloningTT(X(i,:),A(i,:),V(i,:),VarMin,VarMax,1,nVar); %function for the first part of (Eq. 19)
                mCost=CostFunction(X2,F_index,dim);
                XB=[XB;X2];
                VB=[VB;V2];
                AB=[AB;A(i,:)];
                CostsG=[CostsG;mCost];
                
            end
            
        end
        
    end
    [X,A,V,Costs]=OmitExtra(Costs,X,A,V);
    [XB,AB,VB,CostsG]=OmitExtra(CostsG,XB,AB,VB);
    [XL,AL,VL,Costs3]=OmitExtra(Costs3,XL,AL,VL);
    
    NL=round(min((sigma_iter)*nPop,size(XL,1))); % (Eq. 20)
    u=nPop-NL;%max(nPop-NL,1); % (Eq. 21)
    NB=min(round(u*90/100),size(XB,1));% (Eq. 22)
    NE=min(u-NB,size(X,1));% (Eq. 23)
    if NE==0
        if NB>0
            NB=NB-1;
        else
            NL=NL-1;
        end
        NE=1;
    end
    EX=[];ECosts=[];EV=[];EA=[];
    
    for i=1:(nPop-(NB+NE+NL))
        
        EX(i,:) = unifrnd(VarMin, VarMax, VarSize);
        ECosts(i,1) = CostFunction(EX(i,:),F_index,dim);
        mFES=mFES+1;
        EV(i,:)=zeros(1,nVar);
        EA(i,:)=EV(i,:);
    end
    X = [X([1:NE],:); XB([1:NB],:);XL([1:NL],:);EX];
    A = [A([1:NE],:); AB([1:NB],:);AL([1:NL],:);EA];
    V =  [V([1:NE],:);VB([1:NB],:);VL([1:NL],:);EV];
    Costs=[ Costs([1:NE]) ;CostsG([1:NB]);Costs3([1:NL]);ECosts];
    
    XL=[];VL=[];AL=[];Costs3=[];XB=[];AB=[];CostsG=[];VB=[];
    
    % Store Best Cost History
    [Convergence_curve(it) index] = min(Costs);
    
    
    fitness_history(:,it)=Costs;
    if fpt==1
        position_history(:,it,:)=X;
        Trajectories(:,it)=X(:,1);
    end
    FEs_counts(it)=mFES;
    
    Fbest=Convergence_curve(it);
    if DisplayResults==1
        % Display Iteration Information
        disp(['Iteration ' num2str(it-1) ': Best Cost = ' num2str(Convergence_curve(it))]);
    end
    if mFES>=maxFEs,break;end
    if Fbest==bestStar
        break
    end
    
    
    
    it=it+1;
    
end
if it<NumIt
    FEs_counts(it:NumIt)=mFES;
    Convergence_curve(it:NumIt)=Convergence_curve(it);
    for j=1:nPop
        fitness_history(j,it:NumIt)=fitness_history(j,it);
    end
end

%% Results

BestR=X(index,:);

L = length(Convergence_curve);
ind = round((1:NumIt)*L/NumIt);
Convergence_curve = Convergence_curve(ind);
FEs_counts=FEs_counts(ind);

L = size(Trajectories,2);
ind = round((1:NumIt)*L/NumIt);
Trajectories=Trajectories(:,ind);
position_history=position_history(:,ind,:);
fitness_history= fitness_history(:,ind);


end
function MX=mapping(X,dim,ub,lb)
if numel(ub)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
for i=1:dim
    ub_i=ub(i);
    lb_i=lb(i);
    MX(:,i)=X(i).*(ub_i-lb_i)+lb_i; %(Eq. 2)
end
end
function a=TTEq(NF,X,G,iteration,mk);

[N,dim]=size(X);
final_per=2; %In the last iteration, only 2 percent of solutions apply
y=round(N*(2+(1-iteration/mk)*98)/100); % (Eq. 7)
n_iter=max(y,1); %(Eq. 8)

[Ms ds]=sort(NF,'descend');
jj=ds(1:n_iter);
F1=sum(NF(jj)'.*X(jj,:),1)./n_iter; %(Eq. 17)

for i=1:N
    F=rand*F1; %(Eq. 17)
    R=norm(X(i,:)-F,2);% distance with Rnorm=2;
    E2(i,:)=((F-X(i,:)))/(R+eps); %(Eq. 16)
    
end
a=E2.*(20*G);%(Eq. 11)

end
function [X,V]=cloningTT(X,a,V,low,up,N,dim)

V=rand(1,dim).*V+a; %eq. 11. 
X=X+V; %eq. 12.
X=space_bound(X,up,low,N,dim);
end
%This function checks the search space boundaries for solutions.
function  X=space_bound(X,up,low,N,dim)

for i=1:N 
% Solutions that go out of the search space, are reinitialized randomly .
  
Tp=X(i,:)>up;Tm=X(i,:)<low;X(i,:)=(X(i,:).*(~(Tp+Tm)))+((rand(1,dim).*(up-low)+low).*(Tp+Tm));

end
end
function [X,A,V,Costs]=OmitExtra(Costs,X,A,V,meps)

[CG ic]=unique(Costs);
%[CG ic]=unique(X, 'rows');
X=X(ic,:);
V=V(ic,:);
A=A(ic,:);
Costs=Costs(ic,:);

end