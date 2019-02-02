function Sol=MassActionGEM(GEM,x0,K,timespan,Mets2Plot,Rxns2Plot,NetFlux)

%This function solves the dynamic system encoded by Sv=b under mass-action
%kinetics
%
%Inputs:
% GEM: COBRA-like structure containing the genome-scale model
%
% x0: vector of initial concentrations, if empty, values are taken randomly
%
% K: vector of k values (reaction rate constants) for each reaction.
% Reversible reactions are considered as two independent reactions, the
% reverse sense is adjoined at the end of the reaction list in the GEM
%
% timespan: [tinit,tfinal] indicating the time span of the simulation
%
% Mets2Plot, Rxns2Plot: an array with the indeces in the GEM of the
% metabolites of the reactions, respectively, that should be in the plot.
% 
% NetFlux: either T or F, indicating if the flux values of reversible
% reactions should be given for the net contribution of both directions or
% separated per direction, respectively.
%
% Semidan, September 2016

if ~isfield(GEM,'rxnNames'),
    for i=1:size(GEM.S,2),
        Vnam{i,1}=['v',num2str(i)];
    end
    GEM.rxnNames=Vnam;
end
if ~isfield(GEM,'metNames'),
    for i=1:size(GEM.S,1),
        Xnam{i,1}=['x',num2str(i)];
    end
    GEM.metNames=Xnam;
end

S=[GEM.S,-GEM.S(:,GEM.rev==1)];
revNames=GEM.rxnNames(GEM.rev==1);
for i=1:length(revNames),
    revNames{i}=[revNames{i},'_[REV]'];
end
GEM.rxnNames=[GEM.rxnNames;revNames];
Rxns=size(S,2);
Mets=size(S,1);

if nargin<2 || isempty(x0),
    x0=10*rand(Mets,1);
    Sol.x0 = x0;
end
if nargin<3 || isempty(K),
     K=1*rand(Rxns,1);
     Sol.K = K;
end
if nargin<4 || isempty(timespan),
    timespan=[0,10];
    Sol.timespan = timespan;
end
if nargin<5,
    Mets2Plot=1:Mets;
end
if nargin<6,
    Rxns2Plot=[];
end
if nargin<7,
    NetFlux='T';
end


%obtain system of differential equations
v=zeros(Rxns,1);
dx = zeros(Mets,1);
function dx = f(x,t)
   for idx=1:Rxns,
       Subs=find(S(:,idx)<0);
       v(idx)=K(idx)*prod(x(Subs).^abs(S(Subs,idx)));
   end
   dx=S*v;
endfunction
tic
%solve system
[T,Y]=lsode("f",x0,timespan);
%[T,Y]=ode23s("f",x0,timespan);

Sol.T=T;
Sol.Y=Y;
solverTime=toc;
npoints=size(Y,1);

%recover flux values and derivatives
V=zeros(npoints,Rxns);dX=zeros(npoints,Mets);
V2Rec=1:Rxns;
for i=1:npoints,
    for j=1:length(V2Rec),
       Subs=find(S(:,V2Rec(j))<0);
       V(i,V2Rec(j))=K(V2Rec(j))*prod(Y(i,Subs)'.^abs(S(Subs,V2Rec(j))));
    end
    dX(i,:)=V(i,:)*S';
end

if NetFlux=='T',
    V(:,GEM.rev==1)=V(:,GEM.rev==1)-V(:,(1:sum(GEM.rev))+size(GEM.S,2));
    V(:,(1:sum(GEM.rev))+size(GEM.S,2))=[];
end

%plot figures
if ~isempty(Mets2Plot)
    figure(1)
    plot(T,Y(:,Mets2Plot))
    legend(GEM.metNames(Mets2Plot))
    xlabel('time')
    ylabel('concentration')
end
if ~isempty(Rxns2Plot)
    figure(2)
    plot(T,V(:,Rxns2Plot))
    legend(GEM.rxnNames(Rxns2Plot))
    xlabel('time')
    ylabel('flux value')
end



Sol.V=V;
Sol.dX=dX;
Sol.solverTime=solverTime;

end