%  Marine Predators Algorithm source code (Developed in MATLAB R2019b)
% --------------------------------------------
% fobj = @YourCostFunction
% dim = number of your variables
% maxFEs = maximum number of fitness evaluations
% N = number of solutions
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% ---------------------------------------------------------

clear all
clc
%rng(10)   
format long
N=30; % Number of solutions in each iteration

F_index=1;% F1:F29--> 1:29
Function_name=strcat('F',num2str(F_index));

DisplayResults=0;
global maxFEs bestStar 
[lb,ub,dim,bestStar]=test_functions_range(F_index);
fobj=@test_functions;

maxFEs=dim*10000; % Maximum number of fitness evaluations
if F_index<=7
    Parameter.Smax=2;
else 
    Parameter.Smax=40;
end
Parameter.nPop=N;
Parameter.sigma_initial2=0.5;
Parameter.sigma_final2=0.1;
Parameter.FESstar=maxFEs;
Parameter.MaxIt=round(Parameter.FESstar/Parameter.nPop);

if F_index>=24 &&F_index<=29 %composite functions
    global initial_flag
    initial_flag=0;
     clear fun_num func o sigma lamda bias M 
end


[Convergence_curve,FEs_counts,X ,Fbest]=ICO(F_index,N,fobj,lb, ub , dim,Parameter,DisplayResults);

display(['The best solution obtained by ICO is : ', num2str(X,10)]);
display(['The best optimal value of the objective function found by ICO is : ', num2str(Fbest,10)]);
disp(sprintf('--------------------------------------'));

% Convergence curve
semilogy(FEs_counts,Convergence_curve,'Color','r')
title('Objective space')
xlabel('FEs');
ylabel('Best score obtained so far');

