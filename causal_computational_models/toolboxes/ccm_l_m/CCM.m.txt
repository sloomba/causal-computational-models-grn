function [ SC ] = CCM( X , Y , tau , E , LMN , L )

% Plot convergence
% L - vector of computing points

l=length(L); 
SC=zeros(2,l);

for i = 1:l
    [SC(:,i)] = SugiLM(X(1:L(i)) , Y(1:L(i)) , tau , E , LMN );
end
plot(L,SC)
leg=legend('$\hat{X}(t)|M_Y$','$\hat{Y}(t)|M_X$');
set(leg,'Interpreter', 'latex');
xlabel('L','Interpreter', 'latex')
ylabel('$\rho$','Interpreter', 'latex')
end

