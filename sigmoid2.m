function y = sigmoid(x,c,a)

%% Perform mathematics: 
y = 1./(1 + exp(-a.*(x-c)));
end
