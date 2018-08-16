load faces.mat;
test_data=[test_faces' test_nonfaces']; 
train_y = [ones(2429,1)',zeros(2429,1)']; 
y = [ones(472,1)',zeros(472,1)']; 
test_y_pred = zeros(size(test_data,2),1)';
%standardizing training and testing data
n_train=size(training_data',1);
mu=mean(training_data',1);
xb=repmat(mu,n_train,1);
train_x=(training_data'-xb)./std(training_data');

n_test = size(test_data',1);
mu_test=mean(test_data',1);
xb_test=repmat(mu_test,n_test,1);
test_x=(test_data'-xb_test)./std(test_data');

% setting up initial weights for neural network
%weights0=-0.01*ones(361,4); 
weights0 = normrnd(0,0.00001,361,4);
%weights1=-0.01*ones(4,1);
weights1 = normrnd(0,0.00001,4,1);
% define learning rate, rho
rho=0.01;
% setting up feed forward network with back propagaion
% only running 10000 iterations to prevent vanishing gradient issue
for i = 1:10000
    %multiplying train_x by initial weights, then applyiing nonlinearity
    %via sigmoid function to get hidden layer
    hidden_layer = 1./(1+exp(-train_x*weights0));
    %multiplying hidden layer by weight1, then applyiing nonlinearity
    %via sigmoid function to get output layer
    output_layer = 1./(1+exp(-hidden_layer*weights1));
    % back propagation
    output_error=(train_y'-output_layer).*(output_layer.*(1-output_layer));
    hidden_error=output_error*weights1' .* (hidden_layer .* (1-hidden_layer));
    %updating weights
    weights1 = weights1 +rho*hidden_layer'*output_error;
    weights0 = weights0 +rho*train_x'*hidden_error;
end
% training error
train_error= 1-(sum(output_layer(1:2429)>=0.5)+sum(output_layer(2430:end)<0.5))/4858;
train_error
% 0.0068

% applying NN to test set
hidden_layer_test = 1./(1+exp(-test_x*weights0));
output_layer_test = 1./(1+exp(-hidden_layer_test*weights1));

% testing error
test_error= 1-(sum(output_layer_test(1:472)>=0.5)+sum(output_layer_test(473:end)<0.5))/944;  
test_error
%0.2108