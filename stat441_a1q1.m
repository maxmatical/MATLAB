load faces.mat;
% defining training and testing set
training_data=[train_faces' train_nonfaces']; 
test_data=[test_faces' test_nonfaces']; 
train_x = training_data;
n = size(train_x,2);
train_y = [ones(2429,1)',zeros(2429,1)']; 
test_x = test_data;
y = [ones(472,1)',zeros(472,1)']; 
test_y_pred = zeros(size(test_data,2),1)';
% initializing beta
% if beta set = 1, e^(beta*x) becomes too large, so define
% a much smaller value for beta
beta = ones(361,1)/50000;
% define convergence
convergence = 10^(-6);
% defining p_i
p = zeros(n, 1);
for i = 1:n
    p(i) = exp(beta'*train_x(:,i))/(exp(beta'*train_x(:,i)) + 1);
end;
% defining W
w_i = p.*(1-p);
W = diag(w_i);
train_x_t = train_x';
train_y_t = train_y';
% setting up the initial values of beta_t and beta_t+1
beta_current = beta;
beta_new = beta + inv(train_x_t'*W*train_x_t)*train_x_t'*(train_y_t-p);
% using newton-raphson method
while abs(beta_new - beta_current)>=convergence
    beta_t = beta_new;
    for i = 1:n
        p(i) = exp(beta_t'*train_x(:,i))/(1+exp(beta_t'*train_x(:,i)));
    end;
    w = p.*(1-p);
    W = diag(w);    
    beta_current = beta_t;
    beta_new = beta_t + inv(train_x_t'*W*train_x_t)*train_x_t'*(train_y_t-p);
end

% predicting against training set
train_y_pred = zeros(size(train_x,2),1)';
for i = 1:size(train_x,2)
    p = exp(beta_new'*train_x(:,i))/(exp(beta_new'*train_x(:,i)) + 1);
    if(p > 0.5)
        train_y_pred(i) = 1;
    end
end
train_error = size(find(abs(train_y-train_y_pred)),2)/size(train_y,2);
% predicting test set
for i = 1:size(test_x,2)
    p = exp(beta_new'*test_x(:,i))/(exp(beta_new'*test_x(:,i)) + 1);
    if(p > 0.5)
        test_y_pred(i) = 1;
    end
end
test_error = size(find(abs(y - test_y_pred)),2)/size(y,2);
train_error
test_error
