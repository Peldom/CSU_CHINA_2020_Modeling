clear;
clc;

%% ?????
x0 = [1, 1, 1, 1];
tspan = [0, 100];
[T, X] = ode45(@f, tspan, x0);

%% ??
figure;
subplot(4, 1, 1)
plot(T, X(:, 1))
ylabel('X_1');

subplot(4, 1, 2)
plot(T, X(:, 2))
ylabel('X_2');

subplot(4, 1, 3)
plot(T, X(:, 3))
ylabel('X_3');

subplot(4, 1, 4)
plot(T, X(:, 4))
ylabel('X_4');
xlabel('Time');

%% ??????
function fx = f(t, x) %# ok

    % ???fx???????
    fx = zeros(4, 1);

    % ???????
    fx(1) = -1.8637 * x(1) + 1.9749 * x(2) - 0.3442 * x(3) + 0.2657 * x(4);
    fx(2) = -2.1433 * x(2) + 1.8551 * x(1) + 0.4721 * x(3) - 0.1260 * x(4);
    fx(3) = -1.6595 * x(3) - 2.3399 * x(1) + 3.4431 * x(2) + 0.6238 * x(4);
    fx(4) = -2.5938 * x(4) + 4.5728 * x(1) - 3.1455 * x(2) + 1.0323 * x(3);
    
end