%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short g;
C_Target = (62.75+ 68.25)/2;
S = 2228.072;           
K = 2260;               
r = 0.0427;             
p = 10; 
toleranslevel = 0.01;
holiday_vector = ["2023-12-24";"2023-12-25";"2023-12-26";"2023-12-31";"2024-01-01";"2024-01-05";
    "2024-03-29";"2024-04-01";"2024-05-09";"2024-06-06";"2024-06-21"];
startDate = datetime(2023,11,23);
endDate = datetime(2024,03,04);

t = days252bus(startDate, endDate, holiday_vector);
deltaTVar = calculateDeltaT(t, p, 252);

disp("--------")
disp("Task 2a):")
disp("--------")
sigma = ImplicitVolatility(C_Target, S, K, r, t/252, toleranslevel);
d = calculateD(sigma, deltaTVar);
u = calculateU(sigma, deltaTVar);
q = calculateQ(r, d, u, deltaTVar);

disp("--------")
disp("Task 2b):")
disp("--------")
indexPriceMatrix = calcIndexPrice(p, u, d, S);
optionPriceMatrix = zeros(p+1); %have to instantiate here, beacuse recursive algorithm
optionPriceMatrix = calcOptionPrice(optionPriceMatrix, p+1, indexPriceMatrix, K, q, exp(-1*r*deltaTVar));
disp("European binary tree (without dividend payments)");
disp(optionPriceMatrix);

disp("--------")
disp("Task 2c):")
disp("--------")
result = findOptimalPeriods(S,K, r, sigma,t, 200);
disp("Number of periods within 0,5% of BSM = " + result);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImpVol = ImplicitVolatility(C_Target, S, K, r, t, toleranslevel)
    sigma_left = 0;
    sigma_right = 1;
    iteration = 0;
    C = 1000000;
    while(abs(C_Target - C) > toleranslevel)
        C_left = BlackScholes(S,K,r, t, sigma_left);
        C_right = BlackScholes(S,K,r, t, sigma_right);
        diff_C_left = abs(C_Target - C_left);
        diff_C_right = abs(C_Target - C_right);
        if(diff_C_left < diff_C_right)
            sigma_right = (sigma_right + sigma_left)/2;
            C = C_left;
            sigma = sigma_left;
        else
            sigma_left = (sigma_right + sigma_left)/2;
            C = C_right;
            sigma = sigma_right;
        end
        disp("Iteration "+ iteration+ ": " + "Black-Scholes = " + C);
        disp("Iteration "+ iteration+ ": " + "left-side sigma = " + sigma_left*100 + "%");
        disp("Iteration "+ iteration+ ": " + "right-side sigma = " + sigma_right*100 + "%");
        iteration = iteration +1;
    end
    ImpVol = sigma;
    disp("Implicit Volatility = "+ ImpVol*100 + "%");
    disp("Midprice = " + C_Target + ", Derived Black-Scholes = " + C);
    disp("Difference between Midprice and ​Derived Black-Scholes = " + abs(C_Target - C)*100 + "%");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recursive function that iterates thorugh underlying matrix. First call to
% the function the first column is calculated with max otherwise expneg and
% q. Then when we stop the recusrive function when p == 1
function result = calcOptionPrice(optionPriceMatrix, p, indexPriceMatrix, K, q, expNeg)
    %Iterate thorugh column j, recursive is iterate thorugh p
    for j = 1:p
        if(p == length(optionPriceMatrix)) 
            optionPriceMatrix(j, p) = max((indexPriceMatrix(j,p)-K),0);
        else
            optionPriceMatrix(j, p) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
        end
    end
    if(p == 1) 
         optionPriceMatrix(1, 1) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
         result = optionPriceMatrix;
    else
        result = calcOptionPrice(optionPriceMatrix, p-1, indexPriceMatrix, K, q, expNeg);
    end
    
end

function indexPriceMatrix = calcIndexPrice(p, u, d, S)
    indexPriceMatrix = zeros(p+1);
    for i = 0:p %rows
        for j = 0:i %columns
            indexPriceMatrix(j+1, i+1) = S*u^(i-j)*d^j;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = findOptimalPeriods(S,K, r, sigma,t, p)
    vectorPrice = zeros(p);
    counter = 0;
    C=BlackScholes(S, K, r, t/252, sigma);
    for p = 1:p
        deltaTVar = calculateDeltaT(t, p, 252);
        d = calculateD(sigma, deltaTVar);
        u = calculateU(sigma, deltaTVar);
        q = calculateQ(r, d, u, deltaTVar);
        indexPriceMatrix = calcIndexPrice(p, u, d, S);
        optionPriceMatrix = zeros(p+1); %have to instantiate here, beacuse recursive algorithm
        optionPriceMatrix = calcOptionPrice(optionPriceMatrix, p+1, indexPriceMatrix, K, q, exp(-1*r*deltaTVar));
        if(abs(C-optionPriceMatrix(1,1))/C < 0.005)
            vectorPrice(p,p) = optionPriceMatrix(1,1);
            counter = counter + 1;
        end
    end
    result = counter;
    x = 1:length(vectorPrice);
    y = vectorPrice;

    figure;
    plot(x, y); 
    grid on;
    xlabel('Iteration');
    ylabel('Function Output');
    title('Function Output Over Iterations');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Help functions
function deltaT = calculateDeltaT(t, p, daysInAYear)
    deltaT = (t/p)/daysInAYear;
end
function C = BlackScholes(S, K, r, t, sigma)
        d1 = calcD1(S, K, r, sigma, t);
        d2 = calcD2(d1,sigma, t);
        C = calcCallOption(S, K, t, d1, r, d2);
end

function d1 = calcD1(S,K, r, sigma, t)
    d1 = (log(S/K)+(r+sigma^2/2)*(t))/(sigma*sqrt(t));
end
function d2 = calcD2(d1, sigma, t)
    d2 = d1 - (sigma*sqrt(t));
end
function C = calcCallOption(S, K, t, d1, r, d2)
    C = S*normcdf(d1) - K*exp(-1*r*t)*normcdf(d2);
end

function u = calculateU(sigma, deltaT)
    u = exp(sigma*sqrt(deltaT));
end

function d = calculateD(sigma, deltaT)
    d = exp((-1)*sigma*sqrt(deltaT));
end

function q = calculateQ(r, d, u, t)
    q = (exp(r*t)-d)/(u-d);
end


