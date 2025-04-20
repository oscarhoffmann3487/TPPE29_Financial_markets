%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short g; % for displaying purposes
S0 = 62;        % Initial stock price
K = 55;         % Strike price (lösenpris)
T = 1;          % Time to maturity (in years))
r = 0.0427;     % Risk-free rate
sigma = 0.20;   % implicit volatility        
p = 8;          % Number of periods in the binomial tree
startDate = datetime(2023,11,23); % Valuation date
endDate = datetime(2024,07,19); % Exloration day
daysInAYear = 252;
holiday_vector = ["2023-12-24";"2023-12-25";"2023-12-26";"2023-12-31";"2024-01-01";"2024-01-05";
    "2024-03-29";"2024-04-01";"2024-05-09";"2024-06-06";"2024-06-21"];

t = days252bus(startDate, endDate, holiday_vector);
deltaTVar = calculateDeltaT(t, p, daysInAYear);
expPos = calculateExpPositive(r, deltaTVar);
expNeg = calculateExpNegative(r, deltaTVar);
d = calculateD(sigma, deltaTVar);
u = calculateU(sigma, deltaTVar);
q = calculateQ(expPos, d, u);
disp("values to check")
disp("u: " + u);
disp("d: " + d);
disp("q: " + q);
disp("deltaT: " + deltaTVar);
disp("expPos: " + expPos);
disp("expNeg: " + expNeg);

%Underlying
stockPriceMatrix = calcStockPrice(p, u, d, S0);
<<<<<<< HEAD
formattedStockPriceMatrix = formatBinomialTree(stockPriceMatrix); 
disp(formattedStockPriceMatrix); 
%disp(stockPriceMatrix);
=======
disp("Underlying binary tree: ");
disp(stockPriceMatrix); 
>>>>>>> f6a5cca1ef99f0bf34efe450de923afa6ff9fe78

%Option
optionPriceMatrix = zeros(p+1); %have to instantiate here, beacuse recursive algorithm
optionPriceMatrix = calcOptionPrice(optionPriceMatrix, p+1, stockPriceMatrix, K, q, expNeg);
<<<<<<< HEAD
formattedStockPriceMatrix = formatBinomialTree(optionPriceMatrix); 
disp(formattedStockPriceMatrix); 
%disp(optionPriceMatrix);
=======
disp("European binary tree (without dividend payments)");
disp(optionPriceMatrix); 
>>>>>>> f6a5cca1ef99f0bf34efe450de923afa6ff9fe78

function stockPriceMatrix = calcStockPrice(p, u, d, S0)
    stockPriceMatrix = zeros(p+1);
    for i = 0:p %rows
        for j = 0:i %columns
            stockPriceMatrix(j+1, i+1) = S0*u^(i-j)*d^j;
        end
    end
end

% recursive function that iterates thorugh underlying matrix. First call to
% the function the first column is calculated with max otherwise expneg and
% q. Then when we stop the recusrive function when p == 1
function result = calcOptionPrice(optionPriceMatrix, p, stockPriceMatrix, K, q, expNeg)
    %Iterate thorugh column j, recursive is iterate thorugh p
    for j = 1:p
        if(p == length(optionPriceMatrix)) 
            optionPriceMatrix(j, p) = max((stockPriceMatrix(j,p)-K),0);
        else
            optionPriceMatrix(j, p) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
        end
    end
    if(p == 1) 
         optionPriceMatrix(1, 1) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
         result = optionPriceMatrix;
    else
        result = calcOptionPrice(optionPriceMatrix, p-1, stockPriceMatrix, K, q, expNeg);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Help functions
function deltaT = calculateDeltaT(t, p, daysInAYear)
    deltaT = (t/p)/daysInAYear;
end

function expPositive = calculateExpPositive(r, deltaTVar) 
    expPositive = exp(r*deltaTVar);
end

function expNegative = calculateExpNegative(r, deltaTVar) 
    expNegative = exp((-1)*r*deltaTVar);
end

function u = calculateU(sigma, deltaT)
    u = exp(sigma*sqrt(deltaT));
end

function d = calculateD(sigma, deltaT)
    d = exp((-1)*sigma*sqrt(deltaT));
end

function q = calculateQ(expPos, d, u)
    q = (expPos-d)/(u-d);
end

