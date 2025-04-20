%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short g;
Market_Price = (34.75+54.75)/2;
S = 569;           
<<<<<<< HEAD
K = 630;               % Strike price 
sigma = 0.3183; %borde vara 0.3083
r = 0.0427;             % Risk-free rate
p = 200; %8 är optimalt
holiday_vector = ["2023-12-24";"2023-12-25";"2023-12-26";"2023-12-31";"2024-01-01";"2024-01-05";
    "2024-03-29";"2024-04-01";"2024-05-09";"2024-06-06";"2024-06-21"];

=======
K = 630;             
sigma = 0.3083; %Kanske är närmare 0.3183
r = 0.0427;     
DIV = 5.3;
p = 200; % p = 8 är optimalt
holiday_vector = ["2023-12-24";"2023-12-25";"2023-12-26";"2023-12-31";"2024-01-01";"2024-01-05";
    "2024-03-29";"2024-04-01";"2024-05-09";"2024-06-06";"2024-06-21"];
>>>>>>> f6a5cca1ef99f0bf34efe450de923afa6ff9fe78
startDate = datetime(2023,11,23);
divDate = datetime(2024,04,12);
endDate = datetime(2024,09,04);

t = days252bus(startDate, endDate, holiday_vector);
t_div = days252bus(startDate, divDate, holiday_vector)/252;
deltaTVar = calculateDeltaT(t, p, 252);
d = calculateD(sigma, deltaTVar);
u = calculateU(sigma, deltaTVar);
q = calculateQ(r, d, u, deltaTVar);
PDIV = calcPDiv(t_div, r, DIV);
S_star = calcSStar(PDIV, S);
divPeriod = calcDivPeriod(t, (t_div*252), p);
            
stockPriceMatrix = calcStockPrice(p,u,d,S_star);
americanoOptionPriceMatrix = zeros(p+1); 
americanoOptionPriceMatrix = calcOptionPrice(americanoOptionPriceMatrix, p+1, stockPriceMatrix, K, q, exp(-1*r*deltaTVar), DIV, S_star, divPeriod);
disp("Market Price: " + Market_Price);
disp("Value: "+americanoOptionPriceMatrix(1,1));
disp("Quota between Market and Value: " + (abs(Market_Price - americanoOptionPriceMatrix(1,1))/Market_Price)*100 + "%");
disp("Using implicit volatility = " + sigma*100 + "% and p = " + p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%findOptimalPeriods(S,K, r, sigma,t, 1000, Market_Price);
%Help functions
function findOptimalPeriods(S,K, r, sigma,t, p, Market_Price)
    holiday_vector = [
        "2023-12-25";
        "2023-12-26";
        "2024-01-01";
        "2024-01-06";
        "2024-03-29";
        "2024-03-31";
        "2024-04-01";
        "2024-05-01";
        "2024-05-16";
        "2024-06-06";
        "2024-06-21";
        "2024-06-22";
    ];
    vectorPrice = zeros(p);
    startDate = datetime(2023,11,23);
    divDate = datetime(2024,04,12);
    t_div = days252bus(startDate, divDate, holiday_vector)/252;
    for p = 1:p
        deltaTVar = calculateDeltaT(t, p, 252);
        d = calculateD(sigma, deltaTVar);
        u = calculateU(sigma, deltaTVar);
        q = calculateQ(r, d, u, deltaTVar);
        DIV = 5.3;
        PDIV = calcPDiv(t_div, r, DIV);
        S_star = calcSStar(PDIV, S);
        divPeriod = calcDivPeriod(t, (t_div*252), p);
        stockPriceMatrix = calcStockPrice(p,u,d,S_star);
        americanoOptionPriceMatrix = zeros(p+1);
        americanoOptionPriceMatrix = calcOptionPrice(americanoOptionPriceMatrix, p+1, stockPriceMatrix, K, q, exp(-1*r*deltaTVar), DIV, S_star, divPeriod);
        disp("Market_Price: " + Market_Price);
        disp("optionPriceMatrix(1,1): " + americanoOptionPriceMatrix(1,1));
        a = Market_Price-americanoOptionPriceMatrix(1,1);
        disp("Market_Price-result: " + a);
        disp("abs(Market_Price-optionPriceMatrix(1,1)): " + abs(Market_Price-americanoOptionPriceMatrix(1,1)));
        disp("abs(Market_Price-optionPriceMatrix(1,1))/Market_Price: " + abs(Market_Price-americanoOptionPriceMatrix(1,1))/Market_Price);
        if(abs(Market_Price-americanoOptionPriceMatrix(1,1))/Market_Price < 0.005)
            vectorPrice(p,p) = americanoOptionPriceMatrix(1,1);
        end
    end
    x = 1:length(vectorPrice);
    y = vectorPrice;
    
    figure; 
    plot(x, y); 
    grid on; 
    xlabel('Iteration'); 
    ylabel('Function Output');
    title('Function Output Over Iterations'); 
end

function divPeriod = calcDivPeriod(t, t_div, p)
    divPeriod = floor((t_div/t)*p);
end
% recursive function that iterates thorugh underlying matrix. First call to
% the function the first column is calculated with max otherwise expneg and
% q. Then when we stop the recusrive function when p == 1
function result = calcOptionPrice(optionPriceMatrix, p, stockPriceMatrix, K, q, expNeg, DIV, S_star, divPeriod)
    %Iterate thorugh column j, recursive is iterate thorugh p
    for j = 1:p
        if(p == length(optionPriceMatrix)) 
            optionPriceMatrix(j, p) = max((stockPriceMatrix(j,p)-K),0);
        else
            if(p == divPeriod)
                optionPriceMatrix(j, p) = max(expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1)),stockPriceMatrix(j,p)+DIV -K);
            else
                optionPriceMatrix(j, p) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
            end
        end
    end
    if(p == 1) 
         optionPriceMatrix(1, 1) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
         result = optionPriceMatrix;
    else
        result = calcOptionPrice(optionPriceMatrix, p-1, stockPriceMatrix, K, q, expNeg, DIV, S_star, divPeriod);
    end
end

function stockPriceMatrix = calcStockPrice(p, u, d, S_star)
    stockPriceMatrix = zeros(p+1);
    for i = 0:p %rows
        for j = 0:i %columns
            stockPriceMatrix(j+1, i+1) = S_star*u^(i-j)*d^j;
        end
    end
end
function S_star = calcSStar(PDIV, S)
    S_star = S-PDIV;
end
function p_div = calcPDiv(t_div, r, DIV)
    p_div = DIV*exp(-1*r*t_div);
end

function deltaT = calculateDeltaT(t, p, daysInAYear)
    deltaT = (t/p)/daysInAYear;
end

function u = calculateU(sigma, deltaT)
    u = exp(sigma*sqrt(deltaT));
end

function d = calculateD(sigma, deltaT)
    d = exp((-1)*sigma*sqrt(deltaT));
end

function q = calculateQ(r, d, u, deltaT)
    q = (exp(r*deltaT)-d)/(u-d);
end