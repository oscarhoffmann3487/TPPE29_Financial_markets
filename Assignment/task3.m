%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short g;
S = 62;          
K = 55;          
sigma = 0.2;
r = 0.0427;  
DIV = 2.7;
p = 8; 
holiday_vector = ["2023-12-24";"2023-12-25";"2023-12-26";"2023-12-31";"2024-01-01";"2024-01-05";
    "2024-03-29";"2024-04-01";"2024-05-09";"2024-06-06";"2024-06-21"];
<<<<<<< HEAD

=======
>>>>>>> f6a5cca1ef99f0bf34efe450de923afa6ff9fe78
startDate = datetime(2023,11,23);
divDate = datetime(2024,03,21);
endDate = datetime(2024,07,19);

t = days252bus(startDate, endDate, holiday_vector);
t_div = days252bus(startDate, divDate, holiday_vector)/252;
deltaTVar = calculateDeltaT(t, p, 252);
d = calculateD(sigma, deltaTVar);
u = calculateU(sigma, deltaTVar);
q = calculateQ(r, d, u, deltaTVar);
PDIV = calcPDiv(t_div, r, DIV);
S_star = calcSStar(PDIV, S);
divPeriod = calcDivPeriod(t, t_div*252, p);

stockPriceMatrix = calcStockPrice(p,u,d,S_star);
euroOptionPriceMatrix = zeros(p+1);
calcOptionPrice(euroOptionPriceMatrix, p+1, stockPriceMatrix, K, q, exp(-1*r*deltaTVar), false, DIV, S_star, divPeriod);
americanoOptionPriceMatrix = zeros(p+1);
calcOptionPrice(americanoOptionPriceMatrix, p+1, stockPriceMatrix, K, q, exp(-1*r*deltaTVar), true, DIV, S_star, divPeriod);

function divPeriod = calcDivPeriod(t, t_div, p)
    divPeriod = floor((t_div/t)*p);
end
% recursive function that iterates thorugh underlying matrix. First call to
% the function the first column is calculated with max otherwise expneg and
% q. Then when we stop the recusrive function when p == 1
function calcOptionPrice(optionPriceMatrix, p, stockPriceMatrix, K, q, expNeg, americano, DIV, S_star, divPeriod)
    %Iterate thorugh column j, recursive is iterate thorugh p
    for j = 1:p
        if(p == length(optionPriceMatrix)) 
            optionPriceMatrix(j, p) = max((stockPriceMatrix(j,p)-K),0);
        else
            if(americano && p == divPeriod)
                optionPriceMatrix(j, p) = max(expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1)),stockPriceMatrix(j,p)+DIV -K);
            else
                optionPriceMatrix(j, p) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
            end
        end
    end
    if(p == 1) 
         optionPriceMatrix(1, 1) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
         if(americano)
             disp("American binary tree (with dividend payments):");
         else
             disp("European binary tree (with dividend payments):");
         end
         disp(optionPriceMatrix);
    else
        calcOptionPrice(optionPriceMatrix, p-1, stockPriceMatrix, K, q, expNeg, americano, DIV, S_star, divPeriod)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Help functions
function stockPriceMatrix = calcStockPrice(p, u, d, S_star)
    stockPriceMatrix = zeros(p+1);
    for i = 0:p %rows
        for j = 0:i %columns
            stockPriceMatrix(j+1, i+1) = S_star*u^(i-j)*d^j;
        end
    end
    disp("Underlying binary tree:" )
    disp(stockPriceMatrix);
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