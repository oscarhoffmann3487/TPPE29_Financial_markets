%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short g;
S = 2228.072;           
K = 2260;               
r = 0.0427;             
p = 8; 
barrierLevel = 0.05;
sigma = 0.14667;
holiday_vector = ["2023-12-24";"2023-12-25";"2023-12-26";"2023-12-31";"2024-01-01";"2024-01-05";
    "2024-03-29";"2024-04-01";"2024-05-09";"2024-06-06";"2024-06-21"];
startDate = datetime(2023,11,23);
endDate = datetime(2024,03,04);

t = days252bus(startDate, endDate, holiday_vector);
deltaTVar = calculateDeltaT(t, p, 252);
d = calculateD(sigma, deltaTVar);
u = calculateU(sigma, deltaTVar);
q = calculateQ(r, d, u, deltaTVar);
barrierPrice = S*(1+barrierLevel);

indexPriceMatrix = calcIndexPrice(p, u, d, S);
plainVanillaOptionPriceMatrix = zeros(p+1); %have to instantiate here, beacuse recursive algorithm
plainVanillaOptionPriceMatrix = calcOptionPrice(plainVanillaOptionPriceMatrix, p+1, indexPriceMatrix, K, q, exp(-1*r*deltaTVar), nan);
upAndOutOptionPriceMatrix = zeros(p+1); %have to instantiate here, beacuse recursive algorithm
upAndOutOptionPriceMatrix = calcOptionPrice(upAndOutOptionPriceMatrix, p+1, indexPriceMatrix, K, q, exp(-1*r*deltaTVar), barrierPrice);

disp("Plain vanilla: " + plainVanillaOptionPriceMatrix(1,1));
disp("Up-down-and-out: " + upAndOutOptionPriceMatrix(1,1));
upAndIn = plainVanillaOptionPriceMatrix(1,1) - upAndOutOptionPriceMatrix(1,1);
disp("Up-down-and-in: " + upAndIn);
disp("Using barrier level = " + barrierLevel*100 + "%");

% recursive function that iterates thorugh underlying matrix. First call to
% the function the first column is calculated with max otherwise expneg and
% q. Then when we stop the recusrive function when p == 1
function result = calcOptionPrice(optionPriceMatrix, p, indexPriceMatrix, K, q, expNeg, barrierPrice)
    %Iterate thorugh column j, recursive is iterate thorugh p
    for j = 1:p
        if(p == length(optionPriceMatrix)) 
            optionPriceMatrix(j, p) = max((indexPriceMatrix(j,p)-K),0);
        else
            if(indexPriceMatrix(j,p) > barrierPrice) %Up-down-and-out condition
                optionPriceMatrix(j, p) = 0;
            else
                optionPriceMatrix(j, p) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
            end
        end           
    end
    if(p == 1) 
         optionPriceMatrix(1, 1) = expNeg*(q*optionPriceMatrix(j,p+1) + (1-q)*optionPriceMatrix(j+1,p+1));
         result = optionPriceMatrix;
    else
        result = calcOptionPrice(optionPriceMatrix, p-1, indexPriceMatrix, K, q, expNeg, barrierPrice);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Help functions
function indexPriceMatrix = calcIndexPrice(p, u, d, S)
    indexPriceMatrix = zeros(p+1);
    for i = 0:p %rows
        for j = 0:i %columns
            indexPriceMatrix(j+1, i+1) = S*u^(i-j)*d^j;
        end
    end 
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

function q = calculateQ(r, d, u, t)
    q = (exp(r*t)-d)/(u-d);
end


