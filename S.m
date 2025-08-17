function s = S(z)

if z > 0
    s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
elseif z < 0
    s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
else
    s = 1/6;
end

