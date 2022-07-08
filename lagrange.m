% Lagrange
% Washington Yandun

function [P,R,S] = lagrange(X,Y,XX)
X = [1950 1962 1974 1982 1990 2001 2010]
Y = [386520 587835 988306 1382125 1756228 2388817 2576287]

xx = 1950 : 1 : 2010

% check dimension
if size(X,1) > 1
    X = X'
end 

if size(Y,1) > 1
    Y = Y'
end

if size(X,1) > 1 || size(Y,1) > 1 || size(X,2) ~= size(Y,2)
  error('vectores disparejos')
end

% numero de puntos
n = length(X)

pvalores = zeros(n,n)

for i = 1:n
  pp = poly(X( (1:n) ~= i))
  pvalores(i,:) = pp ./ polyval(pp, X(i))
end

P = Y*pvalores

if nargin == 3
  YY = polyval(P,XX)
  P = YY
end

if nargout > 1 
  R = roots( ((n-1):-1:1) .* P(1:(n-1)) )
  if nargout > 2
    S = polyval(P,R)
  end
end
