function [Eex, Hex, GSex] = exact(N, h, JPM, JZ)
   statesNum = 2^N;
   H = zeros(statesNum, statesNum);
   for i = 0 : statesNum - 1
       for j = 0 : i
           if (i == j)
               state = i;
               for bit = 1: N
                   rightSpin = bitand(state, 1) - 0.5;
                   leftSpin = bitand(state, 2) / 2 - 0.5;
                   H(i+1, i+1) = H(i+1, i+1) - h * rightSpin;
                   if (bit < N)
                       H(i+1, i+1) = H(i+1, i+1) + ...
                           JZ * rightSpin  * leftSpin;
                   state = floor(state / 2);
               end
           end
           else
               xor = bitxor(i, j);
               for bit = 1:N
                   if (xor == 3)
                       xor = bitxor(i, j);
                       if (bitand(xor, i) ~= 0 & bitand(xor, j) ~= 0)
                           H(i+1, j+1) = H(i+1, j+1) + JPM / 2;
                           H(j+1, i+1) = H(i+1, j+1);
                       end
                       break;
                   elseif (bitand(xor, 1) == 1)
                       break;
                   end
                   xor = floor(xor / 2);   
               end
           end
       end
   end
%   E = eigs(H,4,'smallestreal');
    [V, E] = eig(H);
    Eex = E(1, 1);
    Hex = H;
    GSex = V(1:2^N, 1);
end