function s2 = stheo2(L1, L2, NA, K, w1, w2)
   c = 1;
    s = (1/12).*K.^(-3/2).*L2.^((-1/2).*K).*log(L2).^(-3/2).*(sqrt(-1).* ...
         2.^(1/2).*exp(1).^((sqrt(-1)*(-1)).*NA.*pi).*(dawson(2.^(-1/2).* ...
          K.^(-1/2).*log(L2).^(-1/2).*(NA.*pi+(sqrt(-1)*(-1)).*K.*log(L2)))+ ...
          (-1).*exp(1).^((sqrt(-1)*2).*NA.*pi).*dawson(2.^(-1/2).*K.^(-1/2) ...
          .*log(L2).^(-1/2).*(NA.*pi+sqrt(-1).*K.*log(L2)))).*(3.*NA.^2.* ...
          pi.^2+K.*((-3)+c.*log(L1)).*log(L2))+6.*K.^(1/2).*log(L2).^( ...
          1/2).*(K.*cos(NA.*pi).*log(L2)+(-1).*NA.*pi.*sin(NA.*pi)));
   s2 = w1.*s + w2.*s2ndOrder(L1, L2, NA, K); 
end

function s2 = s2ndOrder(L1, L2, NA, K)
    c = 1;
    s2 = (1/12).*2.^(-1/2).*exp(1).^((sqrt(-1)*(-1)).*NA.*pi).*K.^(-3/2).* ...
  L2.^((-1/2).*K).*log(L2).^(-3/2).*((sqrt(-1)*(-2)).*dawson(2.^( ...
  -1/2).*K.^(-1/2).*log(L2).^(-1/2).*(NA.*pi+(sqrt(-1)*(-1)).*K.* ...
  log(L2))).*(3.*NA.^2.*pi.^2+K.*((-3)+c.*log(L1)).*log(L2))+(sqrt( ...
  -1)*2).*exp(1).^((sqrt(-1)*2).*NA.*pi).*dawson(2.^(-1/2).*K.^( ...
  -1/2).*log(L2).^(-1/2).*(NA.*pi+sqrt(-1).*K.*log(L2))).*(3.* ...
NA.^2.*pi.^2+K.*((-3)+c.*log(L1)).*log(L2))+(-6).*2.^(1/2).*exp(1) ...
  .^(sqrt(-1).*NA.*pi).*K.^(1/2).*log(L2).^(1/2).*(K.*cos(NA.*pi).* ...
  log(L2)+(-1).*NA.*pi.*sin(NA.*pi))+exp(1).^((sqrt(-1)*(-2)).*NA.* ...
  pi).*L2.^((-4).*K).*((sqrt(-1)*2).*dawson(2.^(-1/2).*K.^(-1/2).* ...
  log(L2).^(-1/2).*(NA.*pi+(sqrt(-1)*(-3)).*K.*log(L2))).*(3.* ...
NA.^2.*pi.^2+K.*((-3)+c.*log(L1)).*log(L2))+(sqrt(-1)*(-2)).*exp( ...
  1).^((sqrt(-1)*6).*NA.*pi).*dawson(2.^(-1/2).*K.^(-1/2).*log(L2) ...
  .^(-1/2).*(NA.*pi+(sqrt(-1)*3).*K.*log(L2))).*(3.*NA.^2.*pi.^2+K.* ...
  ((-3)+c.*log(L1)).*log(L2))+6.*2.^(1/2).*exp(1).^((sqrt(-1)*3).* ...
  NA.*pi).*K.^(1/2).*log(L2).^(1/2).*(3.*K.*cos(3.*NA.*pi).*log(L2)+ ...
  (-1).*NA.*pi.*sin(3.*NA.*pi))));
end