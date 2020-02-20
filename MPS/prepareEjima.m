startWork
Ns = [12 128:32:224];
for i = 1:length(Ns)
    N = Ns(i);
    [gs, ~, ~, ~, EGS, Es] = getGroundState(N, 0, 1, -0.5, 0.25, 0.5*0.25, 0, 'periodic'); 
    [gs1, ~, ~, ~, EGS1, Es1] = getGroundState(N, 0, 1, -0.5, 0.25, 0.5*0.25, -1, 'periodic');
    save(strcat('/home/noa/TN/MPS/ejimaGSJZ-05J2PM025/N', int2str(N)), 'gs', 'EGS', 'Es', 'gs1', 'EGS1', 'Es1');
end

% for i = 1:length(Ns)
% load(strcat('ejimaGSJZ-05J2PM002/N', int2str(Ns(i))));
% eq2(i) = EGS/Ns(i);
% eq3(i) = EGS1 - EGS;
% end
