function kron = kronackerProductQ(op1, op2)
    id = getIdentity(op1, op2);
    kron = contract(id, '12*', contract(op2, 2, contract(op1, 2, id, 1), 2), '12');
end