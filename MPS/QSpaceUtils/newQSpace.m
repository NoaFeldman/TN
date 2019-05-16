function new = newQSpace(Q, data, itags)
    new = QSpace();
    new.info.qtype = '';
    new.info.otype = '';
    new.Q = Q;
    new.data = data;
    new.info.itags = itags;
end
    