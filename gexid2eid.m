function eid = gexid2eid(id)
k = strfind(id, '_');
eid = str2double(id(1:k-1));