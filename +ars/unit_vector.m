%%% turn vector into unit vector
function v = unit_vector(v)
    v = v./norm(v);
end