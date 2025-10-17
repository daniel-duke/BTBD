%%% angle between two vectors
function phi = calcDihedral(v1,v2,v3)
    n1 = cross(v1,v2);
    n2 = cross(v2,v3);
    x = dot(n1,n2);
    y = dot(v2,cross(n1,n2));
    phi = atan2d(y,x);
end