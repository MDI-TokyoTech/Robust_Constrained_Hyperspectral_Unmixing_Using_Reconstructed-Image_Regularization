function[result] = ProjL2Ball(x,b,ep)
   % radius = norm(x(:) - b(:),2);
   % if radius > ep
   %     x = b + (ep/radius)*(x - b);
    l2distance = sqrt(sum((x - b).^2, "all"));
    if l2distance <= ep
        result = x;
    else
        result = b + ep*(x - b)/l2distance;
    end
    end