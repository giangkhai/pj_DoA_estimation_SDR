* X = 1..M, x-axis
* Y = 1..N, y-axis
* card m = card X, card n = card Y

Set 
   X  /1*4/,
   Y  /1*4/,
   m  /1*3/,
   n  /1*3/,
   l  /1*8/;
Alias (X, xm);
Alias (Y, yn);
Alias (X, subxm);
* Parameters
Parameter ord_x(X), ord_y(Y), ord_m(m), ord_n(n), sum_idx(X,m,xm), sum_idy(Y,n,yn), sub_idx(X,m,subxm);

* Number of physical antennas
Scalar p /3/;
ord_x(X) = ord(X);
ord_y(Y) = ord(Y);
ord_m(m) = ord(m);
ord_n(n) = ord(n);

* Precompute index mapping: X + m = xm, Y + n = yn 
Loop((X,m,xm),
    sum_idx(X,m,xm)$(ord(xm) = ord(X) + ord(m) - 1) = 1;
);

Loop((Y,n,yn),
    sum_idy(Y,n,yn)$(ord(yn) = ord(Y) + ord(n) - 1) = 1;
);

Loop((X,m,subxm),
    sub_idx(X,m,subxm)$(ord(subxm) = ord(X) - ord(m) + 1) = 1;
);

Binary Variable a(X,Y,l);
Binary Variable u(l);


* y(X, Y, l)   := a(X, Y, l) * u(l)
* z(X ,m ,Y ,n ,l) := a(X, Y, l) * a(X + m, Y + n, l) * u(l);


Binary Variable yvar(X,Y,l);
Binary Variable zvar(X,Y,m,n,l);
Binary Variable zsub(X,Y,m,n,l);
Variable obj;

* Minimize number of used acquisition card(u(l)) 
Equation obj_def;
obj_def.. obj =e= sum(l, u(l));


Equation
    lin1_y(X,Y,l),
    lin2_y(X,Y,l),
    lin3_y(X,Y,l),
    lin1_z(X,Y,m,n,l),
    lin2_z(X,Y,m,n,l),
    lin3_z(X,Y,m,n,l),
    lin1_zsub(X,Y,m,n,l),
    lin2_zsub(X,Y,m,n,l),
    lin3_zsub(X,Y,m,n,l),
    constraint1(m,n),
    constraint2(l),
    constraintsub(m,n),
    min_L;



* Constraints
lin1_y(X,Y,l)..
    yvar(X,Y,l) =l= a(X,Y,l);

lin2_y(X,Y,l)..
    yvar(X,Y,l) =l= u(l);

lin3_y(X,Y,l)..
    yvar(X,Y,l) =g= a(X,Y,l) + u(l) - 1;

lin1_z(X,Y,m,n,l)..
    zvar(X,Y,m,n,l) =l= yvar(X,Y,l);

lin2_z(X,Y,m,n,l).. 
    zvar(X,Y,m,n,l) =l= sum((xm,yn)$(sum_idx(X,m,xm) and sum_idy(Y,n,yn)), a(xm,yn,l));

lin3_z(X,Y,m,n,l)..
    zvar(X,Y,m,n,l) =g= yvar(X,Y,l) + sum((xm,yn)$(sum_idx(X,m,xm) and sum_idy(Y,n,yn)), a(xm,yn,l)) - 1;


lin1_zsub(X,Y,m,n,l)..
    zsub(X,Y,m,n,l) =l= yvar(X,Y,l);

lin2_zsub(X,Y,m,n,l).. 
    zsub(X,Y,m,n,l) =l= sum((subxm,yn)$(sub_idx(X,m,subxm) and sum_idy(Y,n,yn)), a(subxm,yn,l));

lin3_zsub(X,Y,m,n,l)..
    zsub(X,Y,m,n,l) =g= yvar(X,Y,l) + sum((subxm,yn)$(sub_idx(X,m,subxm) and sum_idy(Y,n,yn)), a(subxm,yn,l)) - 1;


constraint2(l).. sum((X,Y), yvar(X,Y,l)) =e= p*u(l);

constraint1(m,n).. sum((X,Y,l)$(ord(X)+ord(m)-1 <= card(X) and ord(Y)+ord(n)-1 <= card(Y)),
         zvar(X,Y,m,n,l)) =g= 1;

constraintsub(m,n).. sum((X,Y,l)$(ord(X)-ord(m)+1 <= card(X) and ord(Y)+ord(n)-1 <= card(Y)),
         zsub(X,Y,m,n,l)) =g= 1;

min_L.. sum(l,u(l)) =g= (card(X)*card(Y)-1)*2/(p*(p-1)) ;


Model minimize_L /all/;

option optcr = 0;
option reslim = 3000;

Solve minimize_L using mip minimizing obj;

Display obj.l, a.l, u.l, yvar.l, zvar.l, zsub.l;
