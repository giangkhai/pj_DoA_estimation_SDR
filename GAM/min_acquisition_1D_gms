Sets
    m /1*20/,
    n /1*19/,
    l /1*19/;  
Alias (m, mm);
Alias (n, nn);

* Parameters
Parameter ord_m(m), ord_n(n), sum_idx(m,n,mm);
Scalar p /3/;
ord_m(m) = ord(m);
ord_n(n) = ord(n);

* Precompute index mapping: m + n = mm
Loop((m,n,mm),
    sum_idx(m,n,mm)$(ord(mm) = ord(m) + ord(n)) = 1;
);

Binary Variable a(m,l);
Binary Variable u(l);

* z(m,n,l) := a(m,l) * a(m+n,l) * u(l)
* y(m,l) := a(m,l) * u(l)
Binary Variable z(m,n,l);
Binary Variable y(m,l);
Variable obj;

* Objective: minimize sum_l u_l
Equation obj_def;
obj_def.. obj =e= sum(l, u(l));

* Constraints
Equation
    lin1(m,n,l),
    lin2(m,n,l),
    lin3(m,n,l),
    lin5(m,l),
    lin6(m,l),
    lin7(m,l),
    constraint1(n),
    constraint2(l)
    min_L;

lin5(m,l)..
    y(m,l) =l= a(m,l);

lin6(m,l)..
    y(m,l) =l= u(l);

lin7(m,l)..
    y(m,l) =g= a(m,l) + u(l) - 1;

lin1(m,n,l)..
    z(m,n,l) =l= y(m,l);

lin2(m,n,l)$(ord_m(m) + ord_n(n) <= card(m))..
    z(m,n,l) =l= sum(mm$(sum_idx(m,n,mm)), a(mm,l));

lin3(m,n,l)$(ord_m(m) + ord_n(n) <= card(m))..
    z(m,n,l) =g= u(l) + y(m,l) + sum(mm$(sum_idx(m,n,mm)), a(mm,l)) - 2;

constraint2(l).. sum(m, y(m,l)) =e= p*u(l);

constraint1(n).. sum((m,l)$(ord(m) <= card(m) - ord(n)), z(m,n,l)) =g= 1;

min_L.. sum(l,u(l)) =g= (card(m)-1)*2/(p*(p-1)) ;


Model minimize_L /all/;

option optcr = 0;
option reslim = 300;

Solve minimize_L using mip minimizing obj;

Display obj.l, a.l, u.l, y.l, z.l;
