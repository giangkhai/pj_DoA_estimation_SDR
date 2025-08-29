Set 
    m /1*14/, 
    n /1*13/, 
    l /1*2/;

Alias (m, mm); 
Alias (n, nn);

Parameter ord_m(m), ord_n(n), sum_idx(m,n,mm);
Scalar p /5/;

ord_m(m) = ord(m);
ord_n(n) = ord(n);

Loop((m,n,mm),
    sum_idx(m,n,mm)$(ord(mm) = ord(m) + ord(n)) = 1;
);

Binary Variable a(m,l);         
Binary Variable y(n);           
Binary Variable z(m,n,l);       
Positive Variable w(n);        
Variable obj;

Equation obj_def, lin1(m,n,l), lin2(m,n,l), lin3(m,n,l), y_chain(n), y_upper(n), y_lower(n), w_eq(n), count(l), obj_lb, obj_ub;

obj_def.. obj =e= sum(n, y(n));

lin1(m,n,l).. z(m,n,l) =l= a(m,l);
lin2(m,n,l)$(ord_m(m) + ord_n(n) <= card(m)).. z(m,n,l) =l= sum(mm$(sum_idx(m,n,mm)), a(mm,l));
lin3(m,n,l)$(ord_m(m) + ord_n(n) <= card(m)).. z(m,n,l) =g= a(m,l) + sum(mm$(sum_idx(m,n,mm)), a(mm,l)) - 1;

y_chain(n)$(ord_n(n) > 1).. y(n) =l= sum(nn$(ord_n(nn) = ord_n(n) - 1), y(nn));
w_eq(n).. w(n) =e= sum((m,l)$(ord(m) <= card(m) - ord(n)), z(m,n,l));

count(l).. sum(m, a(m,l)) =e= p;
y_upper(n).. w(n) =g= y(n);
y_lower(n).. w(n) =l= card(m) * y(n);

obj_lb.. obj =g= card(m) - 1;
obj_ub.. obj =l= card(l) * card(m) * (card(m)-1);

Model ula_dof_max /all/;

option optcr = 0.01;
option reslim = 15000;

Solve ula_dof_max using mip maximizing obj;

Display obj.l, y.l, w.l, a.l;
