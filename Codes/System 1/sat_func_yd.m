function yd_sat = sat_func_yd(yd,yd_l,yd_u)

    yd_sat = (yd_l).*(yd<yd_l) + (yd).*(yd>=yd_l & yd<=yd_u) + (yd_u).*(yd>yd_u);

end