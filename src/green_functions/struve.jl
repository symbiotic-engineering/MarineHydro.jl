function StruveH0(xx)
    if xx <= 3.0
        yy = (xx / 3.0)^2
        
        P0 = +1.909859164
        P1 = -1.909855001
        P2 = +0.687514637
        P3 = -0.126164557
        P4 = +0.013828813
        P5 = -0.000876918
        
        StruveH0 = P0 + (P1 + (P2 + (P3 + (P4 + P5 * yy) * yy) * yy) * yy) * yy
        StruveH0 *= xx / 3.0
    else
        yy = (3.0 / xx)^2
        
        a0 = 0.99999906
        a1 = 4.77228920
        a2 = 3.85542044
        a3 = 0.32303607

        b1 = 4.88331068
        b2 = 4.28957333
        b3 = 0.52120508

        c1 = 2.0 * (a0 + (a1 + (a2 + a3 * yy) * yy) * yy)
        c2 = pi * xx * (1.0 + (b1 + (b2 + b3 * yy) * yy) * yy)
        
        StruveH0 = c1 / c2 + bessely0(xx)
    end
    
    return StruveH0
end



function StruveH1(xx)
    if xx <= 3.0
        yy = (xx / 3.0)^2
        
        P1 = +1.909859286
        P2 = -1.145914713
        P3 = +0.294656958
        P4 = -0.042070508
        P5 = +0.003785727
        P6 = -0.000207183
      
        StruveH1 = (P1 + (P2 + (P3 + (P4 + (P5 + P6 * yy) * yy) * yy) * yy) * yy) * yy
    else
        yy = (3.0 / xx)^2

        a0 = 1.00000004
        a1 = 3.92205313
        a2 = 2.64893033
        a3 = 0.27450895

        b1 = 3.81095112
        b2 = 2.26216956
        b3 = 0.10885141

        c1 = 2.0 * (a0 + (a1 + (a2 + a3 * yy) * yy) * yy)
        c2 = pi * (1.0 + (b1 + (b2 + b3 * yy) * yy) * yy)

        StruveH1 = c1 / c2 + bessely(1, xx)
    end
    
    return StruveH1
end
