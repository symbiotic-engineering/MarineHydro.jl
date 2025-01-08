struct GFWu <: GreensFunction end

function greens(::GFWu, element_1, element_2, wavenumber)
    x  = element_1.center
    xi = element_2.center
    hh = wavenumber * hypot(x[1] - xi[1], x[2] - xi[2])
    if hh==0.0
        hh = hh + 1e-18
    end
    vv = wavenumber * (x[3] + xi[3])
    dd = sqrt(hh * hh + vv * vv)
    alpha = -vv / dd
    beta = hh / dd
    rho = dd / (1.0 + dd)
    return wavenumber * (-GF_Func_L0(vv, dd, alpha, beta, rho) - GF_Func_W(hh, vv))
end

function gradient_greens(::GFWu, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    x = element_1.center
    xi = element_2.center
    hh = wavenumber * hypot(x[1] - xi[1], x[2] - xi[2])
    if hh==0.0
        hh = hh + 1e-18
    end
    vv = wavenumber * (x[3] + xi[3])
    dd = sqrt(hh * hh + vv * vv)
    alpha = -vv / dd
    beta = hh / dd
    rho = dd / (1.0 + dd)
    dGF_dhh = -GF_Func_Ls(hh, vv, dd, alpha, beta, rho) - GF_Func_Wh(hh, vv)
    dGF_dvv = -GF_Func_L0(vv, dd, alpha, beta, rho) - GF_Func_W(hh, vv) + 2/dd
    if with_respect_to_first_variable
        if abs(hh) > 1e-6
            dhh_dx = wavenumber^2 / hh * [x[1] - xi[1], x[2] - xi[2]]
        else
            dhh_dx = [zero(hh), zero(hh)]
        end
        dvv_dx = wavenumber
        return wavenumber * [dhh_dx * dGF_dhh; dvv_dx * dGF_dvv]
    else
        if abs(hh) > 1e-6
            dhh_dxi = wavenumber^2 / hh * [xi[1] - x[1], xi[2] - x[2]]
        else
            dhh_dxi = [zero(hh), zero(hh)]
        end
        dvv_dxi = wavenumber
        return wavenumber * [dhh_dxi * dGF_dhh; dvv_dxi * dGF_dvv]
    end
end

function integral(g::GFWu, element_1, element_2, wavenumber)
    # One-point approximation of the integral
    return greens(g::GFWu, element_1, element_2, wavenumber) * element_2.area
end

function integral_gradient(g::GFWu, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    # One-point approximation of the integral
    return gradient_greens(g::GFWu, element_1, element_2, wavenumber; with_respect_to_first_variable) * element_2.area
end

######################################################

function GF_dFuncC(tt)
     C0 = 14.19
     C1 = -148.24
     C2 = 847.8
     C3 = -2318.58
     C4 = 3168.35
    C5 = -1590.27

    GF_dFuncC = (C4 + C5 * tt) * tt
    GF_dFuncC = (C3 + GF_dFuncC) * tt
    GF_dFuncC = (C2 + GF_dFuncC) * tt
    GF_dFuncC = (C1 + GF_dFuncC) * tt
    GF_dFuncC = C0 + GF_dFuncC

    return GF_dFuncC
end


function GF_dFuncB(tt)
    B0 = 1.11
    B1 = 2.894
    B2 = -76.765
    B3 = 1565.35
    B4 = -11336.19
    B5 = 44270.15
    B6 = -97014.11
    B7 = 118879.26
    B8 = -76209.82
    B9 = 19923.28

    GF_dFuncB = (B8 + B9 * tt) * tt
    GF_dFuncB = (B7 + GF_dFuncB) * tt
    GF_dFuncB = (B6 + GF_dFuncB) * tt
    GF_dFuncB = (B5 + GF_dFuncB) * tt
    GF_dFuncB = (B4 + GF_dFuncB) * tt
    GF_dFuncB = (B3 + GF_dFuncB) * tt
    GF_dFuncB = (B2 + GF_dFuncB) * tt
    GF_dFuncB = (B1 + GF_dFuncB) * tt
    GF_dFuncB = B0 + GF_dFuncB

    return GF_dFuncB
end

function GF_dFuncA(tt)
    A0 = 2.948
    A1 = -24.53
    A2 = 249.69
    A3 = -754.85
    A4 = -1187.71
    A5 = 16370.75
    A6 = -48811.41
    A7 = 68220.87
    A8 = -46688.0
    A9 = 12622.25

    GF_dFuncA = (A8 + A9 * tt) * tt
    GF_dFuncA = (A7 + GF_dFuncA) * tt
    GF_dFuncA = (A6 + GF_dFuncA) * tt
    GF_dFuncA = (A5 + GF_dFuncA) * tt
    GF_dFuncA = (A4 + GF_dFuncA) * tt
    GF_dFuncA = (A3 + GF_dFuncA) * tt
    GF_dFuncA = (A2 + GF_dFuncA) * tt
    GF_dFuncA = (A1 + GF_dFuncA) * tt
    GF_dFuncA = A0 + GF_dFuncA

    return GF_dFuncA
end

#import Pkg; Pkg.add("SpecialFunctions")
using SpecialFunctions

function GF_Func_Wh(hh, vv)
    H1 = StruveH1(hh)
    J1 = besselj(1, hh)

    GF_Func_Wh = 2.0 * pi * (2.0 / pi - H1 + im * J1) * exp(vv)

    return GF_Func_Wh
end

function GF_Func_Lsp(hh, vv, dd, alpha, beta, rho)
    A = GF_dFuncA(rho)
    B = GF_dFuncB(rho)
    C = GF_dFuncC(rho)

    RR = beta * A
    RR -= (1.0 - alpha) * B
    RR += beta * (1.0 - beta) * rho * (1.0 - 2.0 * rho) * C

    GF_Func_Lsp = rho * (1.0 - rho)^3 * RR

    return GF_Func_Lsp
end

function GF_Func_Ls(hh, vv, dd, alpha, beta, rho)
    PS = (beta + hh) / (dd - vv)
    PS = PS - 2.0 * beta + 2.0 * exp(vv) * dd - hh

    QS = exp(-dd) * (1.0 - beta)
    QS = QS * (1.0 + dd / (1.0 + dd^3))

    Lsp = GF_Func_Lsp(hh, vv, dd, alpha, beta, rho)

    GF_Func_Ls = 2.0 * PS / (1.0 + dd^3) - 4.0 * QS + 2.0 * Lsp

    return GF_Func_Ls
end

function GF_FuncD(tt)
    D0 = 0.632
    D1 = -40.97
    D2 = 667.16
    D3 = -6072.07
    D4 = 31127.39
    D5 = -96293.05
    D6 = 181856.75
    D7 = -205690.43
    D8 = 128170.2
    D9 = -33744.6

    GF_FuncD = (D8 + D9 * tt) * tt
    GF_FuncD = (D7 + GF_FuncD) * tt
    GF_FuncD = (D6 + GF_FuncD) * tt
    GF_FuncD = (D5 + GF_FuncD) * tt
    GF_FuncD = (D4 + GF_FuncD) * tt
    GF_FuncD = (D3 + GF_FuncD) * tt
    GF_FuncD = (D2 + GF_FuncD) * tt
    GF_FuncD = (D1 + GF_FuncD) * tt
    GF_FuncD = D0 + GF_FuncD

    return GF_FuncD
end

function GF_FuncC(tt)
    C0 = 1.268
    C1 = -9.747
    C2 = 209.653
    C3 = -1397.89
    C4 = 5155.67
    C5 = -9844.35
    C6 = 9136.4
    C7 = -3272.62

    GF_FuncC = (C6 + C7 * tt) * tt
    GF_FuncC = (C5 + GF_FuncC) * tt
    GF_FuncC = (C4 + GF_FuncC) * tt
    GF_FuncC = (C3 + GF_FuncC) * tt
    GF_FuncC = (C2 + GF_FuncC) * tt
    GF_FuncC = (C1 + GF_FuncC) * tt
    GF_FuncC = C0 + GF_FuncC

    return GF_FuncC
end

function GF_FuncB(tt)
    B0 = 0.938
    B1 = 5.737
    B2 = -67.92
    B3 = 796.534
    B4 = -4780.77
    B5 = 17137.74
    B6 = -36618.81
    B7 = 44894.06
    B8 = -29030.24
    B9 = 7671.22

    GF_FuncB = (B8 + B9 * tt) * tt
    GF_FuncB = (B7 + GF_FuncB) * tt
    GF_FuncB = (B6 + GF_FuncB) * tt
    GF_FuncB = (B5 + GF_FuncB) * tt
    GF_FuncB = (B4 + GF_FuncB) * tt
    GF_FuncB = (B3 + GF_FuncB) * tt
    GF_FuncB = (B2 + GF_FuncB) * tt
    GF_FuncB = (B1 + GF_FuncB) * tt
    GF_FuncB = B0 + GF_FuncB

    return GF_FuncB
end

function GF_FuncA(tt)
    A0 = 1.21
    A1 = -13.328
    A2 = 215.896
    A3 = -1763.96
    A4 = 8418.94
    A5 = -24314.21
    A6 = 42002.57
    A7 = -41592.9
    A8 = 21859.0
    A9 = -4838.6

    GF_FuncA = (A8 + A9 * tt) * tt
    GF_FuncA = (A7 + GF_FuncA) * tt
    GF_FuncA = (A6 + GF_FuncA) * tt
    GF_FuncA = (A5 + GF_FuncA) * tt
    GF_FuncA = (A4 + GF_FuncA) * tt
    GF_FuncA = (A3 + GF_FuncA) * tt
    GF_FuncA = (A2 + GF_FuncA) * tt
    GF_FuncA = (A1 + GF_FuncA) * tt
    GF_FuncA = A0 + GF_FuncA

    return GF_FuncA
end



function GF_Func_W(hh, vv)
    H0 = StruveH0(hh)
    J0 = besselj(0, hh)

    GF_Func_W = 2.0 * pi * (H0 - im * J0) * exp(vv)

    return GF_Func_W
end

function GF_Func_Lp( alpha, beta, rho)
    A = GF_FuncA(rho)
    B = GF_FuncB(rho)
    C = GF_FuncC(rho)
    D = GF_FuncD(rho)

    RR = (1.0 - beta) * A
    RR -= beta * B
    RR -= alpha * C / (1.0 + 6.0 * alpha * rho * (1.0 - rho))
    RR += beta * (1.0 - beta) * D

    GF_Func_Lp = rho * (1.0 - rho)^3 * RR

    return  GF_Func_Lp
end

function GF_Func_L0(vv, dd, alpha, beta, rho)
    gama = 0.5772156649
    PP = log(0.5 * (dd - vv)) + gama - 2.0 * dd^2
    PP = exp(vv) * PP
    PP = PP + dd^2 - vv

    Lp = GF_Func_Lp( alpha, beta, rho)

    GF_Func_L0 = 2.0 * PP / (1.0 + dd^3) + 2.0 * Lp

    return   GF_Func_L0
end
