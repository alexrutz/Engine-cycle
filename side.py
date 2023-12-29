from math import sin, cos, log, sqrt, pi, e
from cp import cpd, cps

Tref = T0 = 288.15
ptref = p0 = 101325
R = 287.05
y = 45 * 10 ** (-6)
lambda_hg = 60 * 10 ** (-3)
Pr = 0.7

#----------materialvalue

def cv(cp, R):

    return cp - R

def k(cp):

    return cp /(cp - R)

#----------ambient

def Th(h):

    return 288.15 - 0.0065 * h

def ph(Th):

    return 101325 * (Th/288.15) ** 5.255

def rh(ph, R, Th):

    return ph /(R * Th)

#----------temps

def T(Tt, k, Ma):

    return Tt / (1 + ((k - 1) / 2) * Ma ** 2)

def Tt(T, k, Ma):

    return T * (1 + ((k - 1) / 2) * Ma ** 2)

def Tiscomp(T1, Pic, R, cpd):

    return T1 * Pic ** (R / cpd)

def Tiscompdps(Tt1, Pic, R, cpd, nvps, ix):

    return (Tt1 * Pic ** ((R / cpd) * (1 / nvps) ** ix)) * (1 / nvps) ** (1 - ix)

def dTiscompdps(Tt1, Pic, R, cpds, nvps, ix):

    return (Tt1 * (Pic ** ((R / cpds) * (1 / nvps) ** ix) - 1)) * (1 / nvps) ** (1 - ix)

def Tt4max_func(OTDF, Tt3, Tt4):

    return OTDF * (Tt4 - Tt3) + Tt4

def Ttm_func(T0, Ttc, Ttms, cpdm, y_c, mpsus, betaus):

    return ((mpsus * cpd(T0, Ttms, betaus) * (Ttms - T0) + y_c * cpd(T0, Ttc, 0) * (Ttc - T0)) / ((mpsus + y_c) * cpdm)) + T0

def Tt_exp_func(at, T1, cpdt1t2):

    return (at / cpdt1t2) + T1

def T_contour(Tt, c, cpd):

    return Tt - ((c ** 2) / (2 * cpd))

#----------burner/turbine

def Ppe_func(at, Tt1, Tt2, Tt2s, neps, ix, beta):

    return ((at /(cpd(Tt1, Tt2, beta) * Tt1 * neps ** (1 - ix))) + 1) ** (-(cpd(Tt1, Tt2s, beta) / R) * neps ** (-ix))

def betaBK_func(T0, Tt3, Tt4, cpd0t4, beta4, y_kl, Hu, nBK):
    
    return ((1 - y_kl) * (cpd(T0, Tt3, 0) * (Tt3 - T0) - cpd0t4 * (Tt4 - T0))) / ((cpd0t4 * (Tt4 - T0) - (Hu * nBK)))

def beta4_func(betaBK, y_kl):

    return betaBK / (1 - y_kl)

#----------eff

def nvs(cp12, cp12s, Tt1, Pic, R, av):

    return (cp12 * Tt1 * (Pic ** (R / cp12s) - 1)) / av

def nvp(cp12, cp12s, Tt1, Pic, R, av):

    return (R * log(Pic)) / (cp12s * log((av /(cp12 * Tt1)) + 1))

def nes(cp12, cp12s, Tt1, Ppe, R, at):

    return ((cp12 * Tt1 * ((1 / Ppe) ** (R / cp12s) - 1))/ at) ** (-1)

def nep(cp12, cp12s, Tt1, Ppe, R, at):

    return ((R * log(1 / Ppe)) / (cp12s * log((at / (cp12 * Tt1)) + 1))) ** (-1)

#----------perfmetrics

def av(cp12, cp12s, Tt1, Pic, nvps, ix):

    return cp12 * Tt1 * (Pic ** ((R/cp12s)*((1/nvps) ** ix))) * (1/nvps) ** (1 - ix)

def at(avus, mps, nm):

    return -avus /(mps * nm)

def De(cp12s, T1, T2s):

    return cp12s * (T2s - T1)

def Db(av, De):

    return av - De

def P(mp, av):

    return mp * av

#----------ncs

def T4crit(Tt4, cpd):

    return ((1 + (k(cpd) - 1) / 2 ) ** (- 1)) * Tt4

def mpred(Tt4, beta4):

    return sqrt(k(cps(Tt4, beta4)) / R) * (2 / (k(cps(Tt4, beta4)) + 1)) ** ((k(cps(Tt4, beta4)) + 1) / (2 * ((k(cps(Tt4, beta4)) -1))))

def A4crit_func(mpred4, pt4, Tt4, mp4):

    return (mp4 * sqrt(Tt4)) / (mpred4 * pt4)

#----------entro

def Ds(EF12, p1, p2):

    return (EF12 - log(p2 / p1)) * R

def Dh(T1, T2, beta):

    return cpd(T1, T2, beta) * (T2 - T1)

#----------stat

def p(pt, T, Tt, beta):

    return pt * (T / Tt) ** (cpd(Tt, T, beta) / R)

def rho(p, T):

    return p / (R * T)

def Vp(mp, rho):

    return mp / rho

#----------vel

def a(R, T, beta):

    return sqrt(k(cps(T, beta)) * R * T)

def c(Ma, a):

    return Ma * a

def cax(c, alpha):

    return c * cos(alpha)

def Aax(Vp, cax):

    return Vp / cax

#----------comp

def Pbp_func(Tt, T, Ma):

    return (1 + ((k(cpd(T, Tt, 0)) - 1) / 2) * Ma ** 2) ** (cpd(T, Tt, 0) / R)

def PC_func(Cr, PI, Pbp):

    return Cr /(PI * Pbp)

def PLPC_func(PC, P1):

    return sqrt(PC * P1)

def PHPC_func(PC, P1):

    return sqrt(PC / P1)

def PHPC1_func(PHPC, P2):

    return sqrt(PHPC * P2)

def PHPC2_func(PHPC, P2):

    return sqrt(PHPC / P2)

#----------nozzle

def M(Tt, T, k):

    return sqrt(((Tt / T) - 1) * (2 / (k - 1)))

def fs(c0, c9, betaB):

    return (1 + betaB) * c9 - c0

def fsEng_func(fsI, fsII, bypass):

    return (fsI + (bypass * fsII)) / (1 + bypass)

def bsEng_func(betaB, bypass, fsEng):

    return betaB / ((1 + bypass) * fsEng)

def mp1_func(F, fsEng):

    return (F * 1000) / fsEng

def mp2_func(mp1, bypass):

    return mp1 / (1 + bypass)

def mp12_func(mp2, bypass):

    return bypass * mp2

#----------efficiencies

def aFlow_func(c0, c9, c19, betaB, bypass):

    return ((1 + betaB) * (c9 ** 2) - (c0 ** 2) + bypass * ((c19 ** 2) - (c0 ** 2))) / 2

def nth_func(aFlow, betaB, Hu):

    return aFlow/ (betaB * Hu)

def nP_func(c0, c9, c19, betaB, bypass, aFlow):

    return (((1 + betaB) * c9 - c0 + bypass * (c19 - c0)) * c0) / aFlow

def nEng_func(nth, nP):

    return (nth * nP)

def nCore_func(c0, Tt9, T9, betaB, Hu, bypass, nm, aF):

    return ((bypass * aF / nm) + (1 + betaB) * cpd(Tt9, T9, betaB) * (Tt9 - T9) - (1 / 2) * (1 + bypass) * (c0 ** 2)) / (betaB * Hu)

#----------geometry

def Da(Dm, A):

    return Dm + (A /(pi * Dm))

def Di(Dm, A):

    return Dm - (A /(pi * Dm))

def D_euler(Di, Da):

    return sqrt(((Di ** 2) + (Da ** 2)) / 2)

def Nmax(umax, Da):

    return (umax * 60) /(pi * Da)

def A(Di, Da):

    return (pi / 4) * (Da - Di) ** 2

#----------unassigned

def Ttisobaric(T1, s1, s, cps):

    return T1 * e ** ((s - s1) / cps)

#----------compressor/turbine_detail

def stagepressure(Ttz1, Ttz2, cpd, nvp):

    return (Ttz2 / Ttz1) ** ((cpd * nvp) / R)

def w_func(Da1, Da2, Di1, Di2, ar):

    return ((Da2 - Di2 + Da1 - Di1) / 2) * (1 / ar)

def D(D1, D2):

    return (D1 + D2) / 2

def D21_func(D2, D3, wR, wS):

    return ((D2 / wR) + (D3 / wS)) / ((1 / wR) + (1 / wS))

def DmE_func(Da, Di):

    return sqrt(((Da ** 2) + (Di ** 2)) / 2)

def DmE50_func(Da2, Da21, Di2, Di21):

    return sqrt(((((Da21 + Da2) / 2)** 2) + (((Di21 + Di2) / 2) ** 2)) / 2)

def tangamma_func(Dm1, Dm2, w):

    return (Dm2 - Dm1) / (2 * w)

def cu_func(Dhtp, u2, cu2, u21):

    return (Dhtp + (u2 * cu2)) / u21

def cm21_func(cax21, gamma):

    return sqrt((cax21 ** 2) / (1 - (sin(gamma) ** 2)))

def c21_func(cax21, cu21, cr21):

    return sqrt((cax21 ** 2) + (cu21 ** 2) + (cr21 ** 2))

def c3_comp_func(Dhtp, c2, c21, Rct):

    return sqrt(2 * Dhtp + (c2 ** 2) - (2 * Dhtp - ((c21 ** 2) - (c2 ** 2))) / Rct)

def c3_turb_func(Dhtp, c2, c21, Rct):

    return sqrt((2 * Dhtp * (Rct - 1) + Rct * (c2 ** 2) - (c21 ** 2)) / (Rct - 1))

def mainflow(u1x, u2x, u1y, u2y, u1z, u2z):

    return sqrt(((u1x + u2x) ** 2) + ((u1y + u2y) ** 2) + ((u1z + u2z) ** 2)) / 2

def Fc(mp, cu1, cu2):

    return mp * (cu2 - cu1)

def Fax(mp, cax1, cax2, A1, A2, p1, p2):

    return mp * cax2 - mp * cax1 + A2 * p2 - A1 * p1 + (A1 - A2) * (1 / 2) * (p2 + p1)

def Z_func(DmE, t):

    return (pi * DmE) / t

def Zweifel(t, w, rho1, rho2, cm1, cu1, c2, cu2):

    return 2 * (t / w) * (rho1 / rho2) * (cm1 / c2) * (abs(cu2 - cu1) / c2)

#----------cooling

def epsilon_func(Thg, Tm, Tkl1):

    return (Thg - Tm) / (Thg - Tkl1)

def mpf_func(epsilon, nc):

    return epsilon / (nc * (1 - epsilon))

def Re_func(rho, c, s):

    return (rho * c * s) / y

def Nu_func(Re):

    return (0.037 * (Re ** 0.8) - 850) * Pr ** (1 / 3)

def alpha_hg_func(Nu, s):

    return (Nu * lambda_hg) / s

def A_net_func(z, s, Da, Di, w):

    return 2.6 * z * s * ((1 / 2) * (Da - Di)) + pi * w * (Da - Di)

def mpc_func(mpf, alpha_hg, A_net, Tc):
 
    return (mpf * alpha_hg * A_net) / cps(Tc, 0)