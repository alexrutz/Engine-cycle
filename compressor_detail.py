from math import sin, cos, tan, atan, pi
from side import *
from modules import *

class compressor_detail:

    def __init__(self, mp2, Tt2, pt2, Pic, dTt, N, c2, cax21cax2, Dhtp, Rct, nvp, alpha2, alpha3, Dm2, Dm3, Di2, Di3, Da2, Da3, arR, arS, Z, tsR, tsS, name):

        print(name)

        print()

        print("mp2:", mp2, "  Tt2:", Tt2, "  pt2:", pt2, "  Pic:", Pic, "  dTt:", dTt, "  N:", N, "  c2:", c2, "  cax21cax2:", cax21cax2, "  Dhtp:", Dhtp, "  Rct:", Rct, "  nvp:", nvp)
        print("alpha2:", alpha2, "  alpha3:", alpha3, "  Dm2:", Dm2, "  Dm3:", Dm3, "  Di2:", Di2, "  Di3:", Di3, "  Da2:", Da2, "  Da3:", Da3, "  arR:", arR, "  arS", arS, "  Z:", Z)

        print()

        self.inputlist = [mp2, N, Dm2, Dm3, Di2, Di3, Da2, Da3, alpha2, alpha3, Tt2, pt2, Pic, nvp, c2, R, cax21cax2, arR, arS]

        self.alpha2 = alpha2
        self.alpha3 = alpha3

        N = N / 60
        self.N = N
        print("N:         " + str(N))

        Tt3 = Tt21 = Tt2 + dTt
        self.Tt21 = Tt21
        self.Tt3 = Tt3
        print("Tt21/3:    " + str(Tt21))

        pt21 = pt2 * Pic
        self.pt21 = pt21
        print("pt21:      " + str(pt21))

        pt3 = pt21
        self.pt3 = pt3
        print("pt3:       " + str(pt3))

        Di21 = D21_It(Dm2, Dm3, Da2, Da3, Di2, Di3, arR, arS)[0]
        self.Di21 = Di21
        print("Di21:      " + str(Di21))

        Dm21 = D21_It(Dm2, Dm3, Da2, Da3, Di2, Di3, arR, arS)[1]
        self.Dm21 = Dm21
        print("Dm21:      " + str(Dm21))

        Da21 = D21_It(Dm2, Dm3, Da2, Da3, Di2, Di3, arR, arS)[2]
        self.Da21 = Da21
        print("Da21:      " + str(Da21))

        wR = D21_It(Dm2, Dm3, Da2, Da3, Di2, Di3, arR, arS)[3]
        self.wR = wR
        print("wR:        " + str(wR))

        wS = D21_It(Dm2, Dm3, Da2, Da3, Di2, Di3, arR, arS)[4]
        self.wS = wS
        print("wS:        " + str(wS))

        wSt = D21_It(Dm2, Dm3, Da2, Da3, Di2, Di3, arR, arS)[5]
        self.wSt = wSt
        print("wSt:       " + str(wSt))

        tangammaR = tangamma_func(Dm2, Dm21, wR)
        self.tangammaR = tangammaR
        print("tangammaR: " + str(tangammaR))

        tangammaS = tangamma_func(Dm21, Dm3, wS)
        self.tangammaS = tangammaS
        print("tangammaS: " + str(tangammaS))

        tangammaSt = tangamma_func(Dm2, Dm3, wSt)
        self.tangammaSt = tangammaSt
        print("tangammaSt:" + str(tangammaSt))

        DmE2 = DmE_func(Da2, Di2)
        self.DmE2 = DmE2
        print("DmE2:      " + str(DmE2))

        DmER = DmE50_func(Da2, Da21, Di2, Di21)
        self.DmER = DmER
        print("DmER:      " + str(DmER))

        DmE21 = DmE_func(Da21, Di21)
        self.DmE21 = DmE21
        print("DmE21:     " + str(DmE21))

        DmE3 = DmE_func(Da3, Di3)
        self.DmE3 = DmE3
        print("DmE3:      " + str(DmE3))

        u2 = pi * DmE2 * N
        self.u2 = u2
        print("u2:        " + str(u2))

        umER = pi * DmER * N
        self.umER = umER
        print("umER:      " + str(umER))

        u21 = pi * DmE21 * N
        self.u21 = u21
        print("u21:       " + str(u21))

        u3 = pi * DmE3 * N
        self.u3 = u3
        print("u3:        " + str(u3))

        tangammaESt = tangamma_func(DmE2, DmE3, wSt)
        self.tangammaESt = tangammaESt
        print("tangammaESt:" + str(tangammaESt))

        tangammaER = tangamma_func(DmE2, DmE21, wR)
        self.tangammaER = tangammaER
        print("tangammaER:" + str(tangammaER))

        tangammaES = tangamma_func(DmE21, DmE3, wS)
        self.tangammaES = tangammaES
        print("tangammaES:" + str(tangammaES))

        #----------2

        print("c2:        " + str(c2))

        cm2 = c2 * cos(alpha2)
        self.cm2 = cm2
        print("cm2:       " + str(cm2))

        cu2 = cm2 * tan(alpha2)
        self.cu2 = cu2
        print("cu2:       " + str(cu2))

        cr2 = cm2 * sin(atan(tangammaESt))
        self.cr2 = cr2
        print("cr2:       " + str(cr2))

        cax2 = sqrt((cm2 ** 2) - (cr2 ** 2))
        self.cax2 = cax2
        print("cax2:      " + str(cax2))

        #----------21

        cu21 = cu_func(Dhtp, u2, cu2, u21)
        self.cu21 = cu21
        print("cu21:      " + str(cu21))

        cax21 = cax21cax2 * cax2
        self.cax21 = cax21
        print("cax21:     " + str(cax21))

        cm21 = cm21_func(cax21, atan(tangammaESt))
        self.cm21 = cm21
        print("cm21:      " + str(cm21))

        cr21 = cm21 * tangammaESt
        self.cr21 = cr21
        print("cr21:      " + str(cr21))

        c21 = c21_func(cax21, cu21, cr21)
        self.c21 = c21
        print("c21:       " + str(c21))

        alpha21 = atan(cu21 / cm21)
        self.alpha21 = alpha21
        print("alpha21:   " + str((alpha21 * 180) / pi))

        #----------3

        print(Dhtp, c2, c21, Rct)
        print(2 * Dhtp + (c2 ** 2))
        print(2 * Dhtp - ((c21 ** 2) - (c2 ** 2)))
        print(2 * Dhtp + (c2 ** 2) - (2 * Dhtp - ((c21 ** 2) - (c2 ** 2))) / Rct)
        c3 = c3_comp_func(Dhtp, c2, c21, Rct)
        self.c3 = c3
        print("c3:        " + str(c3))

        cm3 = c3 * cos(alpha3)
        self.cm3 = cm3
        print("cm3:       " + str(cm3))

        cr3 = cm3 * sin(atan(tangammaESt))
        self.cr3 = cr3
        print("cr3:       " + str(cr3))

        cax3 = sqrt((cm3 ** 2) - (cr3 ** 2))
        self.cax3 = cax3
        print("cax3:      " + str(cax3))

        cu3 = c3 * sin(alpha3)
        self.cu3 = cu3
        print("cu3:       " + str(cu3))

        #----------relsys

        wax2 = cax2
        self.wax2 = wax2
        print("wax2:      " + str(wax2))

        wu2 = cu2 - u2
        self.wu2 = wu2
        print("wu2:       " + str(wu2))

        wr2 = cr2
        self.wr2 = wr2
        print("wr2:       " + str(wr2))

        w2 = sqrt((wax2 ** 2) + (wr2 ** 2) + (wu2 ** 2))
        self.w2 = w2
        print("w2:        " + str(w2))

        wax21 = cax21
        self.wax21 = wax21
        print("wax21:     " + str(wax21))

        wu21 = cu21 - u21
        self.wu21 = wu21
        print("wu21:      " + str(wu21))

        wr21 = cr21
        self.wr21 = wr21
        print("wr21:      " + str(wr21))

        w21 = sqrt((wax21 ** 2) + (wr21 ** 2) + (wu21 ** 2))
        self.w21 = w21
        print("w21:       " + str(w21))

        w2w21 = w2 / w21
        self.w2w21 = w2w21
        print("w2w21:     " + str(w2w21))

        beta2 = atan(wu2 / cm2)
        self.beta2 = beta2
        print("beta2:     " + str((beta2  * 180) / pi))

        beta21 = atan(wu21 / cm21)
        self.beta21 = beta21
        print("beta21:    " + str((beta21 * 180) / pi))

        #----------contour

        T2 = T_contour_It(Tt2, c2, 0)
        self.T2 = T2
        print("T2:        " + str(T2))

        p2 = p(pt2, T2, Tt2, 0)
        self.p2 = p2
        print("p2:        " + str(p2))

        rho2 = rho(p2, T2)
        self.rho2 = rho2
        print("rho2:      " + str(rho2))

        H2 = mp2 / (rho2 * cax2 * pi * Dm2)
        self.H2 = H2
        print("H2:        " + str(H2))

        Di2 = Dm2 - H2
        self.Di2 = Di2
        print("Di2:       " + str(Di2))

        Da2 = Dm2 + H2
        self.Da2 = Da2
        print("Da2:       " + str(Da2))

        Aax2 = (pi * (Da2 - Di2) ** 2) / 4
        self.Aax2 = Aax2
        print("Aax2:      " + str(Aax2))

        MC2 = c2 / a(R, T2, 0)
        self.MC2 = MC2
        print("MC2:       " + str(MC2))

        MW2 = w2 / a(R, T2, 0)
        self.MW2 = MW2
        print("MW2:       " + str(MW2))

        T21 = T_contour_It(Tt21, c21, 0)
        self.T21 = T21
        print("T21:       " + str(T21))

        p21 = p(pt21, T21, Tt21, 0)
        self.p21 = p21
        print("p21:       " + str(p21))

        rho21 = rho(p21, T21)
        self.rho21 = rho21
        print("rho21:     " + str(rho21))

        H21 = mp2 / (rho21 * cax21 * pi * Dm21)
        self.H21 = H21
        print("H21:       " + str(H21))

        Di21 = Dm21 - H21
        self.Di21 = Di21
        print("Di21:      " + str(Di21))

        Da21 = Dm21 + H21
        self.Da21 = Da21
        print("Da21:      " + str(Da21))

        Aax21 = (pi * (Da21 - Di21) ** 2) / 4
        self.Aax21 = Aax21
        print("Aax21:     " + str(Aax21))

        MC21 = c21 / a(R, T21, 0)
        self.MC21 = MC21
        print("MC21:      " + str(MC21))

        MW21 = w21 / a(R, T21, 0)
        self.MW21 = MW21
        print("MW21:      " + str(MW21))

        T3 = T_contour_It(Tt3, c3, 0)
        self.T3 = T3
        print("T3:        " + str(T3))

        p3 = p(pt3, T3, Tt3, 0)
        self.p3 = p3
        print("p3:        " + str(p3))

        rho3 = rho(p3, T3)
        self.rho3 = rho3
        print("rho3:      " + str(rho3))

        H3 = mp2 / (rho3 * cax3 * pi * Dm3)
        self.H3 = H3
        print("H3:        " + str(H3))

        Di3 = Dm3 - H3
        self.Di3 = Di3
        print("Di3:       " + str(Di3))

        Da3 = Dm3 + H3
        self.Da3 = Da3
        print("Da3:       " + str(Da3))

        Aax3 = (pi * (Da3 - Di3) ** 2) / 4
        self.Aax3 = Aax3
        print("Aax3:      " + str(Aax3))

        MC3 = c3 / a(R, T3, 0)
        self.MC3 = MC3
        print("MC3:       " + str(MC3))

        MFc221 = mainflow(cu2, cu21, cax2, cax21, cr2, cr21)
        self.MFc221 = MFc221
        print("MFc221:    " + str(MFc221))

        MFw221 = mainflow(wu2, wu21, wax2, wax21, wr2, wr21)
        self.MFw221 = MFw221
        print("MFw221:    " + str(MFw221))

        #----------veltri

        """ cu2um = cu2 / umER
        self.cu2um = cu2um
        print("cu2um:     " + str(cu2um))
        
        cu21um = cu21 / umER
        self.cu21um = cu21um
        print("cu21um:    " + str(cu21um))
        
        wu2um = wu2 / umER
        self.wu2um = wu2um
        print("wu2um:     " + str(wu2um))
        
        wu21um = wu21 / umER
        self.wu21um = wu21um
        print("wu21um:    " + str(wu21um))
        
        cax2um = cax2 / umER
        self.cax2um = cax2um
        print("cax2um:    " + str(cax2um))
        
        cax21um = cax21 / umER
        self.cax21um = cax21um
        print("cax21um:   " + str(cax21um))
        
        wax2um = cax2 / umER
        self.wax2um = wax2um
        print("wax2um:    " + str(wax2um))
        
        wax21um = wax21 / umER
        self.wax21um = wax21um
        print("wax21um:   " + str(wax21um)) """

        #----------forces

        Fc221 = Fc(mp2, cu2, cu21)
        self.Fc221 = Fc221
        print("Fc221:     " + str(Fc221))

        Fax221 = Fax(mp2, cax2, cax21, Aax2, Aax21, p2, p21)
        self.Fax221 = Fax221
        print("Fax221:    " + str(Fax221))

        F221 = sqrt((Fc221 ** 2) + (Fax221 ** 2))
        self.F221 = F221
        print("F221:      " + str(F221))

        Fc213 = Fc(mp2, cu21, cu3)
        self.Fc213 = Fc213
        print("Fc213:     " + str(Fc213))

        Fax213 = Fax(mp2, cax21, cax3, Aax21, Aax3, p21, p3)
        self.Fc213 = Fax213
        print("Fc213:     " + str(Fc213))

        F213 = sqrt((Fc213 ** 2) + (Fax213 ** 2))
        self.F213 = F213
        print("F213:      " + str(F213))

        #----------strain

        lambdaR = (1 / 2) * abs(beta21 + beta2)
        self.lambdaR = lambdaR
        print("lambdaR:   " + str(lambdaR))

        bsR = cos(lambdaR)
        self.bsR = bsR
        print("bsR:       " + str(bsR))

        sR = wR / cos(lambdaR)
        self.sR = sR
        print("sR:        " + str(sR))

        tR = sR * tsR
        self.tR = tR
        print("tR:        " + str(tR))

        ZR = Z_func(DmE2, tR)
        self.ZR = ZR
        print("ZR:        " + str(ZR))

        PZ = (Dhtp * mp2) / ZR
        self.PZ = PZ
        print("PZ:        " + str(PZ))

        DeHaller_R = w2 / w21
        self.DeHaller_R = DeHaller_R
        print("DeHaller_R:" + str(DeHaller_R))

        loadFactor_R = (2 * abs(wu2 - wu21)) / ((w2 + w21) / 2)
        self.loadFactor_R = loadFactor_R
        print("loadF_R:   " + str(loadFactor_R))

        Lieblein_R = 1 - (w21 / w2) + abs(wu21 - wu2) / (2 * tsR * w2)
        self.Lieblein_R = Lieblein_R
        print("Lieblein_R:" + str(Lieblein_R))

        lambdaS = (1 / 2) * abs(alpha21 + alpha2)
        self.lambdaS = lambdaS
        print("lambdaS:   " + str(lambdaS))

        bsS = cos(lambdaS)
        self.bsS = bsS
        print("bsS:       " + str(bsS))

        sS = wS / cos(lambdaS)
        self.sS = sS
        print("sS:        " + str(sS))

        tS = sS * tsS
        self.tS = tS
        print("tS:        " + str(tS))

        ZS = Z_func(DmE21, tS)
        self.ZS = ZS
        print("ZS:        " + str(ZS))

        DeHaller_S = c21 / c3
        self.DeHaller_S = DeHaller_S
        print("DeHaller_S:" + str(DeHaller_S))

        loadFactor_S = (2 * abs(cu3 - cu21)) / ((c21 + c3) / 2)
        self.loadFactor_S = loadFactor_S
        print("loadF_S:   " + str(loadFactor_S))

        Lieblein_S = 1 - (c3 / c21) + abs(cu3 - cu21) / (2 * tsS * c3)
        self.Lieblein_S = Lieblein_S
        print("Lieblein_S:" + str(Lieblein_S))

        #----------entrot2

        cpst2 = cps(Tt2, 0)
        self.cpst2 = cpst2
        print("cpst2:     " + str(cpst2))

        EF0t2 = EF(T0, Tt2, 0)
        self.EF0t2 = EF0t2
        print("EF0t2:     " + str(EF0t2))

        Ds0t2 = Ds(EF0t2, p0, pt2)
        self.Ds0t2 = Ds0t2
        print("Ds0t2:     " + str(Ds0t2))

        Dh0t2 = Dh(T0, Tt2, 0)
        self.Dh0t2 = Dh0t2
        print("Dh0t2:     " + str(Dh0t2))

        #----------entrot21

        cpst21 = cps(Tt21, 0)
        self.cpst21 = cpst21
        print("cpst21:    " + str(cpst21))

        EF0t21 = EF(T0, Tt21, 0)
        self.EF0t21 = EF0t21
        print("EF0t21:    " + str(EF0t21))

        Ds0t21 = Ds(EF0t21, p0, pt21)
        self.Ds0t21 = Ds0t21
        print("Ds0t21:    " + str(Ds0t21))

        Dh0t21 = Dh(T0, Tt21, 0)
        self.Dh0t21 = Dh0t21
        print("Dh0t21:    " + str(Dh0t21))

        #----------entrot3

        cpst3 = cps(Tt3, 0)
        self.cpst3 = cpst3
        print("cpst3:     " + str(cpst3))

        EF0t3 = EF(T0, Tt3, 0)
        self.EF0t3 = EF0t3
        print("EF0t3:     " + str(EF0t3))

        Ds0t3 = Ds(EF0t3, p0, pt3)
        self.Ds0t3 = Ds0t3
        print("Ds0t3:     " + str(Ds0t3))

        Dh0t3 = Dh(T0, Tt3, 0)
        self.Dh0t3 = Dh0t3
        print("Dh0t3:     " + str(Dh0t3))
        
        #----------outputlists

        list2 = [pt2, p2, Tt2, T2, umER, c2, MC2, cu2, cax2, cr2, alpha2 * (180 / pi), MFc221, w2, MW2, wu2, wax2, wr2, beta2 * (180 / pi), MFw221]
        list21 = [pt21, p21, Tt21, T21, umER, c21, MC21, cu21, cax21, cr21, alpha21 * (180 / pi), "ct", w21, MW21, wu21, wax2, wr2, beta21 * (180 / pi), "ct"]
        list3 = [pt3, p3, Tt3, T3, umER, c3, MC3, cu3, cax3, cr3, alpha3 * (180 / pi), "ct", "ct", "ct", "ct", "ct", "ct", "ct", "ct"]

        self.st_inputlist = [[mp2, N, Dm2, Dm3, Di2, Di3, Da2, Da3, alpha2 * (180 / pi), alpha3 * (180 / pi), Tt2, pt2, dTt, nvp, c2, Rct, cax21cax2, arR, arS]]
        self.st_outputlist = [[Dm21, Di21, Da21, wR, wS, wSt, tangammaR, tangammaS, tangammaSt, DmE2, DmER, DmE21, DmE3, u2, umER, u21, u3, tangammaER, tangammaES, tangammaESt]]

        self.st_veltri2 = [[0, 0, cax2, 0], [0, umER, cu2, 0]]
        self.st_veltri21 = [[0, 0, cax21, 0], [0, umER, cu21, 0]]

        self.st_canal = [[1 + (2 * Z), 1 + (2 * Z), 2 + (2 * Z), 3 + (2 * Z), 3 + (2 * Z), 3 + (2 * Z), 2 + (2 * Z), 1 + (2 * Z), 1 + (2 * Z), 2 + (2 * Z), 3 + (2 * Z), 2 + (2 * Z), 2 + (2 * Z), 2 + (2 * Z), 2 + (2 * Z), 2 + (2 * Z), 3 + (2 * Z)], [Dm2, Da2, Da21, Da3, Dm3, Di3, Di21, Di2, Dm2, Dm21, Dm3, Dm21, Da21, Dm21, Di21, Dm21, Dm3]]

        self.st_vel = [list2, list21, list3]

        self.st_strain = [[tsR, lambdaR, sR, tR, ZR, PZ, DeHaller_R, loadFactor_R, Lieblein_R], [tsS, lambdaS, sS, tS, ZS, "ct", DeHaller_S, loadFactor_S, Lieblein_S]]

        self.st_T = [Tt2, Tt21, Tt3]
        self.st_s = [Ds0t2, Ds0t21, Ds0t3]
        self.st_p = [pt2, pt21, pt3]
        self.st_v = [1 / rho2, 1 / rho21, 1 / rho3]
        
        print("----------------------------------------")