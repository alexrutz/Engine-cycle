from side import *
from cp import *
from modules import *

class burner_class:
    
    def __init__(self, mps, mp, Tt3, pt3, Pb, nB, Hu, Tt4, y_401, y_41, y_5, OTDF, useAM, Aax3, Aax4, M3, M4, alpha3, alpha4, Dm3, Dm4, name):

        print(name)

        print()

        print("mps: ", mps[0], "  mp:", mp, "  Tt3:", Tt3, "  pt3:", pt3, "  Pb:", Pb, "  nB:", nB, "  Hu:", Hu, "  Tt4:", Tt4, "  y_401:", y_401, "  y_41:", y_41, "  y_5:", y_5, "  OTDF:", OTDF)
        print("useAM:", useAM, "  A3:", Aax3, "  A4:", Aax4, "  M3:", M3, "  M4:", M4, "  alpha3:", alpha3, "  alpha4:", alpha4, "  Dm3:", Dm3, "  Dm4:", Dm4)

        print()
        
        y_kl = y_401 + y_41 + y_5
        self.Aax3 = Aax3
        self.Aax4 = Aax4
        self.M3 = M3
        self.M4 = M4

        self.mps1 = mps[0]
        print("mps1:      " + str(mps[0]))

        mp1 = mp * mps[0]
        self.mp1 = mp1
        print("mp1:       " + str(mp1))

        betaB = beta_func(T0, Tt3, Tt4, y_kl, Hu, nB)
        self.betaB = betaB
        print("betaBK:    " + str(betaB))

        beta4 = beta4_func(betaB, y_kl)
        self.beta4 = beta4
        print("beta4:     " + str(beta4))

        cpd0t4 = cpd(T0, Tt4, beta4)
        self.cpd0t4 = cpd0t4
        print("cpd0t4:    " + str(cpd0t4))

        #----------downstream

        mps[0] = mps[0] + betaB
        self.mps2 = mps[0]
        print("mps2:      " + str(mps[0]))

        mp2 = mp * mps[0]
        self.mp2 = mp2
        print("mp2:       " + str(mp2))

        #----------flows

        beta_4 = betaB /(1 - y_kl)
        self.beta_4 = beta4
        print("beta_4:    " + str(beta_4))

        beta_401m = betaB /(1 - y_kl + y_401)
        self.beta_401m = beta_401m
        print("beta_401m: " + str(beta_401m))

        beta_41m = betaB /(1 - y_kl + y_401 + y_41)
        self.beta_41m = beta_41m
        print("beta_41m:  " + str(beta_41m))

        beta_5m = betaB /(1 - y_kl + y_401 + y_41 + y_5)
        self.beta_5m = beta_5m
        print("beta_5m:   "+ str(beta_5m))

        Tt4max = Tt4max_func(OTDF, Tt3, Tt4)
        self.Tt4max = Tt4max
        print("Tt4max:    " + str(Tt4max))
        
        pt4 = pt3 * Pb
        self.pt4 = pt4
        print("pt4:       " + str(pt4))

        mp3 = massflows_func(mp1, y_kl, betaB)[0]
        self.mp3 = mp3
        print("mp3:       " + str(mp3))

        mp4 = massflows_func(mp1, y_kl, betaB)[1]
        self.mp4 = mp4
        print("mp4:       " + str(mp4))

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

        #----------entrot4

        cpst4 = cps(Tt4, beta4)
        self.cpst4 = cpst4
        print("cpst4:     " + str(cpst4))

        EF0t4 = EF(T0, Tt4, beta4)
        self.EF0t4 = EF0t4
        print("EF0t4:     " + str(EF0t4))

        Ds0t4 = Ds(EF0t4, p0, pt4)
        self.Ds0t4 = Ds0t4
        print("Ds0t4:     " + str(Ds0t4))

        Dh0t4 = Dh(T0, Tt4, beta4)
        self.Dh0t4 = Dh0t4
        print("Dh0t4:     " + str(Dh0t4))

        #----------T3, cpd3

        if useAM == "A":

            M3 = Ma_It(mp3, Tt3, pt3, Aax3, 0)
            self.M3 = M3
            print("M3:        " + str(M3))
        
        T3 = T_func_It(Tt3, M3, 0)
        self.T3 = T3
        print("T3:        " + str(T3))

        cpdt33 = cpd(Tt3, T3, 0)
        self.cpdt33 = cpdt33
        print("cpdt33:    " + str(cpdt33))

        #----------stat3

        p3 = p(pt3, T3, Tt3, 0)
        self.p3 = p3
        print("p3:        " + str(p3))

        rho3 = rho(p3, T3)
        self.rho3 = rho3
        print("rho3:      " + str(rho3))

        Vp3 = Vp(mp3, rho3)
        self.Vp3 = Vp3
        print("Vp3:       " + str(Vp3))

        #----------vel3

        a3 = a(R, T3, 0)
        self.a3 = a3
        print("a3:        " + str(a3))

        c3 = c(M3, a3)
        self.c3 = c3
        print("c3:        " + str(c3))
        
        cax3 = cax(c3, alpha3)
        self.cax3 = cax3
        print("cax3:      " + str(cax3))

        if useAM == "M":

            Aax3 = Aax(Vp3, cax3)
            self.Aax3 = Aax3
            print("Aax3:      " + str(Aax3))

        #----------entro3

        cps3 = cps(T3, 0)
        self.cps3 = cps3
        print("cps3:      " + str(cps3))

        EF03 = EF(T0, T3, 0)
        self.EF03 = EF03
        print("EF03:      " + str(EF03))

        Ds03 = Ds(EF03, p0, p3)
        self.Ds03 = Ds03
        print("Ds03:      " + str(Ds03))

        Dh03 = Dh(T0, T3, 0)
        self.Dh03 = Dh03
        print("Dh03:      " + str(Dh03))

        #----------T4, cpdt44

        if useAM == "A":

            M4 = Ma_It(mp4, Tt4, pt4, Aax4, beta4)
            self.M4 = M4
            print("M4:        " + str(M4))
        
        T4 = T_func_It(Tt4, M4, beta4)
        self.T4 = T4
        print("T4:        " + str(T4))

        cpdt44 = cpd(Tt4, T4, beta4)
        self.cpdt22 = cpdt44
        print("cpdt44:    " + str(cpdt44))

        #----------stat4

        p4 = p(pt4, T4, Tt4, beta4)
        self.p4 = p4
        print("p4:        " + str(p4))

        rho4 = rho(p4, T4)
        self.rho4 = rho4
        print("rho4:      " + str(rho4))

        Vp4 = Vp(mp4, rho4)
        self.Vp4 = Vp4
        print("Vp4:       " + str(Vp4))

        #----------vel4

        a4 = sqrt(k(cpd(Tt4, T4, beta4)) * R * T4)
        self.a4 = a4
        print("a4:        " + str(a4))

        c4 = c(M4, a4)
        self.c4 = c4
        print("c4:        " + str(c4))
        
        cax4 = cax(c4, alpha4)
        self.cax4 = cax4
        print("cax4:      " + str(cax4))

        if useAM == "M":

            Aax4 = Aax(Vp4, cax4)
            self.Aax4 = Aax4
            print("Aax4:      " + str(Aax4))

        #----------entro4

        cps4 = cps(T4, beta4)
        self.cps4 = cps4
        print("cps4:      " + str(cps4))

        EF04 = EF(T0, T4, beta4)
        self.EF04 = EF04
        print("EF04:      " + str(EF04))

        Ds04 = Ds(EF04, p0, p4)
        self.Ds04 = Ds04
        print("Ds04:      " + str(Ds04))

        Dh04 = Dh(T0, T4, beta4)
        self.Dh04 = Dh04
        print("Dh04:      " + str(Dh04))

        #----------critical quantities

        T4crit = T4crit_func_It(Tt4, beta4)
        self.T4crit = T4crit
        print("T4crit:    " + str(T4crit))

        cpd4crit = cpd(Tt4, T4crit, beta4)
        self.cpd4crit = cpd4crit
        print("cpd4crit:  " + str(cpd4crit))

        k4 = k(cpd4crit)
        self.k4 = k4
        print("k4:        " + str(k4))

        mpred4 = mpred(Tt4, beta4)
        self.mpred4 = mpred4
        print("mpred4:    " + str(mpred4))

        A4crit = A4crit_func(mpred4, pt4, Tt4, mp4)
        self.A4crit = A4crit
        print("A4crit:    " + str(A4crit))

        Pt1t2 = P(mp2, Dh(Tt3, Tt4, beta4)) / 1000
        self.Pt1t2 = Pt1t2
        print("Pt1t2:     " + str(Pt1t2))

        #----------geometry

        Da3 = Da(Dm3, Aax3)
        self.Da3 = Da3
        print("Da3:       " + str(Da3))

        Da4 = Da(Dm4, Aax4)
        self.Da4 = Da4
        print("Da4:       " + str(Da4))

        Di3 = Di(Dm3, Aax3)
        self.Di3 = Di3
        print("Di3:       " + str(Di3))

        Di4 = Di(Dm4, Aax4)
        self.Di4 = Di4
        print("Di4:       " + str(Di4))
        
        print("----------------------------------------")

        self.List1 = []

        self.List2 = []

        self.List1.extend([self.mps1, self.mp1, Tt3, "ct", self.T3, "ct", pt3, "ct", "ct", self.cpst3, self.EF0t3, self.Dh0t3, 0, self.M3, self.c3, self.p3, self.rho3, self.Aax3, Dm3, self.Di3, self.Da3])

        self.List2.extend([self.mps2, self.mp2, Tt4, "ct", self.T4, Pb, self.pt4, Dh(Tt3, Tt4, beta4), self.Pt1t2, self.cpst4, self.EF0t4, self.Dh0t4, self.beta4, self.M4, self.c4, self.p4, self.rho4, self.Aax4, Dm4, self.Di4, self.Da4])

#test = burner_class([1], 1, 800, 2800000, 0.948, 0.995, 43124000, 1700, 0.08, 0.05, 0.05, 0.15)