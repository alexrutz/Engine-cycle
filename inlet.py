from modules import *
from side import *
from cp import *

class inlet_class:

    def __init__(self, mps, mp, T0, p0, M0, PI, P1, P2, Cr, bypass, useAM, Aax1, Aax2, M1, M2, alpha1, alpha2, name):

        print(name)

        print()

        print("mps: ", mps[0], "  mp: ", mp, "  T0: ", T0, "  p0:", p0, "  M0:", M0, "  PI:", PI, "  P1:", P1, "  P2:", P2, "  Cr:",Cr, "  bypass:", bypass)
        print("useAM:", useAM, "  A1:", Aax1, "  A2:", Aax2, "  M1:", M1, "  M2:", M2, "  alpha1:", alpha1, "  alpha2: ", alpha2)

        print()

        mps1 = mps[0] + bypass
        self.mps1 = mps1
        print("mps1:      " + str(mps1))

        mp1 = mp * mps1
        self.mps1 = mps1
        print("mps1:      " + str(mp1))

        #----------downstream

        mps2 = mps[0] + bypass
        self.mps2 = mps2
        print("mps2:      " + str(mps2))

        mp2 = mp * mps2
        self.mp2 = mp2
        print("mp2:       " + str(mp2))

        #----------inlet

        if M0 == 0:

            Tt0 = Tt1 = Tt2 = T0
            self.Tt0 = self.Tt1 = self.Tt2 = T0
            print("Tt2:       " + str(Tt2))

            Pbp = 1
            self.Pbp = Pbp
            print("Pbp:       " + str(Pbp))

            pt0 = pt1 = p0
            self.pt0 = pt0
            print("Pt0/pt1:   " + str(pt0))

            pt2 = pt1 * PI
            self.pt2 = pt2
            print("pt2:       " + str(pt2))

        else:

            Tt0 = Tt1 = Tt2 = Tt_func_It(T0, M0, 0)
            self.Tt0 = self.Tt1 = self.Tt2 = Tt0
            print("Tt2:       " + str(Tt0))

            Pbp = Pbp_func(Tt0, T0, M0)
            self.Pbp = Pbp
            print("Pbp:       " + str(Pbp))

            pt0 = pt1 = p0 * Pbp
            self.pt0 = pt0
            print("Pt0/pt1:   " + str(pt0))

            pt2 = pt1 * PI
            self.pt2 = pt2
            print("pt2:       " + str(pt2))

        #----------entro0t1

        cpst1 = cps(Tt1, 0)
        self.cpst1 = cpst1
        print("cpst1:     " + str(cpst1))

        EF0t1 = EF(T0, Tt1, 0)
        self.EF0t1 = EF0t1
        print("EF0t1:     " + str(EF0t1))

        Ds0t1 = Ds(EF0t1, p0, pt1)
        self.Ds0t1 = Ds0t1
        print("Ds0t1:     " + str(Ds0t1))

        Dh0t1 = Dh(T0, Tt1, 0)
        self.Dh0t1 = Dh0t1
        print("Dh0t1:     " + str(Dh0t1))

        #----------entro0t2

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

        #----------T1, cpdt11

        if useAM == "A":

            M1 = Ma_It(mp1, Tt1, pt1, Aax1, 0)
            self.M1 = M1
            print("M1:        " + str(M1))

        T1 = T_func_It(Tt1, M1, 0)
        self.T1 = T1
        print("T1:        " + str(T1))

        cpdt11 = cpd(Tt1, T1, 0)
        self.cpdt11 = cpdt11
        print("cpdt11:    " + str(cpdt11))

        #----------stat1

        p1 = p(pt1, T1, Tt1, 0)
        self.p1 = p1
        print("p1:        " + str(p1))

        rho1 = rho(p1, T1)
        self.rho1 = rho1
        print("rho1:      " + str(rho1))

        Vp1 = Vp(mp1, rho1)
        self.Vp1 = Vp1
        print("Vp1:       " + str(Vp1))

        #----------vel1

        a1 = a(R, T1, 0)
        self.a1 = a1
        print("a1:        " + str(a1))

        c1 = c(M1, a1)
        self.c1 = c1
        print("c1:        " + str(c1))
        
        cax1 = cax(c1, alpha1)
        self.cax1 = cax1
        print("cax1:      " + str(cax1))

        #----------entro1

        cps1 = cps(T1, 0)
        self.cps1 = cps1
        print("cps1:      " + str(cps1))

        EF01 = EF(T0, T1, 0)
        self.EF01 = EF01
        print("EF01:      " + str(EF01))

        Ds01 = Ds(EF01, p0, p1)
        self.Ds01 = Ds01
        print("Ds01:      " + str(Ds01))

        Dh01 = Dh(T0, T1, 0)
        self.Dh01 = Dh01
        print("Dh01:      " + str(Dh01))

        #----------T2, cpdt22

        if useAM == "A":

            M2 = Ma_It(mp2, Tt2, pt2, Aax2, 0)
            self.M2 = M2
            print("M2:        " + str(M2))
        
        T2 = T_func_It(Tt2, M2, 0)
        self.T2 = T2
        print("T2:        " + str(T2))

        cpdt22 = cpd(Tt2, T2, 0)
        self.cpdt22 = cpdt22
        print("cpdt22:    " + str(cpdt22))

        #----------stat2

        p2 = p(pt2, T2, Tt2, 0)
        self.p2 = p2
        print("p2:        " + str(p2))

        rho2 = rho(p2, T2)
        self.rho2 = rho2
        print("rho2:      " + str(rho2))

        Vp2 = Vp(mp2, rho2)
        self.Vp2 = Vp2
        print("Vp2:       " + str(Vp2))

        #----------vel2

        a2 = a(R, T2, 0)
        self.a2 = a2
        print("a2:        " + str(a2))

        c2 = c(M2, a2)
        self.c2 = c2
        print("c2:        " + str(c2))
        
        cax2 = cax(c2, alpha2)
        self.cax2 = cax2
        print("cax2:      " + str(cax2))

        if useAM == "M":

            Aax2 = Aax(Vp2, cax2)
            self.Aax2 = Aax2
            print("Aax2:      " + str(Aax2))

        #----------entro02

        cps2 = cps(T2, 0)
        self.cps2 = cps2
        print("cps2:      " + str(cps2))

        EF02 = EF(T0, T2, 0)
        self.EF02 = EF02
        print("EF0t2:     " + str(EF0t2))

        Ds02 = Ds(EF02, p0, p2)
        self.Ds02 = Ds02
        print("Ds0t2:     " + str(Ds02))

        Dh02 = Dh(T0, T2, 0)
        self.Dh02 = Dh02
        print("Dh02:      " + str(Dh02))

        #----------pressures
            
        PC = PC_func(Cr, PI, Pbp)
        self.PC = PC
        print("PC:        " + str(PC))

        PLPC = PLPC_func(PC, P1)
        self.PLPC = PLPC
        print("PLPC:      " + str(PLPC))

        PHPC = PHPC_func(PC, P1)
        self.PHPC = PHPC
        print("PHPC:      " + str(PHPC))

        PHPC1 = PHPC1_func(PHPC, P2)
        self.PHPC1 = PHPC1
        print("PHPC1:     " + str(PHPC1))

        PHPC2 = PHPC2_func(PHPC, P2)
        self.PHPC2 = PHPC2
        print("PHPC2:     " + str(PHPC2))

        print("----------------------------------------")

        self.List = []

        self.List.extend([mps2, mp2, Tt2, "ct", "ct", PI, pt2])