from cp import *
from modules import *
from side import *

class nozzle_class:

    def __init__(self, mps, mp, Tt5m, pt5m, Tt13, pt13, betaB, Pin, Pon, M0, T0, p0, bypass, Hu, nm, status, F, aF, uoC, uoB, name):

        print(name)

        print()

        print("mps:", mps[0], "  mp:", mp, "  Tt5m:", Tt5m, "  pt5m:", pt5m, "  Tt13:", Tt13, "  pt13:", pt13, "  betaB:", betaB, "  Pin:", Pin, "  Pon:", Pon)
        print("M0:", M0, "  T0:", T0, "  p0:", p0, "  bypass:", bypass, "  Hu:", Hu, "  nm:", nm, "  x:", status, "  F:", F, "  aF:", aF)

        print()

        mps9 = mps[0]
        self.mps9 = mps9
        print("mps9:      " + str(mps9))

        mp9 = mps9 * mp
        self.mp9 = mp9
        print("mp9:       " + str(mp9))

        mps19 = bypass
        self.mps19 = mps19
        print("mps19:     " + str(mps19))

        mp19 = mp * mps19
        self.mp19 = mp19
        print("mp19:      " + str(mp19))

        pt9 = pt8 = pt5m * Pin
        self.pt9 = pt9
        print("pt9:       " + str(pt9))

        Tt9 = Tt5m
        self.Tt9 = Tt9
        print("Tt9:       " + str(Tt9))

        pt19 = pt18 = pt13 * Pon
        self.pt19 = pt19
        print("pt19:      " + str(pt19))

        Tt19 = Tt13
        self.Tt19 = Tt19
        print("Tt19:      " + str(Tt19))

        if status == 0:

            p9 = p0
            self.p9 = p9
            print("p9:        " + str(p9))

            T9 = Tts_func_It(Tt9, p9/pt9, betaB, R)
            self.T9 = T9
            print("T9:        " + str(T9))

            M9 = M(Tt9, T9, k(cpd(Tt9, T9, betaB)))
            self.M9 = M9
            print("M9:        " + str(M9))

            a9 = a(R, T9, betaB)
            self.a9 = a9
            print("a9:        " + str(a9))

            c9 = sqrt(2 * cpd(Tt9, T9, betaB) * (Tt9 - T9))
            self.c9 = c9
            print("c9:        " + str(c9))

            p19 = p0
            self.p19 = p19
            print("p19:       " + str(p19))

            T19 = Tts_func_It(Tt19, p19/pt19, 0, R)
            self.T19 = T19
            print("T19:       " + str(T19))

            M19 = M(Tt19, T19, k(cpd(Tt19, T19, betaB)))
            self.M19 = M19
            print("M19:       " + str(M19))

            a19 = a(R, T19, betaB)
            self.a19 = a19
            print("a19:       " + str(a19))

            c19 = sqrt(2 * cpd(Tt19, T19, 0) * (Tt19 - T19))
            self.c19 = c19
            print("c19:       " + str(c19))

        if status == 1:

            T9 = T_func_It(Tt9, 1, betaB)
            self.T9 = T9
            print("T9:        " + str(T9))

            p9 = pt9 * Pbp_func(Tt9, T9, M0) ** (-1)
            self.p9 = p9
            print("p9:        " + str(p9))

            M9 = M(Tt9, T9, k(cpd(Tt9, T9, betaB)))
            self.M9 = M9
            print("M9:        " + str(M9))

            a9 = a(R, T9, betaB)
            self.a9 = a9
            print("a9:        " + str(a9))

            c9 = sqrt(2 * cpd(Tt9, T9, betaB) * (Tt9 - T9))
            self.c9 = c9
            print("c9:        " + str(c9))

            p19 = p0
            self.p19 = p19
            print("p19:       " + str(p19))

            T19 = T_func_It(Tt19, 1, 0)
            self.T19 = T19
            print("T19:       " + str(T19))

            p19 = pt19 * Pbp_func(Tt19, T19, M0) ** (-1)
            self.p19 = p19
            print("p19:       " + str(p19))

            M19 = M(Tt19, T19, k(cpd(Tt19, T19, 0)))
            self.M19 = M19
            print("M19:       " + str(M19))

            a19 = a(R, T19, 0)
            self.a19 = a19
            print("a19:       " + str(a19))

            c19 = sqrt(2 * cpd(Tt19, T19, 0) * (Tt19 - T19))
            self.c19 = c19
            print("c19:       " + str(c19))

        #----------core

        #----------stat9

        rho9 = rho(p9, T9)
        self.rho9 = rho9
        print("rho9:      " + str(rho9))

        Vp9 = Vp(mp9, rho9)
        self.Vp9 = Vp9
        print("Vp9:       " + str(Vp9))

        #----------vel9

        a9 = a(R, T9, betaB)
        self.a9 = a9
        print("a9:        " + str(a9))

        c9 = c(M9, a9)
        self.c1 = c9
        print("c1:        " + str(c9))
        
        cax9 = cax(c9, 0)
        self.cax9 = cax9
        print("cax9:      " + str(cax9))

        Aax9 = Aax(Vp9, cax9)
        self.Aax9 = Aax9
        print("Aax9:      " + str(Aax9))

        #----------entro9

        cps9 = cps(T9, betaB)
        self.cps9 = cps9
        print("cps9:      " + str(cps9))

        EF09 = EF(T0, T9, betaB)
        self.EF09 = EF09
        print("EF09:      " + str(EF09))

        Ds09 = Ds(EF09, p0, p9)
        self.Ds09 = Ds09
        print("Ds09:      " + str(Ds09))

        Dh09 = Dh(T0, T9, betaB)
        self.Dh09 = Dh09
        print("Dh09:      " + str(Dh09))

        #----------nozzle

        a0 = a(R, T0, 0)
        self.a0 = a0
        print("a0:        " + str(a0))

        c0 = c(M0, a0)
        self.c0 = c0
        print("c0:        " + str(c0))

        fsI = fs(c0, c9, betaB)
        self.fsI = fsI
        print("fsI:       " + str(fsI))

        fsII = fs(c0, c19, 0)
        self.fsII = fsII
        print("fsII:      " + str(fsII))

        fsEng = fsEng_func(fsI, fsII, bypass)
        self.fsEng = fsEng
        print("fsEng:     " + str(fsEng))

        Feng = fsEng * (mp * (1 + bypass))
        self.Feng = Feng
        print("Feng:      " + str(Feng))

        bsEng = bsEng_func(betaB, bypass, fsEng)
        self.bsEng = bsEng
        print("bsEng:     " + str(bsEng))

        mp1 = mp1_func(F, fsEng)
        self.mp1 = mp1
        print("mp1:       " + str(mp1))

        mp2 = mp2_func(mp1, bypass)
        self.mp2 = mp2
        print("mp2:       " + str(mp2))

        mp12 = mp12_func(mp2, bypass)
        self.mp12 = mp12
        print("mp12:      " + str(mp12))

        aFlow = aFlow_func(c0, c9, c19, betaB, bypass)
        self.aFlow = aFlow
        print("aFlow:     " + str(aFlow))

        nth = nth_func(aFlow, betaB, Hu)
        self.nth = nth
        print("nth:       " + str(nth))

        nP = nP_func(c0, c9, c19, betaB, bypass, aFlow)
        self.nP = nP
        print("np:        " + str(nP))

        nCore = nCore_func(c0, Tt9, T9, betaB, Hu, bypass, nm, aF)
        self.nCore = nCore
        print("nCore:     " + str(nCore))

        #----------entroC

        cpst9 = cps(Tt9, betaB)
        self.cpst9 = cpst9
        print("cpst9:     " + str(cpst9))

        EF0t9 = EF(T0, Tt9, betaB)
        self.EF0t9 = EF0t9
        print("EF0t9:     " + str(EF0t9))

        Ds0t9 = Ds(EF0t9, p0, pt9)
        self.Ds0t9 = Ds0t9
        print("Ds0t9:     " + str(Ds0t9))

        Dh0t9 = Dh(T0, Tt9, betaB)
        self.Dh0t9 = Dh0t9
        print("Dh0t9:     " + str(Dh0t9))

        #----------entroB

        cpst19 = cps(Tt19, 0)
        self.cpst19 = cpst19
        print("cpst19:    " + str(cpst19))

        EF0t19 = EF(T0, Tt19, 0)
        self.EF0t19 = EF0t19
        print("EF0t19:    " + str(EF0t19))

        Ds0t19 = Ds(EF0t19, p0, pt19)
        self.Ds0t19 = Ds0t19
        print("Ds0t19:    " + str(Ds0t19))

        Dh0t19 = Dh(T0, Tt19, 0)
        self.Dh0t19 = Dh0t19
        print("Dh0t19:    " + str(Dh0t19))
        
        print("----------------------------------------")

        self.ListC1 = []

        self.ListC2 = []

        self.ListB1 = []

        self.ListB2 = []

        self.ListC1 = [uoC.mps2, uoC.mp2, uoC.Ttm, "ct", "ct", "ct", uoC.ptm, "ct", "ct", uoC.cpstm, uoC.EF0tm, uoC.Dh0tm, uoC.betam, "ct", self.c9]

        self.ListC2 = [self.mps9, self.mp9, self.Tt9, "ct", self.T9, Pin, self.pt9, "ct", "ct", self.cpst9, self.EF0t9, self.Dh0t9, betaB, "ct", self.c9]

        self.ListB1 = [uoB.mps2, uoB.mp2, uoB.Tt2, "ct", "ct", "ct", uoB.pt2, "ct", "ct", uoB.cpst2, uoB.EF0t2, uoB.Dh0t2, 0, "ct", self.c19]

        self.ListB2 = [self.mps19, self.mp19, self.Tt19, "ct", self.T19, Pon, self.pt19, "ct", "ct", self.cpst19, self.EF0t19, self.Dh0t19, 0, "ct", self.c19]

#test = nozzle_class([0.022065], 1, 925.30, 209391.54, 337.52, 168629.19, 0.022065, 0.99, 0.98, 0, 288.15, 101325, 4.78, 43124000, 0.995, 1, 61570, 49597)