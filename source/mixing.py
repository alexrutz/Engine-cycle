from side import *
from cp import *
from modules import *

class mixing_class():

    def __init__(self, mps, mp1, Ttus, ptus, Pm, Ttc, betaus, betam, y_c, useAM, Aaxus, Aaxm, Mus, Mm, alphaus, alpham, name):

        print(name)

        print()

        print("mps:", mps[0], "  mp1:", mp1, "  Ttus:", Ttus, "  ptus:", ptus, "  Pm:", Pm, "  Ttc:", Ttc, "  betaus:", betaus, "  betam:", betam, "  y_c:", y_c)
        print("useAM:", useAM, "  Aus:", Aaxus, "  Am:", Aaxm, "  Mus:", Mus, "  Mm:", Mm, "  alphaus:", alphaus, "  alpham:", alpham)

        print()

        self.Aaxus = Aaxus
        self.Aaxm = Aaxm
        self.Mus = Mus
        self.Mm = Mm

        mps1 = mps[0]
        self.mps1 = mps1
        print("mps1:      " + str(mps[0]))

        mp1 = mp1 * mps[0]
        self.mp1 = mp1
        print("mp1:       " + str(mp1))

        self.Ttus = Ttus

        #----------mixing
                
        Ttm = Ttm_func_It(T0, Ttc, Ttus, y_c, mps[0], betaus, betam)
        self.Ttm = Ttm
        print("Ttm:       " + str(Ttm))

        cpd0tm = cpd(T0, Ttm, betam)
        self.cpd0tm = cpd0tm
        print("cpd0tm:    " + str(cpd0tm))

        ptm = ptus * Pm
        self.ptm = ptm
        print("ptm:       " + str(ptm))

        self.betam = betam
        print("betam:     " + str(betam))

        #----------downstream

        mps[0] = mps[0] + y_c
        self.mps2 = mps[0]
        print("mps2:      " + str(mps[0]))

        mp2 = mp1 * mps[0]
        self.mp2 = mp2
        print("mp2:       " + str(mp2))

        #----------entrotus

        cpst1 = cps(Ttus, betaus)
        self.cpst1 = cpst1
        print("cpst1:     " + str(cpst1))

        EF0t1 = EF(T0, Ttus, betaus)
        self.EF0t1 = EF0t1
        print("EF0t1:     " + str(EF0t1))

        Ds0t1 = Ds(EF0t1, p0, ptus)
        self.Ds0t1 = Ds0t1
        print("Ds0t1:     " + str(Ds0t1))

        Dh0t1 = Dh(T0, Ttus, betaus)
        self.Dh0t1 = Dh0t1
        print("Dh0tm:     " + str(Dh0t1))

        #----------entrotm

        cpstm = cps(Ttm, betam)
        self.cpstm = cpstm
        print("cpstm:     " + str(cpstm))

        EF0tm = EF(T0, Ttm, betam)
        self.EF0tm = EF0tm
        print("EF0tm:     " + str(EF0tm))

        Ds0tm = Ds(EF0tm, p0, ptm)
        self.Ds0tm = Ds0tm
        print("Ds0tm:     " + str(Ds0tm))

        Dh0tm = Dh(T0, Ttm, (betam))
        self.Dh0tm = Dh0tm
        print("Dh0tm:     " + str(Dh0tm))

        #----------Tus, cpdtusus

        if useAM == "A":

            Mus = Ma_It(mp1, Ttus, ptus, Aaxus, betaus)
            self.Mus = Mus
            print("Mus:       " + str(Mus))
        
        Tus = T_func_It(Ttus, Mus, betaus)
        self.Tus = Tus
        print("Tus:       " + str(Tus))

        cpdtusus = cpd(Ttus, Tus, 0)
        self.cpdtusus = cpdtusus
        print("cpdtusus:  " + str(cpdtusus))

        #----------status

        pus = p(ptus, Tus, Ttus, 0)
        self.pus = pus
        print("pus:       " + str(pus))

        rhous = rho(pus, Tus)
        self.rhous = rhous
        print("rhous:     " + str(rhous))

        Vpus = Vp(mp1, rhous)
        self.Vpus = Vpus
        print("Vpus:      " + str(Vpus))

        #----------velus

        aus = a(R, Tus, betaus)
        self.aus = aus
        print("aus:       " + str(aus))

        cus = c(Mus, aus)
        self.cus = cus
        print("cus:       " + str(cus))
        
        caxus = cax(cus, alphaus)
        self.caxus = caxus
        print("caxus:     " + str(caxus))

        if useAM == "M":

            Aaxus = Aax(Vpus, caxus)
            self.Aaxus = Aaxus
            print("Aaxus:     " + str(Aaxus))

        #----------entrous

        cpsus = cps(Tus, betaus)
        self.cpsus = cpsus
        print("cpsus:     " + str(cpsus))

        EF0us = EF(T0, Tus, betaus)
        self.EF0us = EF0us
        print("EF0us:     " + str(EF0us))

        Ds0us = Ds(EF0us, p0, pus)
        self.Ds0us = Ds0us
        print("Ds0us:     " + str(Ds0us))

        Dh0us = Dh(T0, Tus, betaus)
        self.Dh02 = Dh0us
        print("Dh0us:     " + str(Dh0us))

        #----------Tm, cpdtmm

        if useAM == "A":

            Mm = Ma_It(mp2, Ttm, ptm, Aaxm, betam)
            self.Mm = Mm
            print("Mm:        " + str(Mm))
        
        Tm = T_func_It(Ttm, Mm, betam)
        self.Tm = Tm
        print("Tm:        " + str(Tm))

        cpdtmm = cpd(Ttm, Tm, betam)
        self.cpdtmm = cpdtmm
        print("cpdtmm:    " + str(cpdtmm))

        #----------statm

        pm = p(ptm, Tm, Ttm, betam)
        self.pm = pm
        print("pm:        " + str(pm))

        rhom = rho(pm, Tm)
        self.rhom = rhom
        print("rhom:      " + str(rhom))

        Vpm = Vp(mp2, rhom)
        self.Vpm = Vpm
        print("Vpm:       " + str(Vpm))

        #----------velm

        am = a(R, Tm, betam)
        self.am = am
        print("am:        " + str(am))

        cm = c(Mm, am)
        self.cm = cm
        print("cm:        " + str(cm))
        
        caxm = cax(cm, alpham)
        self.caxm = caxm
        print("caxm:      " + str(caxm))

        if useAM == "M":

            Aaxm = Aax(Vpm, caxm)
            self.Aaxm = Aaxm
            print("Aaxm:      " + str(Aaxm))

        #----------entro0m

        cpsm = cps(Tm, betam)
        self.cpsm = cpsm
        print("cpsm:      " + str(cpsm))

        EF0m = EF(T0, Tm, betam)
        self.EF0m = EF0m
        print("EF0m:      " + str(EF0m))

        Ds0m = Ds(EF0m, p0, pm)
        self.Ds0m = Ds0m
        print("Ds0m:      " + str(Ds0m))

        Dh0m = Dh(T0, Tm, betam)
        self.Dh0m = Dh0m
        print("Dh0m:      " + str(Dh0m))
        
        print("----------------------------------------")

        self.List1 = []

        self.List2 = []

        self.List1.extend([self.mps1, self.mp1, Ttus, "ct", self.Tus, "ct", ptus, "ct", "ct", self.cpst1, self.EF0t1, self.Dh0t1, betaus, self.Mus, self.caxus, self.pus, self.rhous, self.Aaxus])

        self.List2.extend([self.mps2, self.mp2, self.Ttm, "ct", self.Tm, Pm, self.ptm, "ct", "ct", self.cpstm, self.EF0tm, self.Dh0tm, betam, self.Mm, self.caxm, self.pm, self.rhom, self.Aaxm])

#Stat5 = mixing_class([1 - 0.18 + 0.022169], 1, 1700, 2654399.86, 1, 800, 0.027035, 0.024632, 0.08)