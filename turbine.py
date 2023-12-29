from side import *
from cp import *
from modules import *

class turbine_class:

    def __init__(self, mps, mp, Tt1, pt1, beta, neps, ix, useAM, Aax1, Aax2 ,M1, M2, alpha1, alpha2, atsum, nm, Dm1, Dm2, umax, name):

        print(name)

        print()

        print("mps: ", mps[0], "  mps:", mp, "  Tt1:", Tt1, "  pt1:", pt1, "  beta:", beta, "  neps:", neps, "  ix:", ix)
        print("useAM:", useAM, "  A1:", Aax1, "  A2:", Aax2, "  M1:", M1, "  M2:", M2, "  alpha1:", alpha1, "  alpha2:", alpha2, "  atsum:", atsum, "  nm:", nm, "  Dm1:", Dm1, "  Dm2:", Dm2, "  umax:", umax)

        print()

        self.Tt1 = Tt1
        self.pt1 = pt1
        self.Aax1 = Aax1
        self.Aax2 = Aax2
        self.M1 = M1
        self.M2 = M2

        mps1 = mps[0]
        self.mps1 = mps1
        print("mps1:      " + str(mps[0]))

        mp1 = mp * mps[0]
        self.mp1 = mp1
        print("mp1:       " + str(mp1))

        #----------downstream

        mps2 = mps[0]
        self.mps2 = mps2
        print("mps2:      " + str(mps[0]))

        mp2 = mp * mps[0]
        self.mp2 = mp2
        print("mp2:       " + str(mp2))

        #----------at, Ppe, Tt2s
        
        att1t2 = at(atsum, mps[0], nm)
        self. att1t2 = att1t2
        print("att1t2:    " + str(att1t2))

        Tt2 = Tt_exp_func_It(Tt1, att1t2, beta)
        self.Tt2 = Tt2
        print("Tt2:       " + str(Tt2))

        cpdt1t2 = cpd(Tt1, Tt2, beta)
        self.cpdt1t2 = cpdt1t2
        print("cpdt1t2:   " + str(cpdt1t2))

        Dht1t2 = Dh(Tt1, Tt2, beta)
        self.Dht1t2 = Dht1t2
        print("Dht1t2:    " + str(Dht1t2))

        Ppe = Ppe_func(att1t2, Tt1, Tt2, Tt2, neps, ix, beta)
        self.Ppe = Ppe
        print("Ppe:       " + str(Ppe))

        Tt2s = Tiscomp(Tt1, (1 / Ppe), R, cpd(Tt1, Tt2, beta))
        self.Tt2s = Tt2s
        print("Tt2s:      " + str(Tt2s))

        cpdt1t2s = cpd(Tt1, Tt2s, beta)
        self.cpdt1t2s = cpdt1t2s
        print("cpdt1t2s:  " + str(cpdt1t2s))

        Ppes = Ppe_func_It(att1t2, Tt1, Tt2, neps, ix, beta, (1 / Ppe))
        self.Ppes = Ppes
        print("Ppes:      " + str(Ppes))

        pt2 = pt1 * (1 / Ppes)
        self.pt2 = pt2
        print("pt2:       " + str(pt2))

        #----------entrot1

        cpst1 = cps(Tt1, beta)
        self.cpst1 = cpst1
        print("cpst1:     " + str(cpst1))

        EF0t1 = EF(T0, Tt1, beta)
        self.EF0t1 = EF0t1
        print("EF0t1:     " + str(EF0t1))

        Ds0t1 = Ds(EF0t1, p0, pt1)
        self.Ds0t1 = Ds0t1
        print("Ds0t1:     " + str(Ds0t1))

        Dh0t1 = Dh(T0, Tt1, beta)
        self.Dh0t1 = Dh0t1
        print("Dh0t1:     " + str(Dh0t1))

        #----------entrot2

        cpst2 = cps(Tt2, beta)
        self.cpst2 = cpst2
        print("cpst2:     " + str(cpst2))

        EF0t2 = EF(T0, Tt2, beta)
        self.EF0t2 = EF0t2
        print("EF0t2:     " + str(EF0t2))

        Ds0t2 = Ds(EF0t2, p0, pt2)
        self.Ds0t2 = Ds0t2
        print("Ds0t2:     " + str(Ds0t2))

        Dh0t2 = Dh(T0, Tt2, beta)
        self.Dh0t2 = Dh0t2
        print("Dh0t2:     " + str(Dh0t2))

        #----------T1, cpdt11

        if useAM == "A":

            M1 = Ma_It(mp1, Tt1, pt1, Aax1, beta)
            self.M1 = M1
            print("M1:        " + str(M1))

        T1 = T_func_It(Tt1, M1, beta)
        self.T1 = T1
        print("T1:        " + str(T1))

        cpdt11 = cpd(Tt1, T1, beta)
        self.cpdt11 = cpdt11
        print("cpdt11:    " + str(cpdt11))

        #----------stat1

        p1 = p(pt1, T1, Tt1, beta)
        self.p1 = p1
        print("ps1:       " + str(p1))

        rho1 = rho(p1, T1)
        self.rho1 = rho1
        print("rho1:      " + str(rho1))

        Vp1 = Vp(mp1, rho1)
        self.Vp1 = Vp1
        print("Vp1:       " + str(Vp1))

        #----------vel1

        a1 = a(R, T1, beta)
        self.a1 = a1
        print("a1:        " + str(a1))

        c1 = c(M1, a1)
        self.c1 = c1
        print("c1:        " + str(c1))
        
        cax1 = cax(c1, alpha1)
        self.cax1 = cax1
        print("cax1:      " + str(cax1))

        if useAM == "M":

            Aax1 = Aax(Vp1, cax1)
            self.Aax1 = Aax1
            print("Aax1:      " + str(Aax1))

        #----------entro1

        cps1 = cps(T1, beta)
        self.cps1 = cps1
        print("cps1:      " + str(cps1))

        EF01 = EF(T0, T1, beta)
        self.EF01 = EF01
        print("EF01:      " + str(EF01))

        Ds01 = Ds(EF01, p0, p1)
        self.Ds01 = Ds01
        print("Ds01:      " + str(Ds01))

        Dh01 = Dh(T0, T1, beta)
        self.Dh01 = Dh01
        print("Dh01:      " + str(Dh01))

        #----------T2, cpdt22, mp2

        if useAM == "A":

            M2 = Ma_It(mp2, Tt2, pt2, Aax2, beta)
            self.M2 = M2
            print("M2:        " + str(M2))
        
        T2 = T_func_It(Tt2, M2, beta)
        self.T2 = T2
        print("T2:        " + str(T2))

        cpdt22 = cpd(Tt2, T2, beta)
        self.cpdt22 = cpdt22
        print("cpdt22:    " + str(cpdt22))

        mp2 = mp1
        self.mp2 = mp2
        print("mp2:       " + str(mp2))

        #----------stat2

        p2 = p(pt2, T2, Tt2, beta)
        self.p2 = p2
        print("ps2:       " + str(p2))

        rho2 = rho(p2, T2)
        self.rho2 = rho2
        print("rho2:      " + str(rho2))

        Vp2 = Vp(mp2, rho2)
        self.Vp2 = Vp2
        print("Vp2:       " + str(Vp2))

        #----------vel2

        a2 = a(R, T2, beta)
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

        #----------entro2

        cps2 = cps(T2, beta)
        self.cps2 = cps2
        print("cps2:      " + str(cps2))

        EF02 = EF(T0, T2, beta)
        self.EF02 = EF02
        print("EF02:      " + str(EF02))

        Ds02 = Ds(EF02, p0, p2)
        self.Ds02 = Ds02
        print("Ds02:      " + str(Ds02))

        Dh02 = Dh(T0, T2, beta)
        self.Dh02 = Dh02
        print("Dh02:      " + str(Dh02))

        #----------powermetrics

        att1t2 = Dh(Tt1, Tt2, beta)
        self.att1t2 = att1t2
        print("att1t2:    " + str(att1t2))

        nest1t2 = nes(cpdt1t2, cpdt1t2s, Tt1, Ppes, R, att1t2)
        self.nest1t2 = nest1t2
        print("nest1t2:   " + str(nest1t2))

        nept1t2 = nep(cpdt1t2, cpdt1t2s, Tt1, Ppes, R, att1t2)
        self.nept1t2 = nept1t2
        print("nept1t2:   " + str(nept1t2))

        EFt1t2 = EF(Tt1, Tt2, beta)
        self.EFt1t2 = EFt1t2
        print("EFt1t2:    " + str(EFt1t2))

        Dst1t2 = Ds(EFt1t2, pt1, pt2)
        self.Dst1t2 = Dst1t2
        print("Dst1t2:    " + str(Dst1t2))

        Det1t2s = De(cpdt1t2s, Tt1, Tt2s)
        self.Det1t2s = Det1t2s
        print("Det1t2s:   " + str(Det1t2s))

        Det1t2sDht1t2 = Dht1t2 / Det1t2s
        self.Det1t2s = Det1t2sDht1t2
        print("De/Dh:     " + str(Det1t2sDht1t2))

        Dbt1t2 = Db(att1t2, Det1t2s)
        self.Dbt1t2 = Dbt1t2
        print("Dbt1t2:    " + str(Dbt1t2))

        Pt1t2 = P(mp1, att1t2) / 1000
        self.Pt1t2 = Pt1t2
        print("Pt1t2:     " + str(Pt1t2))

        #----------geometry

        Da1 = Da(Dm1, Aax1)
        self.Da1 = Da1
        print("Da1:       " + str(Da1))

        Da2 = Da(Dm2, Aax2)
        self.Da2 = Da2
        print("Da2:       " + str(Da2))

        Di1 = Di(Dm1, Aax1)
        self.Di1 = Di1
        print("Di1:       " + str(Di1))

        Di2 = Di(Dm2, Aax2)
        self.Di2 = Di2
        print("Di2:       " + str(Di2))

        Nmax1 = Nmax(umax, Da1)
        self.Nmax1 = Nmax1
        print("Nmax1:     " + str(Nmax1))

        Nmax2 = Nmax(umax, Da2)
        self.Nmax2 = Nmax2
        print("Nmax2:     " + str(Nmax2))

        print("----------------------------------------")

        self.List1 = []

        self.List2 = []

        self.List1.extend([self.mps1, self.mp1, Tt1, Tt1, self.T1, "ct", pt1, "ct", "ct", self.cpst1, self.EF0t1, self.Dh0t1, beta, self.M1, self.c1, self.p1, self.rho1, self.Aax1, Dm1, self.Di1, self.Da1, self.Nmax1])

        self.List2.extend([self.mps2, self.mp2, self.Tt2, self.Tt2s, self.T2, self.Ppes, self.pt2, self.att1t2, self.Pt1t2, self.cpst2, self.EF0t2, self.Dh0t2, beta, self.M2, self.c2, self.p2, self.rho2, self.Aax2, Dm2, self.Di2, self.Da2, self.Nmax2])
 
#test = turbine_class([1], 1, 1600, 2600000, 0.0233, 0.907, 1, 0.5, 0.5, 0, 0, 300000, 1)