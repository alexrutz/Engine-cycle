aaa = [0.992313, 0.236688, -1.852148, 6.083152, -8.893933, 7.097112, -3.234725, 0.794571, -0.081873]
bbb = [-0.718874, 8.747481, -15.863157, 17.254096, -10.233795, 3.081778, -0.361112, -0.003919]

def cps(T, beta):

    cpsa = 0
    cpsb = 0

    for i, ainds in enumerate(aaa):
        cpsaadd = ainds * (T / 1000) ** i
        cpsa = cpsa + cpsaadd

    for j, binds in enumerate(bbb):
        cpsbadd = binds * (T / 1000) ** j
        cpsb = cpsb + cpsbadd

    cps = cpsa + (beta / (1 + beta)) * cpsb

    return 1000 * cps

def cpd(T1, T2, beta):

    if T1 == T2:

        return cps(T1, beta)
    
    cpdau = 0
    cpdbu = 0
    cpdal = 0
    cpdbl = 0

    for idu, ainddu in enumerate(aaa):
        
        cpdaaddu = (ainddu / (idu + 1)) * (T2 / 1000) ** (idu + 1)
        cpdau = cpdau + cpdaaddu
    
    for jdu, binddu in enumerate(bbb):
        
        cpdbaddu = (binddu / (jdu + 1)) * (T2 / 1000) ** (jdu + 1)
        cpdbu = cpdbu + cpdbaddu
    
    for idl, ainddl in enumerate(aaa):
        
        cpdaaddl = (ainddl / (idl + 1)) * (T1 / 1000) ** (idl + 1)
        cpdal = cpdal + cpdaaddl

    for jdl, binddl in enumerate(bbb):
        
        cpdbaddl = (binddl / (jdl + 1)) * (T1 / 1000) ** (jdl + 1)
        cpdbl = cpdbl + cpdbaddl
        
    cpd = (((cpdau + (beta / (1 + beta)) * cpdbu) - (cpdal + (beta / (1 + beta)) * cpdbl)) * 1000) / (T2 - T1)

    return 1000 * cpd