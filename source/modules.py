from PyQt5.QtCore import Qt
from side import *
from cp import *
from math import sqrt, atan, atan2, asin, acos
from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem
from scipy.optimize import fsolve

from PyQt5 import QtWidgets as qtw

import sys
import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

def Tts_func_It(Tt1, Pic, beta, R):

    for i in range(0, 6):
        if i == 0:
            Tt2s = Tiscomp(Tt1, Pic, R, cps(Tt1, beta))
        else:
            Tt2s = Tiscomp(Tt1, Pic, R, cpd12s)
        cpd12s = cpd(Tt1, Tt2s, beta)

    return Tt2s

def T_func_It(Tt, Ma, beta):

    for i in range(0, 6):
        if i == 0:
            T_IT = T(Tt, k(cps(Tt, beta)), Ma)
        else:
            T_IT = T(Tt, k(cpdlul), Ma)
        cpdlul = cpd(Tt, T_IT, beta)

    return T_IT

def Tt_func_It(T, Ma, beta):

    for i in range(0, 6):
        if i == 0:
            Tt_IT = Tt(T, k(cps(T, 0)), Ma)
        else:
            Tt_IT = Tt(T, k(cpdlul), Ma)
        cpdlul = cpd(T, Tt_IT, beta)

    return Tt_IT

def beta_func(T0, Tt3, Tt4, y_kl, Hu, nBK):

    for i in range(0, 6):
        if i == 0:
            betaBK = betaBK_func(T0, Tt3, Tt4, cps(Tt3, 0), 0.03, y_kl, Hu, nBK)
        else:
            betaBK = betaBK_func(T0, Tt3, Tt4, cpd0t42, beta42, y_kl, Hu, nBK)
        beta42 = beta4_func(betaBK, y_kl)
        cpd0t42 = cpd(T0, Tt4, beta42)

    return betaBK

def T4crit_func_It(Tt4, beta4):

    for i in range(0, 6):
        if i == 0:
            T4crit1 = T4crit(Tt4, cps(Tt4, beta4))
        else:
            T4crit1 = T4crit(Tt4, cpd4crit)
        cpd4crit = cpd(Tt4, T4crit1, beta4)

    return T4crit1
        
def Ttm_func_It(T0, Ttc, Ttms, y_c, mpsus, betaus, betam):

    for i in range(0, 6):
        if i == 0:
            Ttm = Ttm_func(T0, Ttc, Ttms, cpd(T0, Ttms, betam), y_c, mpsus, betaus)
        else:
            Ttm = Ttm_func(T0, Ttc, Ttms, cpdm, y_c, mpsus, betaus)
        cpdm = cpd(T0, Ttm, betam)

    return Ttm
        
def Tt_exp_func_It(Tt1, at, beta):

    for i in range(0, 6):
        if i == 0:
            Tt2 = Tt_exp_func(at, Tt1, cps(Tt1, beta))
        else:
            Tt2 = Tt_exp_func(at, Tt1, cpdt)
        cpdt = cpd(Tt1, Tt2, beta)

    return Tt2
        
def Ppe_func_It(at, Tt1, Tt2, neps, ix, beta, Ppe):

    for i in range(0, 6):
        if i == 0:
            Tts = Tiscomp(Tt1, Ppe, R, cpd(Tt1, Tt2, beta))
        else:
            Tts = Tiscomp(Tt1, (1 / Ppe), R, cpd(Tt1, Tts, beta))
        Ppe = Ppe_func(at, Tt1, Tt2, Tts, neps, ix, beta)

    return Ppe

def T_contour_It(Tt, c, beta):

    for i in range(0, 6):
        if i == 0:
            T = T_contour(Tt, c, cps(Tt, beta))
        else:
            T = T_contour(Tt, c, cpdc)
        cpdc = cpd(Tt, T, beta)

    return T

def Ma_It(mp, Tt, pt, A, beta):

    def rm(Ma):        
        return sqrt(k(cps(Tt, beta)) / R) * Ma * (1 + ((k(cps(Tt, beta)) - 1) / 2) * (Ma **2)) ** (- (k(cps(Tt, beta)) + 1) / (2 * (k(cps(Tt, beta)) - 1))) - ((mp * sqrt(Tt)) / (A * pt))

    Ma = fsolve(rm, 0.1)

    return Ma[0]

def EF(T1, T2, beta):

    Dsaadd = 0
    Dsbadd = 0

    for i, value in enumerate(aaa):
        if i != 0:
            Dsaadd = Dsaadd + ((value / i) * (((T2 / 1000) ** i) - ((T1 / 1000) ** i)))

    for i, value in enumerate(bbb):
        if i != 0:
            Dsbadd = Dsbadd + ((value / i) * (((T2 / 1000) ** i) - ((T1 / 1000) ** i)))

    return 1000 * ((aaa[0] * log(T2 / T1) + Dsaadd + (beta / (1 + beta)) * (bbb[0] * log(T2 / T1) + Dsbadd)) / R)

def fuelflows_func(betaBK, y_kl, y_401, y_41, y_5):

    beta_4 = betaBK /(1 - y_kl)
    beta_401m = betaBK /(1 - y_kl + y_401)
    beta_41m = betaBK /(1 - y_kl + y_401 + y_41)
    beta_5m = betaBK /(1 - y_kl + y_401 + y_41 + y_5)

    return [beta_4, beta_401m, beta_41m, beta_5m]

def massflows_func(mp2, y_kl, betaBK):

    mp3 = mp2 * (1 - y_kl)
    mp4 = mp2 * (1 - y_kl + betaBK)

    return [mp3, mp4]

class filler():
    pass

def readGeo(table, type):

    Geo = filler()

    if type == "Ma1":
        i = 0

    if type == "Ma2":
        i = 1

    if type == "A1":
        i = 2

    if type == "A2":
        i = 3

    if type == "alpha1":
        i = 4

    if type == "alpha2":
        i = 5

    if type == "Dm1":
        i = 6

    if type == "Dm2":
        i = 7

    Geo.Inlet = float(table.item(i, 0).text())
    Geo.Fan = float(table.item(i, 1).text())
    Geo.LPC = float(table.item(i, 2).text())
    Geo.HPC1 = float(table.item(i, 3).text())
    Geo.HPC2 = float(table.item(i, 4).text())
    Geo.Burner = float(table.item(i, 5).text())
    Geo.mixHPTV1 = float(table.item(i, 6).text())
    Geo.HPT = float(table.item(i, 7).text())
    Geo.mixHPT = float(table.item(i, 8).text())
    Geo.LPT = float(table.item(i, 9).text())
    Geo.mixLPT = float(table.item(i, 10).text())
    Geo.nozzle_C = float(table.item(i, 11).text())
    Geo.nozzle_B = float(table.item(i, 12).text())

    return Geo

def setTable(components, table, x, l):

    for m, component in enumerate(components):

        for n, value in enumerate(component):

            if value == "ct":
                continue

            cell = QTableWidgetItem(str(component[n])[:l])
            cell.setTextAlignment(0)
            table.setItem(m, n + x, cell)
            cell.setTextAlignment(Qt.AlignCenter)

def setTable_inv(components, table,l ):

    for m, component in enumerate(components):

        for n, value in enumerate(component):

            if value == "ct":
                continue

            cell = QTableWidgetItem(str(component[n])[:l])
            cell.setTextAlignment(0)
            table.setItem(n, m, cell)
            cell.setTextAlignment(Qt.AlignCenter)

def setCooling(form, pt22, pt3, pt401, pt41, pt5):

    form.pt22.setText(str(pt22)[:10])
    form.pt5.setText(str(pt5)[:10])
    form.pt22pt5.setText(str(pt22 / pt5)[:10])
    form.pt3.setText(str(pt3)[:10])
    form.pt41.setText(str(pt41)[:10])
    form.pt3pt41.setText(str(pt3 / pt41)[:10])
    form.pt3_2.setText(str(pt3)[:10])
    form.pt401.setText(str(pt401)[:10])
    form.pt3pt401.setText(str(pt3 / pt401)[:10])

def setEff(form, nozzle):

    form.kT0.setText(str(k(cps(T0, 0)))[:8])
    form.c0.setText(str(nozzle.c0)[:8])
    form.c9.setText(str(nozzle.c9)[:8])
    form.c19.setText(str(nozzle.c19)[:8])
    form.aFlow.setText(str(nozzle.aFlow)[:8])
    form.fsI.setText(str(nozzle.fsI)[:8])
    form.fsII.setText(str(nozzle.fsII)[:8])
    form.fsEng.setText(str(nozzle.fsEng)[:8])
    form.Feng.setText(str(nozzle.Feng)[:8])
    form.bsEng.setText(str(nozzle.bsEng)[:8])
    form.nth.setText(str(nozzle.nth)[:8])
    form.nP.setText(str(nozzle.nP)[:8])
    form.nCore.setText(str(nozzle.nCore)[:8])

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=10, height=10, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

def Tt_isobaric_It(T1, s1, s, beta):

    for i in range(0, 3):
        if i == 0:
            #Tt = Ttisobaric(T1, s1, s, cps(T1, beta))
            Tt = Ttisobaric(T1, s1, s, 1200)
        else:
            #Tt = Ttisobaric(T1, s1, s, cps(Tt, beta))
            Tt = Ttisobaric(T1, s1, s, 1200)

    return Tt

def isobaric_table(s1, T1, beta, smax):

    s = range(-100, smax, 10)
    T = []
    pairs = [s, T]

    for x in s:
        T.append(Tt_isobaric_It(T1, s1, x, beta))

    return pairs

def setDiags_Tspv(form, Tt_diag, s_diag, beta_diag, p_diag, v_diag, anno, Ts_offset, pv_offset, Ts_canvas, pv_canvas, dpiTs, dpipv):

        smax = int(max(s_diag) + 100)
        Ts = MplCanvas(form, width=10, height=10, dpi=dpiTs)
        Ts.axes.plot(s_diag, Tt_diag)

        for i in range(0, len(Tt_diag)):
            data = isobaric_table(s_diag[i], Tt_diag[i], beta_diag[i], smax)
            Ts.axes.plot(data[0], data[1], "k", linewidth=0.5)
            Ts.axes.annotate(anno[i], [s_diag[i] + Ts_offset[i], Tt_diag[i]])

        pv = MplCanvas(form, width=10, height=10, dpi=dpipv)
        pv.axes.plot(v_diag, p_diag)

        for i in range(0, len(p_diag)):
            pv.axes.annotate(anno[i], [v_diag[i] + pv_offset[i], p_diag[i]])

        Ts_scene = qtw.QGraphicsScene()
        Ts_scene.addWidget(Ts)
        Ts_canvas.setScene(Ts_scene)

        pv_scene = qtw.QGraphicsScene()
        pv_scene.addWidget(pv)
        pv_canvas.setScene(pv_scene)

def createStagesLists_t(form, Z_turb):

    form.st_t = list(range(0, Z_turb))
    form.st_t_list = []
    form.st_t_d = list(range(0, Z_turb))
    form.st_t_d_vel = []
    form.st_t_mix = list(range(0, Z_turb))
    form.st_t_cool = list(range(0, Z_turb))
    form.st_t_cool_output = []
    form.st_t_d_input = []
    form.st_t_d_output = []
    form.st_t_d_veltri2 = []
    form.st_t_d_veltri21 = []
    form.st_t_d_canal = []
    form.st_t_d_strain = []
    form.st_t_d_beta = []
    form.st_t_d_T = []
    form.st_t_d_s = []
    form.st_t_d_p = []
    form.st_t_d_v = []
    form.t_beta = []
    form.t_T = []
    form.t_s = []
    form.t_p = []
    form.t_v = []

def createStagesLists_c(form, Z_comp):

    form.st_c = list(range(0, Z_comp))
    form.st_c_list = []
    form.st_c_d = list(range(0, Z_comp))
    form.st_c_d_vel = []
    form.st_c_d_input = []
    form.st_c_d_output = []
    form.st_c_d_veltri2 = []
    form.st_c_d_veltri21 = []
    form.st_c_d_canal = []
    form.st_c_d_strain = []
    form.st_c_d_beta = []
    form.st_c_d_T = []
    form.st_c_d_s = []
    form.st_c_d_p = []
    form.st_c_d_v = []
    form.c_beta = []
    form.c_T = []
    form.c_s = []
    form.c_p = []
    form.c_v = []

def setDiags_veltri(form, veltri1_2_list, veltri1_21_list, canvas):

    veltri = MplCanvas(form, width=4.5, height=4.5, dpi=88)
    veltri.axes.plot(veltri1_2_list[0], veltri1_2_list[1], label="2")
    veltri.axes.plot(veltri1_21_list[0], veltri1_21_list[1], label="21")
    veltri_scene = qtw.QGraphicsScene()
    veltri_scene.addWidget(veltri)
    canvas.setScene(veltri_scene)

def setDiags_canal(form, canal_list, canvas):

    canal = MplCanvas(form, width=4.5, height=4.5, dpi=88)
    canal.axes.plot(canal_list[0], canal_list[1])
    canal_scene = qtw.QGraphicsScene()
    canal_scene.addWidget(canal)
    canvas.setScene(canal_scene)

def setDiags_canal_full(form, canal_list_full, canvas):

    canal_full = MplCanvas(form, width=6, height=4.5, dpi=88)
    
    for stage in canal_list_full:

        canal_full.axes.plot(stage[0], stage[1], "tab:blue")

    canal_full.axes.plot(0, 0)
    canal_full_scene = qtw.QGraphicsScene()
    canal_full_scene.addWidget(canal_full)
    canvas.setScene(canal_full_scene)

def stages_area(Da1, Da2, Di1, Di2, Z):

    areas = []
    dDa = (Da2 - Da1) / Z
    dDi = (Di2 - Di1) / Z

    for i in range(0, Z):

        interstage_area = [0, 0]
        interstage_area[0] = (pi / 4) * (Da1 + (dDa * i) - (Di1 + (dDi * i))) ** 2
        interstage_area[1] = (pi / 4) * (Da1 + (dDa * (i + 1)) - (Di1 + (dDi * (i + 1)))) ** 2
        areas.append(interstage_area)

    return areas

def stages_inputlist(Tt1, Tt2, Ma1, Ma2, alpha1, alpha2, Dm1, Dm2, Z, cpdt1t2s, nvp):

    inputs = []
    dTt = (Tt2 - Tt1) / Z
    dMa = (Ma2 - Ma1) / Z
    dalpha = (alpha2 - alpha1) / Z
    dDm = (Dm2 - Dm1) / Z

    inputs.extend([[dTt, dMa, dalpha, dDm]])

    for value in range(0, Z):
        interstage_Ma = [0, 0]
        interstage_alpha = [0, 0]
        interstage_Dm = [0, 0]
        group = []

        if value == 0:
            Ttz1 = Tt1
            Maz = Ma1 + (dMa / 2)
            alphaz = alpha1 + (dalpha / 2)
            Dmz = Dm1 + (dDm / 2)

        Ttz2 = Ttz1 + dTt
        group.append(stagepressure(Ttz1, Ttz2, cpdt1t2s, nvp))
        Ttz1 = Ttz1 + dTt

        interstage_Ma[0] = Maz - (dMa / 2)
        interstage_Ma[1] = Maz + (dMa / 2)
        group.append(interstage_Ma)
        Maz = Maz + dMa

        interstage_alpha[0] = alphaz - (dalpha / 2)
        interstage_alpha[1] = alphaz + (dalpha / 2)
        group.append(interstage_alpha)
        alphaz = alphaz + dalpha

        interstage_Dm[0] = Dmz - (dDm / 2)
        interstage_Dm[1] = Dmz + (dDm / 2)
        group.append(interstage_Dm)
        Dmz = Dmz + dDm

        inputs.append(group)

    return inputs

def D21_It(Dm2, Dm3, Da2, Da3, Di2, Di3, arS, arR):

    D21 = []

    for i in range(0, 6):

        if i == 0:

            Di21 = D(Di2, Di3)
            Dm21 = D(Dm2, Dm3)
            Da21 = D(Da2, Da3)

        else:

            Di21 = D21_func(Di2, Di3, wR, wS)
            Dm21 = D21_func(Dm2, Dm3, wR, wS)
            Da21 = D21_func(Da2, Da3, wR, wS)

        wS = w_func(Da21, Da3, Di21, Di3, arS)
        wR = w_func(Da2, Da21, Di2, Di21, arR)
    
    D21.extend([Di21, Dm21, Da21, wR, wS, wR + wS])

    return D21

def setDiags_geo(form, beta1, beta2, alpha1, alpha2, canvas):

    geometry = MplCanvas(form, width=4.5, height=4.5, dpi=88)
    geometry.axes.plot([1, 1],[0, 6], "tab:blue")
    geometry.axes.plot([3, 3],[0, 6], "tab:blue")
    geometry.axes.plot([4, 4],[0, 6], "tab:blue")
    geometry.axes.plot([6, 6],[0, 6], "tab:blue")
    skellet_alpha = createSkellet([1, 3], beta1, beta2)
    skellet_beta = createSkellet([4, 3], alpha1, alpha2)
    geometry.axes.plot(skellet_alpha[0], skellet_alpha[1], "k")
    geometry.axes.plot(skellet_beta[0], skellet_beta[1], "k")
    geometry_scene = qtw.QGraphicsScene()
    geometry_scene.addWidget(geometry)
    canvas.setScene(geometry_scene)


def circle(start, end, root, radius):

    start = int(10000 * start)
    end = int(10000 * end)
    lul = [start, end]
    startlul = min(lul)
    endlul = max(lul)
    print("circle start: ", start, end, root ,radius)

    x_coords = []
    y_coords = []

    for degree in range(startlul, endlul + 1):

        rad = degree / 10000

        x = (cos(rad) * radius) + root[0]
        y = (sin(rad) * radius) + root[1]

        x_coords.append(x)
        y_coords.append(y)

    coords = [x_coords, y_coords]

    print("resolution of circle:   ", len(x_coords))

    print("**************************************************")

    return coords

def createSkellet(anchor, beta1, beta2):

    if beta1 == beta2:

        beta1 = beta1 + 0.5

    if beta1 < beta2:

        print("skellet starting values:", beta1 * (180 / pi), beta2 * (180 / pi), anchor)

        r = 0
        dX = sin(beta1) * r

        while r <= dX + 2.001:

            dX = sin(beta1) * r
            r = r + 0.001

        beta2_calc = asin((sin(beta1) * r + 2) / r)

        print("initial r found: ", r)
        print("beta2:           ", beta2 * (180 / pi), beta2)
        print("beta2_calc_start:", beta2_calc * (180 / pi))
        
        while beta2_calc > beta2:

            beta2_calc = asin((sin(beta1) * r + 2) / r)
            r = r + 0.001

        dY = cos(beta1) * r
        dX = sin(beta1) * r

        print("beta2_calc_final:", beta2_calc * (180 / pi))
        print("dY, dX:", dY, dX)
        print("r:     ", r)
        
        return circle((pi / 2) - beta1, (pi / 2) - beta2, [anchor[0] - dX, anchor[1] - dY], r)

    if beta1 > beta2:

        print("skellet starting values:", beta1 * (180 / pi), beta2 * (180 / pi), anchor)

        r = 0
        dX = - sin(beta1) * r

        while r <= dX + 2.001:

            dX = - sin(beta1) * r
            r = r + 0.001

        beta2_calc = acos(- (sin(beta1) * r + 2) / r) - (pi / 2)

        print("initial r found: ", r)
        print("beta2:           ", beta2 * (180 / pi), beta2)
        print("beta2_calc_start:", beta2_calc * (180 / pi), beta2_calc)
        
        while beta2_calc < beta2:

            beta2_calc = acos((- (sin(beta1) * r) + 2) / r) - (pi / 2)
            r = r + 0.001

        dY = - cos(beta1) * r
        dX = - sin(beta1) * r

        print("beta2_calc_final:", beta2_calc * (180 / pi))
        print("dY, dX:", dY, dX)
        print("r:     ", r)
        
        return circle((3 / 2) * pi - beta1, (3 / 2) * pi - beta2, [anchor[0] - dX, anchor[1] - dY], r)

#createSkellet([1, 3], 20, 70)    

#----------inconsistencies: Ma and M, ps and p, av -> ac, Pic and Pc, compressor and c, delete R from all arguments