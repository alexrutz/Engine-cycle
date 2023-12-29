from inlet import *
from compressor import *
from compressor_detail import *
from burner import *
from mixing import *
from turbine import *
from turbine_detail import *
from cooling import *
from nozzle import *
from cp import *
from modules import *
from side import *
from presets import preset_list

from astwgui_new import Ui_ASTW

from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc

class ASTW(qtw.QMainWindow, Ui_ASTW):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setupUi(self)

        self.magic.clicked.connect(self.magie)
        self.stageBox_compressor_1.currentIndexChanged.connect(self.readStage_1_compressor)
        self.stageBox_compressor_2.currentIndexChanged.connect(self.readStage_2_compressor)
        self.stageBox_turbine_1.currentIndexChanged.connect(self.readStage_1_turbine)
        self.stageBox_turbine_2.currentIndexChanged.connect(self.readStage_2_turbine)
        self.presets.currentIndexChanged.connect(self.loadPreset)

    def magie(self):

        print("****************************************************************************************************start")

        #----------readValues

        #Tamb = float(self.Tamb.text())
        #pamb = float(self.pamb.text())
        #Tref = float(self.Tref.text())
        #pref = float(self.pref.text())
        R = float(self.R.text())
        Hu = float(self.Hu.text())
        nm = float(self.nm.text())
        M0 = float(self.M0.text())
        H0 = float(self.H0.text())
        P1 = float(self.P1.text())
        P2 = float(self.P2.text())
        Pb = float(self.Pb.text())
        Pin = float(self.Pin.text())
        Pon = float(self.Pon.text())
        PI = float(self.PI.text())
        PmV1 = float(self.PmV1.text())
        PmHPT = float(self.PmHPT.text())
        PmLPT = float(self.PmLPT.text())
        F = float(self.F.text())
        mp2 = float(self.mp2.text())
        bypass = float(self.bypass.text())
        PF = float(self.PF.text())
        y_401 = float(self.y_401.text())
        y_41 = float(self.y_41.text())
        y_5 = float(self.y_5.text())
        nF = float(self.nF.text())
        nLPC = float(self.nLPC.text())
        nHPC = float(self.nHPC.text())
        nHPT = float(self.nHPT.text())
        nLPT = float(self.nLPT.text())
        nB = float(self.nB.text())
        Tt4 = float(self.Tt4.text())
        Cr = float(self.Cr.text())
        OTDF = float(self.OTDF.text())
        uHP = float(self.uHP.text())
        uLP = float(self.uLP.text())
        Z_comp = int(self.Z_comp.text())
        cax21cax2_comp = float(self.cax21cax2_comp.text())
        Rct_comp = float(self.Rct_comp.text())
        arR_comp = float(self.arR_comp.text())
        arS_comp = float(self.arS_comp.text())
        tsR_comp = float(self.tsR_comp.text())
        tsS_comp = float(self.tsS_comp.text())
        Z_turb = int(self.Z_turb.text())
        cax21cax2_turb = float(self.cax21cax2_turb.text())
        Rct_turb = float(self.Rct_turb.text())
        arR_turb = float(self.arR_turb.text())
        arS_turb = float(self.arS_turb.text())
        tsR_turb = float(self.tsR_turb.text())
        tsS_turb = float(self.tsS_turb.text())
        Tm_S = float(self.Tm_S.text())
        Tm_R = float(self.Tm_R.text())
        nS = float(self.nS.text())
        nR = float(self.nR.text())
        mpsleak = float(self.mpsleak.text())

        nIt_c = int(self.nIt_c.text())
        nIt_t = int(self.nIt_t.text())
        nSt_c = int(self.nSt_c.text())
        nSt_t = int(self.nSt_t.text())
        if self.c_d_active.isChecked() == True:
            status_c_d = True
        if self.c_d_inactive.isChecked() == True:
            status_c_d = False
        if self.t_d_active.isChecked() == True:
            status_t_d = True
        if self.t_d_inactive.isChecked() == True:
            status_t_d = False
        if self.adapt.isChecked() == True:
            status = 0
        if self.notadapt.isChecked() == True:
            status = 1
        if self.Ma_button.isChecked() == True:
            useAM = "M"
        if self.A_button.isChecked() == True:
            useAM = "A"
        if self.mpCore.isChecked() == True:
            pass
        if self.mpEng.isChecked() == True:
            mp2 = mp2 / (1 + bypass)
        
        Ma1 = readGeo(self.geo, "Ma1")
        Ma2 = readGeo(self.geo, "Ma2")
        alpha1 = readGeo(self.geo, "alpha1")
        alpha2 = readGeo(self.geo, "alpha2")
        Dm1 = readGeo(self.geo, "Dm1")
        Dm2 = readGeo(self.geo, "Dm2")
        A1 = readGeo(self.geo, "A1")
        A2 = readGeo(self.geo, "A2")


        #coolingLabels = [self.pt22, self.pt5, self.pt22pt5, self.pt3, self.pt41, self.pt3pt41, self.pt3_2, self.pt401, self.pt3pt401]

        T0 = Th(H0)
        p0 = ph(T0)
        rho0 = rh(p0, R, T0)

        amb = ["ct", "ct", "ct", "ct", T0, "ct", p0, "ct", "ct", "ct", "ct", "ct", "ct", "ct", "ct", "ct", p0]

        #for value in range(0, 2):

        mps = [1]
        Inlet = inlet_class(mps, mp2, T0, p0, M0, PI, P1, P2, Cr, bypass, useAM, A1.Inlet, A2.Inlet, Ma1.Inlet, Ma2.Inlet, alpha1.Inlet, alpha2.Inlet, "Inlet")
        Fan = compressor_class([bypass], mp2, Inlet.Tt2, Inlet.Tt2, Inlet.pt2, PF, 0, 0, 0, 0, nF, 1, useAM, A1.Fan, A2.Fan, Ma1.Fan, Ma2.Fan, alpha1.Fan, alpha2.Fan, Dm1.Fan, Dm2.Fan, uLP, 0, "Fan")       
        LPC = compressor_class(mps, mp2, Inlet.Tt2, Inlet.Tt2, Inlet.pt2, Inlet.PLPC, 0, 0, 0, 0, nLPC, 1, useAM, A1.LPC, A2.LPC, Ma1.LPC, Ma2.LPC, alpha1.LPC, alpha2.LPC, Dm1.LPC, Dm2.LPC, uLP, 0, "LPC")
        
        Tt3s = Tts_func_It(LPC.Tt2, Inlet.PHPC, 0, R)
        print("Tt3s:      " + str(Tt3s))

        cpdt21t3s = cpd(LPC.Tt2, Tt3s, 0)
        print("cpdt21t3s: " + str(cpdt21t3s))

        Tt3 = Tiscompdps(LPC.Tt2, Inlet.PHPC, R, cpdt21t3s, nHPC, 1)
        print("Tt3:       " + str(Tt3))

        cpdt21t3 = cpd(LPC.Tt2, Tt3, 0)
        print("cpdt21t3:  " + str(cpdt21t3))

        print("----------------------------------------")

        HPC1 = compressor_class(mps, mp2, LPC.Tt2, LPC.Tt2s, LPC.pt2, Inlet.PHPC1, 0, cpdt21t3s, cpdt21t3, 1, nHPC, 1, useAM, A1.HPC1, A2.HPC2, Ma1.HPC1, Ma2.HPC1, alpha1.HPC1, alpha2.HPC1, Dm1.HPC1, Dm2.HPC1, uHP, y_5, "HPC1")
        HPC2 = compressor_class(mps, mp2, HPC1.Tt2, HPC1.Tt2s, HPC1.pt2, Inlet.PHPC2, 0, cpdt21t3s, cpdt21t3, 1, nHPC, 1, useAM, A1.HPC2, A2.HPC2, Ma1.HPC2, Ma2.HPC2, alpha1.HPC2, alpha2.HPC2, Dm1.HPC2, Dm2.HPC2, uHP, y_41 + y_401, "HPC2")

        burner = burner_class(mps, mp2, HPC2.Tt2, HPC2.pt2, Pb, nB, Hu, Tt4, y_401, y_41, y_5, OTDF, useAM, A1.Burner, A2.Burner, Ma1.Burner, Ma2.Burner, alpha1.Burner, alpha2.Burner, Dm1.Burner, Dm2.Burner, "Burner")
        mixing_HPTV1 = mixing_class(mps, mp2, Tt4, burner.pt4, PmV1, HPC2.Tt2, burner.beta4, burner.beta_401m, y_401, useAM, A1.mixHPTV1, A2.mixHPTV1, Ma1.mixHPTV1, Ma2.mixHPTV1, alpha1.mixHPTV1, alpha2.mixHPTV1, "mixing_HPTV1")

        HPT = turbine_class(mps, mp2, mixing_HPTV1.Ttm, mixing_HPTV1.ptm, burner.beta_401m, nHPT, 1, useAM, A1.HPT, A2.HPT, Ma1.HPT, Ma2.HPT, alpha1.HPT, alpha2.HPT, HPC1.avt1t2 + HPC2.avt1t2 * (1 - y_5), nm, Dm1.HPT, Dm2.HPT, uHP, "HPT")
        mixing_HPT = mixing_class(mps, mp2, HPT.Tt2, HPT.pt2, PmHPT, HPC2.Tt2, burner.beta_401m, burner.beta_41m, y_41, useAM, A1.mixHPT, A2.mixHPT, Ma1.mixHPT, Ma2.mixHPT, alpha1.mixHPT, alpha2.mixHPT, "mixing_HPT")

        LPT = turbine_class(mps, mp2, mixing_HPT.Ttm, mixing_HPT.ptm, burner.beta_41m, nLPT, 1, useAM, A1.LPT, A2.LPT, Ma1.LPT, Ma2.LPT, alpha1.LPT, alpha2.LPT, LPC.avt1t2 + Fan.avt1t2 * bypass, nm, Dm1.LPT, Dm2.LPT, uLP, "LPT")
        mixing_LPT = mixing_class(mps, mp2, LPT.Tt2, LPT.pt2, PmLPT, LPC.Tt2, burner.beta_41m, burner.beta_5m, y_5, useAM, A1.mixLPT, A2.mixLPT, Ma1.mixLPT, Ma2.mixLPT, alpha1.mixLPT, alpha2.mixLPT, "mixing_LPT")

        nozzle = nozzle_class(mps, mp2, mixing_LPT.Ttm, mixing_LPT.ptm, Fan.Tt2, Fan.pt2, burner.betaB, Pin, Pon, M0, T0, p0, bypass, Hu, nm, status, F, Fan.avt1t2, mixing_LPT, Fan, "Nozzle")
        #mp2 = nozzle.mp2

        perf_components = [amb, Inlet.List, Fan.List1, Fan.List2, LPC.List1, LPC.List2, HPC1.List1, HPC1.List2, HPC2.List1, HPC2.List2, burner.List1, burner.List2, mixing_HPTV1.List1, mixing_HPTV1.List2, HPT.List1, HPT.List2, mixing_HPT.List1, mixing_HPT.List2, LPT.List1, LPT.List2, mixing_LPT.List1, mixing_LPT.List2, nozzle.ListC2, nozzle.ListB1, nozzle.ListB2]

        NHP = min(HPC1.Nmax1, HPC1.Nmax2, HPC2.Nmax1, HPC2.Nmax2, HPT.Nmax1, HPT.Nmax2)
        NLP = min(Fan.Nmax1, Fan.Nmax2, LPC.Nmax1, LPC.Nmax2, LPT.Nmax1, LPT.Nmax2)

        setTable(perf_components, self.perf_table, 1, 9)

        #----------process

        setCooling(self, HPC1.pt2, HPC2.pt2, burner.pt4, HPT.pt2, LPT.pt2)

        setEff(self, nozzle)

        Tt_diag = [Inlet.Tt1, Inlet.Tt2, LPC.Tt2, HPC1.Tt2, HPC2.Tt2, Tt4, mixing_HPTV1.Ttm, HPT.Tt2, mixing_HPT.Ttm, LPT.Tt2, mixing_LPT.Ttm, nozzle.Tt9]
        s_diag = [Inlet.Ds0t1, Inlet.Ds0t2, LPC.Ds0t2, HPC1.Ds0t2, HPC2.Ds0t2, burner.Ds0t4, mixing_HPTV1.Ds0tm, HPT.Ds0t2, mixing_HPT.Ds0tm, LPT.Ds0t2, mixing_LPT.Ds0tm, nozzle.Ds0t9]
        beta_diag = [0, 0, 0, 0, 0, burner.beta4, mixing_HPTV1.betam, mixing_HPTV1.betam, mixing_HPT.betam, mixing_HPT.betam, mixing_LPT.betam, burner.betaB] 
        #beta_diag = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        p_diag = [Inlet.pt0, Inlet.pt2 ,LPC.pt2, HPC1.pt2, HPC2.pt2, burner.pt4, mixing_HPTV1.ptm, HPT.pt2, mixing_HPT.ptm, LPT.pt2, mixing_LPT.ptm, nozzle.pt9]
        v_diag = [1/rho0, 1/Inlet.rho2, 1/LPC.rho2, 1/HPC1.rho2, 1/HPC2.rho2, 1/burner.rho4, 1/mixing_HPTV1.rhom, 1/HPT.rho2, 1/mixing_HPT.rhom, 1/LPT.rho2, 1/mixing_LPT.rhom, 1/nozzle.rho9]

        anno_main = ["0", "2", "21", "22", "3", "401", "401m", "41", "41m", "5", "5m", "9"]
        Ts_offset_main = [-50, -50, -50, -50, -50, +30, +30, +30, -70, +30, -70, -70]
        pv_offset_main = [-0.05, 0.05, 0.05, 0.05, -0.05, 0.05, -0.05, 0.05, -0.05, 0.05, -0.05, 0.05]

        setDiags_Tspv(self, Tt_diag, s_diag, beta_diag, p_diag, v_diag, anno_main, Ts_offset_main, pv_offset_main, self.Ts_C, self.pv_C, 68, 68)

        #----------compressor detail

        createStagesLists_c(self, Z_comp)

        if status_c_d == True:

            stages_c = stages_inputlist(LPC.Tt2, HPC2.Tt2, Ma1.HPC1, Ma2.HPC2, alpha1.HPC1, alpha2.HPC2, Dm1.HPC1, Dm2.HPC2, Z_comp, cpdt21t3s, nHPC)
            stages_c_ar = stages_area(HPC1.Da1, HPC2.Da2, HPC1.Di1, HPC2.Di2, Z_comp)

            for i, value in enumerate(self.st_c):

                if i >= nSt_c:
                    continue

                if i == 0:

                    for g in range(0, nIt_c):

                        self.st_c[i] = compressor_class([1], mp2, LPC.Tt2, LPC.Tt2s, LPC.pt2, stages_c[i + 1][0], 0, cpdt21t3s, cpdt21t3, 1, nHPC, 1, useAM, stages_c_ar[i][0], stages_c_ar[i][1], stages_c[i + 1][1][0], stages_c[i + 1][1][1], stages_c[i + 1][2][0], stages_c[i + 1][2][1], stages_c[i + 1][3][0], stages_c[i + 1][3][1], uHP, 0, "compressor stage: " + str(i + 1) + "  Iteration: " + str(g + 1))
                        self.st_c_d[i] = compressor_detail(mp2, self.st_c[i].Tt1, self.st_c[i].pt1, self.st_c[i].Pic, stages_c[0][0], NHP, self.st_c[i].c1, cax21cax2_comp, self.st_c[i].avt1t2, Rct_comp, nHPC, stages_c[i + 1][2][0], stages_c[i + 1][2][1], stages_c[i + 1][3][0], stages_c[i + 1][3][1], self.st_c[i].Di1, self.st_c[i].Di2, self.st_c[i].Da1, self.st_c[i].Da2, arR_comp, arS_comp, i, tsR_comp, tsS_comp, "compressor stage detail: " + str(i + 1) + "  Iteration: " + str(g + 1))
                        #self.st_c_d[i] = compressor_detail(mp2, self.st_c[i].Tt1, self.st_c[i].pt1, self.st_c[i].Pic, stages_c[0][0], NHP, self.st_c[i].c1, cax21cax2_comp, self.st_c[i].avt1t2, Rct_comp, nHPC, 20 * (pi / 180), 20 * (pi / 180), stages_c[i + 1][3][0], stages_c[i + 1][3][1], self.st_c[i].Di1, self.st_c[i].Di2, self.st_c[i].Da1, self.st_c[i].Da2, arR_comp, arS_comp, i, tsR_comp, tsS_comp, "compressor stage detail: " + str(i + 1) + "  Iteration: " + str(g + 1))
                        stages_c[i + 1][1][1] = self.st_c_d[i].MC3
                        stages_c[i + 2][1][0] = self.st_c_d[i].MC3

                else:

                    for g in range(0, nIt_c):

                        self.st_c[i] = compressor_class([1], mp2, st_Tt1, st_Tt1s, st_pt1, stages_c[i + 1][0], 0, cpdt21t3s, cpdt21t3, 1, nHPC, 1, useAM, stages_c_ar[i][0], stages_c_ar[i][1], stages_c[i + 1][1][0], stages_c[i + 1][1][1], stages_c[i + 1][2][0], stages_c[i + 1][2][1], stages_c[i + 1][3][0], stages_c[i + 1][3][1], uHP, 0, "compressor stage: " + str(i + 1) + "  Iteration: " + str(g + 1))
                        self.st_c_d[i] = compressor_detail(mp2, self.st_c[i].Tt1, self.st_c[i].pt1, self.st_c[i].Pic, stages_c[0][0], NHP, self.st_c[i].c1, cax21cax2_comp, self.st_c[i].avt1t2, Rct_comp, nHPC, stages_c[i + 1][2][0], stages_c[i + 1][2][1], stages_c[i + 1][3][0], stages_c[i + 1][3][1], self.st_c[i].Di1, self.st_c[i].Di2, self.st_c[i].Da1, self.st_c[i].Da2, arR_comp, arS_comp, i, tsS_comp, tsS_comp, "compressor stage detail: " + str(i + 1) + "  Iteration: " + str(g + 1))
                        #self.st_c_d[i] = compressor_detail(mp2, self.st_c[i].Tt1, self.st_c[i].pt1, self.st_c[i].Pic, stages_c[0][0], NHP, self.st_c[i].c1, cax21cax2_comp, self.st_c[i].avt1t2, Rct_comp, nHPC, 20 * (pi / 180), 20 * (pi / 180), stages_c[i + 1][3][0], stages_c[i + 1][3][1], self.st_c[i].Di1, self.st_c[i].Di2, self.st_c[i].Da1, self.st_c[i].Da2, arR_comp, arS_comp, i, tsR_comp, tsS_comp, "compressor stage detail: " + str(i + 1) + "  Iteration: " + str(g + 1))
                        stages_c[i + 1][1][1] = self.st_c_d[i].MC3
                        stages_c[i + 2][1][0] = self.st_c_d[i].MC3
                
                st_Tt1 = self.st_c[i].Tt2
                st_Tt1s = self.st_c[i].Tt2s
                st_pt1 = self.st_c[i].pt2
                self.st_c_d_input.append(self.st_c_d[i].st_inputlist)
                self.st_c_d_output.append(self.st_c_d[i].st_outputlist)
                self.st_c_d_vel.append(self.st_c_d[i].st_vel)
                self.st_c_d_veltri2.append(self.st_c_d[i].st_veltri2)
                self.st_c_d_veltri21.append(self.st_c_d[i].st_veltri21)
                self.st_c_d_canal.append(self.st_c_d[i].st_canal)
                self.st_c_d_strain.append(self.st_c_d[i].st_strain)
                self.st_c_list.append(self.st_c[i].st_List)
                self.st_c_d_T.append(self.st_c_d[i].st_T)
                self.st_c_d_s.append(self.st_c_d[i].st_s)
                self.st_c_d_p.append(self.st_c_d[i].st_p)
                self.st_c_d_v.append(self.st_c_d[i].st_v)
                self.c_T.append(self.st_c[i].Tt2)
                self.c_s.append(self.st_c[i].Ds0t2)
                self.c_p.append(self.st_c[i].pt2)
                self.c_v.append(1 / self.st_c[i].rho2)
                self.st_c_d_beta.append([0, 0, 0])
                self.c_beta.append(0)
                self.stageBox_compressor_1.addItem(str(i + 1))
                self.stageBox_compressor_2.addItem(str(i + 1))

            setTable(self.st_c_list, self.stages_compressor_table, 0, 9)

            setDiags_canal_full(self, self.st_c_d_canal, self.compressor_canal_all)

            anno_compressor = range(1, Z_comp + 1)
            Ts_offset_compressor = [-10 for n in range(0, Z_comp)]
            pv_offset_compressor = [0.01 for n in range(0, Z_comp)]

            setDiags_Tspv(self, self.c_T, self.c_s, self.c_beta, self.c_p, self.c_v, anno_compressor, Ts_offset_compressor, pv_offset_compressor, self.compressor_Ts, self.compressor_pv, 50, 50)

        #----------turbine detail

        createStagesLists_t(self, Z_turb)

        if status_t_d == True:

            stagepressure_turbine = HPT.Ppes ** (1 / Z_turb)

            stages_t = stages_inputlist(100, 1000, Ma1.HPT, Ma2.HPT, alpha1.HPT, alpha2.HPT, Dm1.HPT, Dm2.mixHPT, Z_turb, 1000, 1)
            stages_t_ar = stages_area(HPT.Da1, HPT.Da2, HPT.Di1, HPT.Di2, Z_turb)

            for i, value in enumerate(self.st_t):

                if i >= nSt_t:
                    continue

                if i == 0:

                    for g in range(0, nIt_t):

                        self.st_t_mix[i] = mixing_class([burner.mps2], mp2, Tt4, burner.pt4, PmV1, HPC2.Tt2, burner.beta4, burner.beta_401m, y_401, useAM, stages_t_ar[i][0], stages_t_ar[i][1], Ma1.mixHPTV1, Ma2.mixHPTV1, alpha1.mixHPTV1, alpha2.mixHPTV1, "turbine stage mixing: " + str(i + 1) + "  Iteration: " + str(g + 1))
                        self.st_t[i] = compressor_class([mixing_HPTV1.mps2], mixing_HPTV1.mp2, mixing_HPTV1.Ttm, 1000, mixing_HPTV1.ptm, 1 / stagepressure_turbine, burner.beta_401m, HPT.cpdt1t2s, HPT.cpdt1t2, 1, 1 / nHPT, 1, useAM, stages_t_ar[i][0], stages_t_ar[i][1], stages_t[i + 1][1][0], stages_t[i + 1][1][1], stages_t[i + 1][2][0], stages_t[i + 1][2][1], stages_t[i + 1][3][0], stages_t[i + 1][3][1], NHP, 0, "turbine stage: " + str(i + 1) + "  Iteration: " + str(g + 1))
                        self.st_t_d[i] = turbine_detail(self.st_t[i].mp1, self.st_t_mix[i].Ttus, self.st_t_mix[i].Ttm, self.st_t[i].pt1, self.st_t[i].Pic, self.st_t[i].Tt2 - self.st_t[i].Tt1, NHP, self.st_t_mix[i].cus, cax21cax2_turb, self.st_t[i].avt1t2, Rct_turb, nHPT, stages_t[i + 1][2][0], stages_t[i + 1][2][1], stages_t[i + 1][3][0], stages_t[i + 1][3][1], self.st_t[i].Di1, self.st_t[i].Di2, self.st_t[i].Da1, self.st_t[i].Da2, arR_turb, arS_turb, Z_turb, tsR_turb, tsS_turb, burner.beta4, burner.beta_401m, "turbine stage detail: " + str(i + 1) + "  Iteration: " + str(g + 1))
                        self.st_t_cool[i] = cooling_class(mp2, burner.Tt4max, mixing_HPTV1.Ttm, mixing_HPTV1.Tm, Tm_S, Tm_R, HPC2.Tt2, burner.beta_401m, nS, nR, self.st_t_d[i].rho21, self.st_t_d[i].rho3, self.st_t_d[i].c21, self.st_t_d[i].w21, self.st_t_d[i].w3, self.st_t_d[i].sS, self.st_t_d[i].sR, self.st_t_d[i].ZS, self.st_t_d[i].ZR, self.st_t_d[i].Di2, self.st_t_d[i].Di3, self.st_t_d[i].Da2, self.st_t_d[i].Da3, self.st_t_d[i].wS, self.st_t_d[i].wR, mpsleak, "cooling stage: " + str(i + 1) + "  Iteration: " + str(g + 1))

                else:
                    pass

                st_Tt1 = self.st_t[i].Tt2
                st_Tt1s = self.st_t[i].Tt2s
                st_pt1 = self.st_t[i].pt2
                self.st_t_d_input.append(self.st_t_d[i].st_inputlist)
                self.st_t_d_output.append(self.st_t_d[i].st_outputlist)
                self.st_t_d_vel.append(self.st_t_d[i].st_vel)
                self.st_t_d_veltri2.append(self.st_t_d[i].st_veltri2)
                self.st_t_d_veltri21.append(self.st_t_d[i].st_veltri21)
                self.st_t_d_canal.append(self.st_t_d[i].st_canal)
                self.st_t_d_strain.append(self.st_t_d[i].st_strain)
                self.st_t_list.append(self.st_t[i].st_List)
                self.st_t_d_T.append(self.st_t_d[i].st_T)
                self.st_t_d_s.append(self.st_t_d[i].st_s)
                self.st_t_d_p.append(self.st_t_d[i].st_p)
                self.st_t_d_v.append(self.st_t_d[i].st_v)
                self.t_T.append(self.st_t[i].Tt2)
                self.t_s.append(self.st_t[i].Ds0t2)
                self.t_p.append(self.st_t[i].pt2)
                self.t_v.append(1 / self.st_t[i].rho2)
                self.st_t_d_beta.append([0, 0, 0])
                self.t_beta.append(0)
                self.st_t_cool_output.append(self.st_t_cool[i].List)
                self.stageBox_turbine_1.addItem(str(i + 1))
                self.stageBox_turbine_2.addItem(str(i + 1))

            setTable(self.st_t_list, self.stages_turbine_table, 0, 9)

            setDiags_canal_full(self, self.st_t_d_canal, self.turbine_canal_all)

            anno_compressor = range(1, Z_comp + 1)
            Ts_offset_compressor = [-10 for n in range(0, Z_comp)]
            pv_offset_compressor = [0.01 for n in range(0, Z_comp)]

            setDiags_Tspv(self, self.t_T, self.t_s, self.t_beta, self.t_p, self.t_v, anno_compressor, Ts_offset_compressor, pv_offset_compressor, self.turbine_Ts, self.turbine_pv, 50, 50)

        print("****************************************************************************************************end")

    def readStage_1_compressor(self):

        z1 = int(self.stageBox_compressor_1.currentText()) - 1

        setTable_inv(self.st_c_d_vel[z1], self.stages_compressor_vel_1, 6)
        setTable_inv(self.st_c_d_input[z1], self.stages_compressor_input_1, 6)
        setTable_inv(self.st_c_d_output[z1], self.stages_compressor_output_1, 6)
        setDiags_veltri(self, self.st_c_d_veltri2[z1], self.st_c_d_veltri21[z1], self.stages_compressor_tri_1)
        setDiags_canal(self, self.st_c_d_canal[z1], self.stages_compressor_canal)
        setTable_inv(self.st_c_d_strain[z1], self.stages_compressor_strain, 6)
        #setDiags_geo(self, self.st_c_d[z1].beta2, self.st_c_d[z1].beta21, self.st_c_d[z1].alpha21, self.st_c_d[z1].alpha3, self.stages_compressor_geo)
        anno_stage = [2, 21, 3]
        Ts_offset_stage = [-10 for n in range(0, 3)]
        pv_offset_stage = [0.01 for n in range(0, 3)]
        setDiags_Tspv(self, self.st_c_d_T[z1], self.st_c_d_s[z1], self.st_c_d_beta[z1], self.st_c_d_p[z1], self.st_c_d_v[z1], anno_stage, Ts_offset_stage, pv_offset_stage, self.stages_compressor_Ts, self.stages_compressor_pv, 50, 50)
        self.stages_compressor_tri_1_anno.setText("Geschwindigkeitsdreieck Stufe " + str(z1 + 1))
        self.stages_compressor_canal_anno.setText("Kanalgeometrie Stufe " + str(z1 + 1))
        self.stages_compressor_geo_anno.setText("Schaufelgeometrie Stufe " + str(z1 + 1))
        self.stages_compressor_Ts_anno.setText("T-s-Diagramm Stufe " + str(z1 + 1))
        self.stages_compressor_pv_anno.setText("p-v-Diagramm Stufe  " + str(z1 + 1))

    def readStage_1_turbine(self):

        z1 = int(self.stageBox_turbine_1.currentText()) - 1

        setTable_inv(self.st_t_d_vel[z1], self.stages_turbine_vel_1, 6)
        setTable_inv(self.st_t_d_input[z1], self.stages_turbine_input_1, 6)
        setTable_inv(self.st_t_d_output[z1], self.stages_turbine_output_1, 6)
        setDiags_veltri(self, self.st_t_d_veltri2[z1], self.st_t_d_veltri21[z1], self.stages_turbine_tri_1)
        setDiags_canal(self, self.st_t_d_canal[z1], self.stages_turbine_canal)
        setTable_inv(self.st_t_d_strain[z1], self.stages_turbine_strain, 6)
        setTable_inv(self.st_t_cool_output[z1], self.stages_turbine_cooling, 6)
        #setDiags_geo(self, self.st_c_d[z1].beta2, self.st_c_d[z1].beta21, self.st_c_d[z1].alpha21, self.st_c_d[z1].alpha3, self.stages_compressor_geo)
        anno_stage = [2, 21, 3]
        Ts_offset_stage = [-10 for n in range(0, 3)]
        pv_offset_stage = [0.01 for n in range(0, 3)]
        setDiags_Tspv(self, self.st_t_d_T[z1], self.st_t_d_s[z1], self.st_t_d_beta[z1], self.st_t_d_p[z1], self.st_t_d_v[z1], anno_stage, Ts_offset_stage, pv_offset_stage, self.stages_turbine_Ts, self.stages_turbine_pv, 50, 50)
        self.stages_turbine_tri_1_anno.setText("Geschwindigkeitsdreieck Stufe " + str(z1 + 1))
        self.stages_turbine_canal_anno.setText("Kanalgeometrie Stufe " + str(z1 + 1))
        self.stages_turbine_geo_anno.setText("Schaufelgeometrie Stufe " + str(z1 + 1))
        self.stages_turbine_Ts_anno.setText("T-s-Diagramm Stufe " + str(z1 + 1))
        self.stages_turbine_pv_anno.setText("p-v-Diagramm Stufe  " + str(z1 + 1))

    def readStage_2_compressor(self):

        z2 = int(self.stageBox_compressor_2.currentText()) - 1

        setTable_inv(self.st_c_d_vel[z2], self.stages_compressor_vel_2, 6)
        setTable_inv(self.st_c_d_input[z2], self.stages_compressor_input_2, 6)
        setTable_inv(self.st_c_d_output[z2], self.stages_compressor_output_2, 6)
        setDiags_veltri(self, self.st_c_d_veltri2[z2], self.st_c_d_veltri21[z2], self.stages_compressor_tri_2)
        self.stages_compressor_tri_2_anno.setText("Geschwindigkeitsdreieck Stufe " + str(z2 + 1))

    def readStage_2_turbine(self):

        z2 = int(self.stageBox_turbine_2.currentText()) - 1

        setTable_inv(self.st_t_d_vel[z2], self.stages_turbine_vel_2, 6)
        setTable_inv(self.st_t_d_input[z2], self.stages_turbine_input_2, 6)
        setTable_inv(self.st_t_d_output[z2], self.stages_turbine_output_2, 6)
        setDiags_veltri(self, self.st_t_d_veltri2[z2], self.st_t_d_veltri21[z2], self.stages_turbine_tri_2)
        self.stages_turbine_tri_2_anno.setText("Geschwindigkeitsdreieck Stufe " + str(z2 + 1))

    def loadPreset(self):

        item = preset_list[self.presets.currentIndex()]
        self.R.setText(str(item.R))
        self.Hu.setText(str(item.Hu))
        self.nm.setText(str(item.nm))
        self.M0.setText(str(item.M0))
        self.H0.setText(str(item.H0))
        self.P1.setText(str(item.P1))
        self.P2.setText(str(item.P2))
        self.Pb.setText(str(item.Pb))
        self.Pin.setText(str(item.Pin))
        self.Pon.setText(str(item.Pon))
        self.PI.setText(str(item.PI))
        self.PmV1.setText(str(item.PmV1))
        self.PmHPT.setText(str(item.PmHPT))
        self.PmLPT.setText(str(item.PmLPT))
        self.F.setText(str(item.F))
        self.mp2.setText(str(item.mp2))
        self.bypass.setText(str(item.bypass))
        self.PF.setText(str(item.PF))
        self.y_401.setText(str(item.y_401))
        self.y_41.setText(str(item.y_41))
        self.y_5.setText(str(item.y_5))
        self.nF.setText(str(item.nF))
        self.nLPC.setText(str(item.nLPC))
        self.nHPC.setText(str(item.nHPC))
        self.nHPT.setText(str(item.nHPT))
        self.nLPT.setText(str(item.nLPT))
        self.nB.setText(str(item.nB))
        self.Tt4.setText(str(item.Tt4))
        self.Cr.setText(str(item.Cr))
        self.OTDF.setText(str(item.OTDF))
        self.uHP.setText(str(item.uHP))
        self.uLP.setText(str(item.uLP))
        self.Z_comp.setText(str(item.Z_comp))
        self.cax21cax2_comp.setText(str(item.cax21cax2_comp))
        self.Rct_comp.setText(str(item.Rct_comp))
        self.arR_comp.setText(str(item.arR_comp))
        self.arS_comp.setText(str(item.arS_comp))
        self.tsR_comp.setText(str(item.tsR_comp))
        self.tsS_comp.setText(str(item.tsS_comp))
        self.Z_turb.setText(str(item.Z_turb))
        self.cax21cax2_turb.setText(str(item.cax21cax2_turb))
        self.Rct_turb.setText(str(item.Rct_turb))
        self.arR_turb.setText(str(item.arR_turb))
        self.arS_turb.setText(str(item.arS_turb))
        self.tsR_turb.setText(str(item.tsR_turb))
        self.tsS_turb.setText(str(item.tsS_turb))
        self.Tm_S.setText(str(item.Tm_S))
        self.Tm_R.setText(str(item.Tm_R))
        self.nS.setText(str(item.nS))
        self.nR.setText(str(item.nR))
        self.mpsleak.setText(str(item.mpsleak))
        setTable(item.geo_list, self.geo, 0, 7)

if __name__ == '__main__':

    app = qtw.QApplication([])

    widget = ASTW()
    widget.show()

    app.exec_()