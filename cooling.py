from side import *
from cp import *


class cooling_class:
    
    def __init__(self, mp2, Ttmax, Tthg_R, Thg_R, Tm_S, Tm_R, Tc, betam, nc_S, nc_R, rho_S, rho_R, c_S, w_R1, w_R2, s_S, s_R, z_S, z_R, Di_S, Di_R, Da_S, Da_R, b_S, b_R, mpsleak, name):

        print(name)

        print()

        print("mp2:", mp2, "  Ttmax:", Ttmax, "  Tthg_R:" , Tthg_R, "  Thg_R:", Thg_R, "  Tm_S:", Tm_S, "  Tm_R:", Tm_R, "  Tc:", Tc)
        print("betam:", betam, "  nc_S:", nc_S,  "  nc_R:", nc_R, "  rho_S:", rho_S, "  rho_R:", rho_R, "  c_S:", c_S, "  w_R1:", w_R1, "  w_R2:", w_R2)
        print("s_S:", s_S, "  s_R:", s_R, "  z_S:", z_S, "  z_R:", z_R, "  Di_S:", Di_S, "  Di_R:", Di_R, "  Da_S:", Da_S, "  Da_R:", Da_R, "  b_S:", b_S, "  b_R:", b_R, "  mpsleak:", mpsleak)

        print()

        self.mpsleak = mpsleak

        #----------stator 

        epsilon_S = epsilon_func(Ttmax, Tm_S, Tc)
        self.epsilon_S = epsilon_S
        print("epsilon_S: " + str(epsilon_S))

        mpf_S = mpf_func(epsilon_S, nc_S)
        self.mpf_S = mpf_S
        print("mpf_S:     " + str(mpf_S))

        Re_S = Re_func(rho_S, c_S, s_S)
        self.Re_S = Re_S
        print("Re_S:      " + str(Re_S))

        Nu_S = Nu_func(Re_S)
        self.Nu_S = Nu_S
        print("Nu_S:      " + str(Nu_S))

        alpha_hg_S = alpha_hg_func(Nu_S, s_S)
        self.alpha_hg_S = alpha_hg_S
        print("alpha_hg_S:" + str(alpha_hg_S))

        A_net_S = A_net_func(z_S, s_S, Da_S, Di_S, b_S)
        self.A_net_S = A_net_S
        print("A_net_S:   " + str(A_net_S))

        mpc_S = mpc_func(mpf_S, alpha_hg_S, A_net_S, Tc)
        self.mpc_S = mpc_S
        print("mpc_S:     " + str(mpc_S))

        #----------rotor

        Thg_Rw = Thg_R + (w_R1 ** 2) / (2 * cpd(Tthg_R, Thg_R, betam))
        self.Thg_Rw = Thg_Rw
        print("Thg_Rw:    " + str(Thg_Rw))

        epsilon_R = epsilon_func(Thg_Rw, Tm_R, Tc)
        self.epsilon_R = epsilon_R
        print("epsilon_R: " + str(epsilon_R))

        mpf_R = mpf_func(epsilon_R, nc_R)
        self.mpf_R = mpf_R
        print("mpf_R:     " + str(mpf_R))

        Re_R = Re_func(rho_R, w_R2, s_R)
        self.Re_R = Re_R
        print("Re_R:      " + str(Re_R))

        Nu_R = Nu_func(Re_R)
        self.Nu_R = Nu_R
        print("Nu_R:      " + str(Nu_R))

        alpha_hg_R = alpha_hg_func(Nu_R, s_R)
        self.alpha_hg_R = alpha_hg_R
        print("alpha_hg_R:" + str(alpha_hg_R))

        A_net_R = A_net_func(z_R, s_R, Da_R, Di_R, b_R)
        self.A_net_R = A_net_R
        print("A_net_R:   " + str(A_net_R))

        mpc_R = mpc_func(mpf_R, alpha_hg_R, A_net_R, Tc)
        self.mpc_R = mpc_R
        print("mpc_R:     " + str(mpc_R))

        #----------summary

        mpleak = mpsleak * mp2
        self.mpleak = mpleak
        print("mpleak:    " + str(mpleak))

        mpsc_S = mpc_S / mp2
        self.mpsc_S = mpsc_S
        print("mpsc_S:    " + str(mpsc_S))

        mpsc_R = mpc_R / mp2
        self.mpsc_R = mpsc_R
        print("mpsc_R:    " + str(mpsc_R))

        mpc = mpc_R + mpc_S + mpleak
        self.mpc = mpc
        print("mpc:       " + str(mpc))

        mpsc = mpsc_R + mpsc_S + mpsleak
        self.mpsc = mpsc
        print("mpsc:      " + str(mpsc))
        
        print("----------------------------------------")

        self.List = []

        self.List.extend([[Ttmax, epsilon_S, mpf_S, Re_S, Nu_S, alpha_hg_S, A_net_S, mpc_S, mpsc_S, mpc], [Thg_Rw, epsilon_R, mpf_R, Re_S, Nu_R, alpha_hg_R, A_net_R, mpc_R, mpsc_R, mpsc]])