import numpy as np

from MATRICES.ETAXJA import ETAXJA

def FORCE(PHI, POS, NCONEC, SURFACE, KCONEC, N, NE, Region):

    # CALCULATION OF FORCE COMPONENTS USING GAUSS QUADRATURE
    # FUNCTION FOR IMPERMEABLE STRUCTURES

    # Adjust for local or global domain
    NE_prev = sum(NE[0:Region])
    N_prev  = sum(N[0:Region])

    # For a total of 36 points in a 6x6 grid
    NPG = 6
    NG  = NPG*NPG

    # Gauss weights
    OME = np.array([0.467913934572691,0.467913934572691,
                    0.360761573048139,0.360761573048139,
                    0.171324492379170,0.171324492379170])

    # Element nodal coordinates
    XJ1 = np.zeros(NCONEC) 
    XJ2 = np.zeros(NCONEC) 
    XJ3 = np.zeros(NCONEC)

    # Nodal values of PHI
    PHI_ELEM = np.zeros(NCONEC, dtype=complex)

    # Values of PHI at Gauss points
    PHI_G   = np.zeros(NG,dtype=complex)

    # Declare and reset forces components
    FX = 0+0j 
    FY = 0+0j 
    FZ = 0+0j

    for IEJ in SURFACE + NE_prev:

        # Reset Gauss arrays per element
        PHI_G[:] = 0
        
        # Element nodes
        for M in range(NCONEC):
            KK = KCONEC[M,IEJ]
            XJ1[M] = POS[KK - N_prev,0]
            XJ2[M] = POS[KK - N_prev,1]
            XJ3[M] = POS[KK - N_prev,2]
            PHI_ELEM[M] = PHI[KK]

        # Gauss geometry for regular integration
        ETA1, ETA2, ETA3, XJA, XG1, XG2, XG3, FF = \
            ETAXJA(XJ1, XJ2, XJ3, NCONEC)       
        
        for i in range(NPG):
            for j in range(NPG):

                # Order number of Gauss point
                K = NPG*i + j

                # Interpolate potential
                for ik in range(NCONEC):
                    PHI_G[K] += FF[ik,K]*PHI_ELEM[ik]

                # Compute normal × Jacobian
                nJ1 = ETA1[K]*XJA[K]
                nJ2 = ETA2[K]*XJA[K]
                nJ3 = ETA3[K]*XJA[K]

                # Weight
                weight = OME[i]*OME[j]

                # Force components without rho*omega
                FX += PHI_G[K]*nJ1*weight
                FY += PHI_G[K]*nJ2*weight
                FZ += PHI_G[K]*nJ3*weight                

    return FX, FY, FZ