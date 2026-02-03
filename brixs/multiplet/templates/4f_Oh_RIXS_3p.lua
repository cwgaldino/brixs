--------------------------------------------------------------------------------
-- Quanty input file generated using Crispy. If you use this file please cite
-- the following reference: http://dx.doi.org/10.5281/zenodo.1008184.
--
-- elements: 4f
-- symmetry: Oh
-- experiment: RIXS
-- edge: M2,3 (3p)
--------------------------------------------------------------------------------
Verbosity($Verbosity)

--------------------------------------------------------------------------------
-- Initialize the Hamiltonians.
--------------------------------------------------------------------------------
H_i = 0
H_m = 0
H_f = 0

--------------------------------------------------------------------------------
-- Toggle the Hamiltonian terms.
--------------------------------------------------------------------------------
H_atomic = $H_atomic
H_crystal_field = $H_crystal_field
H_4f_ligands_hybridization_lmct = $H_4f_ligands_hybridization_lmct
H_magnetic_field = $H_magnetic_field
H_exchange_field = $H_exchange_field

--------------------------------------------------------------------------------
-- Define the number of electrons, shells, etc.
--------------------------------------------------------------------------------
NBosons = 0
NFermions = 20

NElectrons_3p = 6
NElectrons_4f = $NElectrons_4f

IndexDn_3p = {0, 2, 4}
IndexUp_3p = {1, 3, 5}
IndexDn_4f = {6, 8, 10, 12, 14, 16, 18}
IndexUp_4f = {7, 9, 11, 13, 15, 17, 19}

if H_4f_ligands_hybridization_lmct == 1 then
    NFermions = 34

    NElectrons_L1 = 14

    IndexDn_L1 = {20, 22, 24, 26, 28, 30, 32}
    IndexUp_L1 = {21, 23, 25, 27, 29, 31, 33}
end

--------------------------------------------------------------------------------
-- Define the atomic term.
--------------------------------------------------------------------------------
N_3p = NewOperator('Number', NFermions, IndexUp_3p, IndexUp_3p, {1, 1, 1})
     + NewOperator('Number', NFermions, IndexDn_3p, IndexDn_3p, {1, 1, 1})

N_4f = NewOperator('Number', NFermions, IndexUp_4f, IndexUp_4f, {1, 1, 1, 1, 1, 1, 1})
     + NewOperator('Number', NFermions, IndexDn_4f, IndexDn_4f, {1, 1, 1, 1, 1, 1, 1})

if H_atomic == 1 then
    F0_4f_4f = NewOperator('U', NFermions, IndexUp_4f, IndexDn_4f, {1, 0, 0, 0})
    F2_4f_4f = NewOperator('U', NFermions, IndexUp_4f, IndexDn_4f, {0, 1, 0, 0})
    F4_4f_4f = NewOperator('U', NFermions, IndexUp_4f, IndexDn_4f, {0, 0, 1, 0})
    F6_4f_4f = NewOperator('U', NFermions, IndexUp_4f, IndexDn_4f, {0, 0, 0, 1})

    F0_3p_4f = NewOperator('U', NFermions, IndexUp_3p, IndexDn_3p, IndexUp_4f, IndexDn_4f, {1, 0}, {0, 0})
    F2_3p_4f = NewOperator('U', NFermions, IndexUp_3p, IndexDn_3p, IndexUp_4f, IndexDn_4f, {0, 1}, {0, 0})
    G2_3p_4f = NewOperator('U', NFermions, IndexUp_3p, IndexDn_3p, IndexUp_4f, IndexDn_4f, {0, 0}, {1, 0})
    G4_3p_4f = NewOperator('U', NFermions, IndexUp_3p, IndexDn_3p, IndexUp_4f, IndexDn_4f, {0, 0}, {0, 1})

    U_4f_4f_i = $U(4f,4f)_i_value
    F2_4f_4f_i = $F2(4f,4f)_i_value * $F2(4f,4f)_i_scale
    F4_4f_4f_i = $F4(4f,4f)_i_value * $F4(4f,4f)_i_scale
    F6_4f_4f_i = $F6(4f,4f)_i_value * $F6(4f,4f)_i_scale
    F0_4f_4f_i = U_4f_4f_i + 4 / 195 * F2_4f_4f_i + 2 / 143 * F4_4f_4f_i + 100 / 5577 * F6_4f_4f_i

    U_4f_4f_m = $U(4f,4f)_m_value
    F2_4f_4f_m = $F2(4f,4f)_m_value * $F2(4f,4f)_m_scale
    F4_4f_4f_m = $F4(4f,4f)_m_value * $F4(4f,4f)_m_scale
    F6_4f_4f_m = $F6(4f,4f)_m_value * $F6(4f,4f)_m_scale
    F0_4f_4f_m = U_4f_4f_m + 4 / 195 * F2_4f_4f_m + 2 / 143 * F4_4f_4f_m + 100 / 5577 * F6_4f_4f_m
    U_3p_4f_m = $U(3p,4f)_m_value
    F2_3p_4f_m = $F2(3p,4f)_m_value * $F2(3p,4f)_m_scale
    G2_3p_4f_m = $G2(3p,4f)_m_value * $G2(3p,4f)_m_scale
    G4_3p_4f_m = $G4(3p,4f)_m_value * $G4(3p,4f)_m_scale
    F0_3p_4f_m = U_3p_4f_m + 3 / 70 * G2_3p_4f_m + 2 / 63 * G4_3p_4f_m

    U_4f_4f_f = $U(4f,4f)_f_value
    F2_4f_4f_f = $F2(4f,4f)_f_value * $F2(4f,4f)_f_scale
    F4_4f_4f_f = $F4(4f,4f)_f_value * $F4(4f,4f)_f_scale
    F6_4f_4f_f = $F6(4f,4f)_f_value * $F6(4f,4f)_f_scale
    F0_4f_4f_f = U_4f_4f_f + 4 / 195 * F2_4f_4f_f + 2 / 143 * F4_4f_4f_f + 100 / 5577 * F6_4f_4f_f


    H_i = H_i + Chop(
          F0_4f_4f_i * F0_4f_4f
        + F2_4f_4f_i * F2_4f_4f
        + F4_4f_4f_i * F4_4f_4f
        + F6_4f_4f_i * F6_4f_4f)

    H_m = H_m + Chop(
          F0_4f_4f_m * F0_4f_4f
        + F2_4f_4f_m * F2_4f_4f
        + F4_4f_4f_m * F4_4f_4f
        + F6_4f_4f_m * F6_4f_4f
        + F0_3p_4f_m * F0_3p_4f
        + F2_3p_4f_m * F2_3p_4f
        + G2_3p_4f_m * G2_3p_4f
        + G4_3p_4f_m * G4_3p_4f)
        
    H_f = H_f + Chop(
          F0_4f_4f_f * F0_4f_4f
        + F2_4f_4f_f * F2_4f_4f
        + F4_4f_4f_f * F4_4f_4f
        + F6_4f_4f_f * F6_4f_4f)

    ldots_4f = NewOperator('ldots', NFermions, IndexUp_4f, IndexDn_4f)

    ldots_3p = NewOperator('ldots', NFermions, IndexUp_3p, IndexDn_3p)

    zeta_4f_i = $zeta(4f)_i_value * $zeta(4f)_i_scale

    zeta_4f_m = $zeta(4f)_m_value * $zeta(4f)_m_scale
    zeta_3p_m = $zeta(3p)_m_value * $zeta(3p)_m_scale

    zeta_4f_f = $zeta(4f)_f_value * $zeta(4f)_f_scale

    H_i = H_i + Chop(
          zeta_4f_i * ldots_4f)

    H_m = H_m + Chop(
          zeta_4f_m * ldots_4f
        + zeta_3p_m * ldots_3p)
    
    H_f = H_f + Chop(
          zeta_4f_f * ldots_4f)

end

--------------------------------------------------------------------------------
-- Define the crystal field term.
--------------------------------------------------------------------------------
if H_crystal_field == 1 then
    -- PotentialExpandedOnClm('Oh', 3, {Ea2u, Et1u, Et2u})
    -- Ea2u_4f = NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, PotentialExpandedOnClm('Oh', 3, {1, 0, 0}))
    -- Et2u_4f = NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, PotentialExpandedOnClm('Oh', 3, {0, 1, 0}))
    -- Et1u_4f = NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, PotentialExpandedOnClm('Oh', 3, {0, 0, 1}))

    B40_4f_i = $B40(4f)_i_value
    B60_4f_i = $B60(4f)_i_value

    Akm_4f_i = {
        {4,  0, B40_4f_i},
        {4, -4, math.sqrt(5/14) * B40_4f_i},
        {4,  4, math.sqrt(5/14) * B40_4f_i},
        {6,  0, B60_4f_i},
        {6, -4, -math.sqrt(7/2) * B60_4f_i},
        {6,  4, -math.sqrt(7/2) * B60_4f_i},
    }

    io.write('Energies of the 4f orbitals in the initial Hamiltonian (crystal field term only):\n')
    io.write('================\n')
    io.write('Irrep.        E\n')
    io.write('================\n')
    io.write(string.format('a2u     %8.3f\n', -4 / 11 * B40_4f_i +  80 / 143 * B60_4f_i))
    io.write(string.format('t1u     %8.3f\n',  2 / 11 * B40_4f_i + 100 / 429 * B60_4f_i))
    io.write(string.format('t2u     %8.3f\n', -2 / 33 * B40_4f_i -  60 / 143 * B60_4f_i))
    io.write('================\n')
    io.write('\n')

    B40_4f_m = $B40(4f)_m_value
    B60_4f_m = $B60(4f)_m_value

    Akm_4f_m = {
        {4,  0, B40_4f_m},
        {4, -4, math.sqrt(5/14) * B40_4f_m},
        {4,  4, math.sqrt(5/14) * B40_4f_m},
        {6,  0, B60_4f_m},
        {6, -4, -math.sqrt(7/2) * B60_4f_m},
        {6,  4, -math.sqrt(7/2) * B60_4f_m},
    }


    B40_4f_f = $B40(4f)_f_value
    B60_4f_f = $B60(4f)_f_value

    Akm_4f_f = {
        {4,  0, B40_4f_f},
        {4, -4, math.sqrt(5/14) * B40_4f_f},
        {4,  4, math.sqrt(5/14) * B40_4f_f},
        {6,  0, B60_4f_i},
        {6, -4, -math.sqrt(7/2) * B60_4f_f},
        {6,  4, -math.sqrt(7/2) * B60_4f_f},
    }

    H_i = H_i + Chop(NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, Akm_4f_i))

    H_m = H_m + Chop(NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, Akm_4f_m))

    H_f = H_f + Chop(NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, Akm_4f_f))

end

--------------------------------------------------------------------------------
-- Define the 4f-ligands hybridization term.
--------------------------------------------------------------------------------

if H_4f_ligands_hybridization_lmct == 1 then
    N_L1 = NewOperator('Number', NFermions, IndexUp_L1, IndexUp_L1, {1, 1, 1, 1, 1, 1, 1})
         + NewOperator('Number', NFermions, IndexDn_L1, IndexDn_L1, {1, 1, 1, 1, 1, 1, 1})

    Delta_4f_L1_i = $Delta(4f,L1)_i_value
    e_4f_i = (28 * Delta_4f_L1_i - 27 * U_4f_4f_i * NElectrons_4f - U_4f_4f_i * NElectrons_4f^2) / (2 * (14 + NElectrons_4f))
    e_L1_i = NElectrons_4f * (-2 * Delta_4f_L1_i + U_4f_4f_i * NElectrons_4f + U_4f_4f_i) / (2 * (NElectrons_4f + 14))

    Delta_4f_L1_m = $Delta(4f,L1)_m_value
    e_4f_m = (28 * Delta_4f_L1_m - U_4f_4f_m * NElectrons_4f^2 - 39 * U_4f_4f_m * NElectrons_4f - 228 * U_3p_4f_m) / (2 * (NElectrons_4f + 20))
    e_3p_m = (28 * Delta_4f_L1_m + U_4f_4f_m * NElectrons_4f^2 + U_4f_4f_m * NElectrons_4f - 2 * U_3p_4f_m * NElectrons_4f^2 - 30 * U_3p_4f_m * NElectrons_4f - 28 * U_3p_4f_m) / (2 * (NElectrons_4f + 20))
    e_L1_m = (-2 * Delta_4f_L1_m * NElectrons_4f - 12 * Delta_4f_L1_m + U_4f_4f_m * NElectrons_4f^2 + U_4f_4f_m * NElectrons_4f + 12 * U_3p_4f_m * NElectrons_4f + 12 * U_3p_4f_m) / (2 * (NElectrons_4f + 20))

    Delta_4f_L1_f = $Delta(4f,L1)_f_value
    e_4f_f = (28 * Delta_4f_L1_f - 27 * U_4f_4f_f * NElectrons_4f - U_4f_4f_f * NElectrons_4f^2) / (2 * (14 + NElectrons_4f))
    e_L1_f = NElectrons_4f * (-2 * Delta_4f_L1_f + U_4f_4f_f * NElectrons_4f + U_4f_4f_f) / (2 * (NElectrons_4f + 14))


    H_i = H_i + Chop(
          e_4f_i * N_4f
        + e_L1_i * N_L1)

    H_m = H_m + Chop(
          e_4f_m * N_4f
        + e_3p_m * N_3p
        + e_L1_m * N_L1)

    H_f = H_f + Chop(
          e_4f_f * N_4f
        + e_L1_f * N_L1)


    B40_L1_i = $B40(L1)_i_value
    B60_L1_i = $B60(L1)_i_value

    Akm_L1_i = {
        {4,  0, B40_L1_i},
        {4, -4, math.sqrt(5/14) * B40_L1_i},
        {4,  4, math.sqrt(5/14) * B40_L1_i},
        {6,  0, B60_L1_i},
        {6, -4, -math.sqrt(7/2) * B60_L1_i},
        {6,  4, -math.sqrt(7/2) * B60_L1_i},
    }

    H_i = H_i + Chop(NewOperator('CF', NFermions, IndexUp_L1, IndexDn_L1, Akm_L1_i))

    B40_L1_m = $B40(L1)_m_value
    B60_L1_m = $B60(L1)_m_value

    Akm_L1_m = {
        {4,  0, B40_L1_m},
        {4, -4, math.sqrt(5/14) * B40_L1_m},
        {4,  4, math.sqrt(5/14) * B40_L1_m},
        {6,  0, B60_L1_m},
        {6, -4, -math.sqrt(7/2) * B60_L1_m},
        {6,  4, -math.sqrt(7/2) * B60_L1_m},
    }

    H_m = H_m + Chop(NewOperator('CF', NFermions, IndexUp_L1, IndexDn_L1, Akm_L1_m))


    B40_L1_f = $B40(L1)_f_value
    B60_L1_f = $B60(L1)_f_value

    Akm_L1_f = {
        {4,  0, B40_L1_f},
        {4, -4, math.sqrt(5/14) * B40_L1_f},
        {4,  4, math.sqrt(5/14) * B40_L1_f},
        {6,  0, B60_L1_f},
        {6, -4, -math.sqrt(7/2) * B60_L1_f},
        {6,  4, -math.sqrt(7/2) * B60_L1_f},
    }

    H_f = H_f + Chop(NewOperator('CF', NFermions, IndexUp_L1, IndexDn_L1, Akm_L1_f))


    -- Mixing of the f-orbitals with the ligands.
    Va2u_4f_L1 = NewOperator('CF', NFermions, IndexUp_L1, IndexDn_L1, IndexUp_4f, IndexDn_4f, PotentialExpandedOnClm('Oh', 3, {1, 0, 0}))
               + NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm('Oh', 3, {1, 0, 0}))

    Vt2u_4f_L1 = NewOperator('CF', NFermions, IndexUp_L1, IndexDn_L1, IndexUp_4f, IndexDn_4f, PotentialExpandedOnClm('Oh', 3, {0, 1, 0}))
               + NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm('Oh', 3, {0, 1, 0}))

    Vt1u_4f_L1 = NewOperator('CF', NFermions, IndexUp_L1, IndexDn_L1, IndexUp_4f, IndexDn_4f, PotentialExpandedOnClm('Oh', 3, {0, 0, 1}))
               + NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm('Oh', 3, {0, 0, 1}))

    Va2u_4f_L1_i = $Va2u(4f,L1)_i_value
    Vt2u_4f_L1_i = $Vt2u(4f,L1)_i_value
    Vt1u_4f_L1_i = $Vt1u(4f,L1)_i_value

    Va2u_4f_L1_m = $Va2u(4f,L1)_m_value
    Vt2u_4f_L1_m = $Vt2u(4f,L1)_m_value
    Vt1u_4f_L1_m = $Vt1u(4f,L1)_m_value

    Va2u_4f_L1_f = $Va2u(4f,L1)_f_value
    Vt2u_4f_L1_f = $Vt2u(4f,L1)_f_value
    Vt1u_4f_L1_f = $Vt1u(4f,L1)_f_value


    H_i = H_i + Chop(
        Va2u_4f_L1_i * Va2u_4f_L1
      + Vt2u_4f_L1_i * Vt2u_4f_L1
      + Vt1u_4f_L1_i * Vt1u_4f_L1)

    H_m = H_m + Chop(
        Va2u_4f_L1_m * Va2u_4f_L1
      + Vt2u_4f_L1_m * Vt2u_4f_L1
      + Vt1u_4f_L1_m * Vt1u_4f_L1)

    H_f = H_f + Chop(
        Va2u_4f_L1_f * Va2u_4f_L1
      + Vt2u_4f_L1_f * Vt2u_4f_L1
      + Vt1u_4f_L1_f * Vt1u_4f_L1)
end


--------------------------------------------------------------------------------
-- Define the magnetic field and exchange field terms.
--------------------------------------------------------------------------------
Sx_4f = NewOperator('Sx', NFermions, IndexUp_4f, IndexDn_4f)
Sy_4f = NewOperator('Sy', NFermions, IndexUp_4f, IndexDn_4f)
Sz_4f = NewOperator('Sz', NFermions, IndexUp_4f, IndexDn_4f)
Ssqr_4f = NewOperator('Ssqr', NFermions, IndexUp_4f, IndexDn_4f)
Splus_4f = NewOperator('Splus', NFermions, IndexUp_4f, IndexDn_4f)
Smin_4f = NewOperator('Smin', NFermions, IndexUp_4f, IndexDn_4f)

Lx_4f = NewOperator('Lx', NFermions, IndexUp_4f, IndexDn_4f)
Ly_4f = NewOperator('Ly', NFermions, IndexUp_4f, IndexDn_4f)
Lz_4f = NewOperator('Lz', NFermions, IndexUp_4f, IndexDn_4f)
Lsqr_4f = NewOperator('Lsqr', NFermions, IndexUp_4f, IndexDn_4f)
Lplus_4f = NewOperator('Lplus', NFermions, IndexUp_4f, IndexDn_4f)
Lmin_4f = NewOperator('Lmin', NFermions, IndexUp_4f, IndexDn_4f)

Jx_4f = NewOperator('Jx', NFermions, IndexUp_4f, IndexDn_4f)
Jy_4f = NewOperator('Jy', NFermions, IndexUp_4f, IndexDn_4f)
Jz_4f = NewOperator('Jz', NFermions, IndexUp_4f, IndexDn_4f)
Jsqr_4f = NewOperator('Jsqr', NFermions, IndexUp_4f, IndexDn_4f)
Jplus_4f = NewOperator('Jplus', NFermions, IndexUp_4f, IndexDn_4f)
Jmin_4f = NewOperator('Jmin', NFermions, IndexUp_4f, IndexDn_4f)

Tx_4f = NewOperator('Tx', NFermions, IndexUp_4f, IndexDn_4f)
Ty_4f = NewOperator('Ty', NFermions, IndexUp_4f, IndexDn_4f)
Tz_4f = NewOperator('Tz', NFermions, IndexUp_4f, IndexDn_4f)

Sx = Sx_4f
Sy = Sy_4f
Sz = Sz_4f

Lx = Lx_4f
Ly = Ly_4f
Lz = Lz_4f

Jx = Jx_4f
Jy = Jy_4f
Jz = Jz_4f

Tx = Tx_4f
Ty = Ty_4f
Tz = Tz_4f

Ssqr = Sx * Sx + Sy * Sy + Sz * Sz
Lsqr = Lx * Lx + Ly * Ly + Lz * Lz
Jsqr = Jx * Jx + Jy * Jy + Jz * Jz

if H_magnetic_field == 1 then
    Bx_i = $Bx_i_value
    By_i = $By_i_value
    Bz_i = $Bz_i_value

    Bx_m = $Bx_m_value
    By_m = $By_m_value
    Bz_m = $Bz_m_value

    Bx_f = $Bx_f_value
    By_f = $By_f_value
    Bz_f = $Bz_f_value


    H_i = H_i + Chop(
          Bx_i * (2 * Sx + Lx)
        + By_i * (2 * Sy + Ly)
        + Bz_i * (2 * Sz + Lz))

    H_m = H_m + Chop(
          Bx_m * (2 * Sx + Lx)
        + By_m * (2 * Sy + Ly)
        + Bz_m * (2 * Sz + Lz))
        
    H_f = H_f + Chop(
          Bx_f * (2 * Sx + Lx)
        + By_f * (2 * Sy + Ly)
        + Bz_f * (2 * Sz + Lz))
end

if H_exchange_field == 1 then
    Hx_i = $Hx_i_value
    Hy_i = $Hy_i_value
    Hz_i = $Hz_i_value

    Hx_m = $Hx_m_value
    Hy_m = $Hy_m_value
    Hz_m = $Hz_m_value

    Hx_f = $Hx_f_value
    Hy_f = $Hy_f_value
    Hz_f = $Hz_f_value


    H_i = H_i + Chop(
          Hx_i * Sx
        + Hy_i * Sy
        + Hz_i * Sz)
    
    H_m = H_m + Chop(
          Hx_m * Sx
        + Hy_m * Sy
        + Hz_m * Sz)

    H_f = H_f + Chop(
          Hx_f * Sx
        + Hy_f * Sy
        + Hz_f * Sz)
end

NConfigurations = $NConfigurations

--------------------------------------------------------------------------------
-- Define the restrictions and set the number of initial states.
--------------------------------------------------------------------------------
InitialRestrictions = {NFermions, NBosons, {'111111 00000000000000', NElectrons_3p, NElectrons_3p},
                                           {'000000 11111111111111', NElectrons_4f, NElectrons_4f}}

IntermediateRestrictions = {NFermions, NBosons, {'111111 00000000000000', NElectrons_3p - 1, NElectrons_3p - 1},
                                         {'000000 11111111111111', NElectrons_4f + 1, NElectrons_4f + 1}}

FinalRestrictions = InitialRestrictions


if H_4f_ligands_hybridization_lmct == 1 then
    InitialRestrictions = {NFermions, NBosons, {'111111 00000000000000 00000000000000', NElectrons_3p, NElectrons_3p},
                                               {'000000 11111111111111 00000000000000', NElectrons_4f, NElectrons_4f},
                                               {'000000 00000000000000 11111111111111', NElectrons_L1, NElectrons_L1}}
    
    IntermediateRestrictions = {NFermions, NBosons, {'111111 00000000000000 00000000000000', NElectrons_3p - 1, NElectrons_3p - 1},
                                             {'000000 11111111111111 00000000000000', NElectrons_4f + 1, NElectrons_4f + 1},
                                             {'000000 00000000000000 11111111111111', NElectrons_L1, NElectrons_L1}}

    FinalRestrictions = InitialRestrictions
                                                   
    CalculationRestrictions = {NFermions, NBosons, {'000000 00000000000000 11111111111111', NElectrons_L1 - (NConfigurations - 1), NElectrons_L1}}
end

T = $T * EnergyUnits.Kelvin.value

-- Approximate machine epsilon for single precision arithmetics.
epsilon = 1.19e-07

NPsis = $NPsis
NPsisAuto = $NPsisAuto

dZ = {}

if NPsisAuto == 1 and NPsis ~= 1 then
    NPsis = 4
    NPsisIncrement = 8
    NPsisIsConverged = false

    while not NPsisIsConverged do
        if CalculationRestrictions == nil then
            Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis)
        else
            Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis, {{'restrictions', CalculationRestrictions}})
        end

        if not (type(Psis_i) == 'table') then
            Psis_i = {Psis_i}
        end

        E_gs_i = Psis_i[1] * H_i * Psis_i[1]

        Z = 0

        for i, Psi in ipairs(Psis_i) do
            E = Psi * H_i * Psi

            if math.abs(E - E_gs_i) < epsilon then
                dZ[i] = 1
            else
                dZ[i] = math.exp(-(E - E_gs_i) / T)
            end

            Z = Z + dZ[i]

            if (dZ[i] / Z) < math.sqrt(epsilon) then
                i = i - 1
                NPsisIsConverged = true
                NPsis = i
                Psis_i = {unpack(Psis_i, 1, i)}
                dZ = {unpack(dZ, 1, i)}
                break
            end
        end

        if NPsisIsConverged then
            break
        else
            NPsis = NPsis + NPsisIncrement
        end
    end
else
    if CalculationRestrictions == nil then
        Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis)
    else
        Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis, {{'restrictions', CalculationRestrictions}})
    end

    if not (type(Psis_i) == 'table') then
        Psis_i = {Psis_i}
    end
        E_gs_i = Psis_i[1] * H_i * Psis_i[1]

    Z = 0

    for i, Psi in ipairs(Psis_i) do
        E = Psi * H_i * Psi

        if math.abs(E - E_gs_i) < epsilon then
            dZ[i] = 1
        else
            dZ[i] = math.exp(-(E - E_gs_i) / T)
        end

        Z = Z + dZ[i]
    end
end

-- Normalize dZ to unity.
for i in ipairs(dZ) do
    dZ[i] = dZ[i] / Z
end

--------------------------------------------------------------------------------
-- Define some helper function for the spectra calculation.
--------------------------------------------------------------------------------

function ValueInTable(value, table)
    -- Check if a value is in a table.
    for k, v in ipairs(table) do
        if value == v then
            return true
        end
    end
    return false
end

--------------------------------------------------------------------------------
-- Define the transition operators.
--------------------------------------------------------------------------------
t = math.sqrt(1/2)

Txy_3p_4f   = NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, IndexUp_3p, IndexDn_3p, {{2, -2, t * I}, {2, 2, -t * I}})
Txz_3p_4f   = NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, IndexUp_3p, IndexDn_3p, {{2, -1, t    }, {2, 1, -t    }})
Tyz_3p_4f   = NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, IndexUp_3p, IndexDn_3p, {{2, -1, t * I}, {2, 1,  t * I}})
Tx2y2_3p_4f = NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, IndexUp_3p, IndexDn_3p, {{2, -2, t    }, {2, 2,  t    }})
Tz2_3p_4f   = NewOperator('CF', NFermions, IndexUp_4f, IndexDn_4f, IndexUp_3p, IndexDn_3p, {{2,  0, 1    }                })

Txy_4f_3p   = NewOperator('CF', NFermions, IndexUp_3p, IndexDn_3p, IndexUp_4f, IndexDn_4f, {{2, -2, t * I}, {2, 2, -t * I}})
Txz_4f_3p   = NewOperator('CF', NFermions, IndexUp_3p, IndexDn_3p, IndexUp_4f, IndexDn_4f, {{2, -1, t    }, {2, 1, -t    }})
Tyz_4f_3p   = NewOperator('CF', NFermions, IndexUp_3p, IndexDn_3p, IndexUp_4f, IndexDn_4f, {{2, -1, t * I}, {2, 1,  t * I}})
Tx2y2_4f_3p = NewOperator('CF', NFermions, IndexUp_3p, IndexDn_3p, IndexUp_4f, IndexDn_4f, {{2, -2, t    }, {2, 2,  t    }})
Tz2_4f_3p   = NewOperator('CF', NFermions, IndexUp_3p, IndexDn_3p, IndexUp_4f, IndexDn_4f, {{2,  0, 1    }                })

k1  = $k1
ev1 = $eps11
eh1 = $eps12

k2  = $k2
ev2 = $eps21
eh2 = $eps22

-- Calculate the right and left polarization vectors.
er1 = {t * (eh1[1] - I * ev1[1]),
       t * (eh1[2] - I * ev1[2]),
       t * (eh1[3] - I * ev1[3])}

el1 = {-t * (eh1[1] + I * ev1[1]),
       -t * (eh1[2] + I * ev1[2]),
       -t * (eh1[3] + I * ev1[3])}

er2 = {t * (eh2[1] - I * ev2[1]),
       t * (eh2[2] - I * ev2[2]),
       t * (eh2[3] - I * ev2[3])}

el2 = {-t * (eh2[1] + I * ev2[1]),
       -t * (eh2[2] + I * ev2[2]),
       -t * (eh2[3] + I * ev2[3])}


function CalculateT_3p_4f(e, k)
    -- Calculate the transition operator for arbitrary polarization.
    T = (e[1] * k[2] + e[2] * k[1]) * Txy_3p_4f
      + (e[1] * k[3] + e[3] * k[1]) * Txz_3p_4f
      + (e[2] * k[3] + e[3] * k[2]) * Tyz_3p_4f
      + (e[1] * k[1] + e[2] * k[2]) * Tx2y2_3p_4f
      + e[3] * k[3] * Tz2_3p_4f
    return Chop(T)
end

function CalculateT_4f_3p(e, k)
    -- Calculate the transition operator for arbitrary polarization.
    T = (e[1] * k[2] + e[2] * k[1]) * Txy_4f_3p
      + (e[1] * k[3] + e[3] * k[1]) * Txz_4f_3p
      + (e[2] * k[3] + e[3] * k[2]) * Tyz_4f_3p
      + (e[1] * k[1] + e[2] * k[2]) * Tx2y2_4f_3p
      + e[3] * k[3] * Tz2_4f_3p
    return Chop(T)
end

Tv_3p_4f = CalculateT_3p_4f(ev1, k1)
Th_3p_4f = CalculateT_3p_4f(eh1, k1)
Tr_3p_4f = CalculateT_3p_4f(er1, k1)
Tl_3p_4f = CalculateT_3p_4f(el1, k1)
Tk_3p_4f = CalculateT_3p_4f(k1, k1)

Tv_4f_3p = CalculateT_4f_3p(ev2, k2)
Th_4f_3p = CalculateT_4f_3p(eh2, k2)
Tr_4f_3p = CalculateT_4f_3p(er2, k2)
Tl_4f_3p = CalculateT_4f_3p(el2, k2)
Tk_4f_3p = CalculateT_4f_3p(k2, k2)


-- List with the user selected spectra.
spectra = {$spectra}

--------------------------------------------------------------------------------
-- Calculate and save the spectrum.
--------------------------------------------------------------------------------
k  = $k
Sk = Chop(k[1] * Sx + k[2] * Sy + k[3] * Sz)
Lk = Chop(k[1] * Lx + k[2] * Ly + k[3] * Lz)
Jk = Chop(k[1] * Jx + k[2] * Jy + k[3] * Jz)
Tk = Chop(k[1] * Tx + k[2] * Ty + k[3] * Tz)

Operators = {H_i, Ssqr, Lsqr, Jsqr, Sk, Lk, Jk, Tk, ldots_4f, N_3p, N_4f, 'dZ'}
header = 'Analysis of the initial Hamiltonian:\n'
header = header .. '=================================================================================================================================\n'
header = header .. 'State         <E>     <S^2>     <L^2>     <J^2>      <Sk>      <Lk>      <Jk>      <Tk>     <l.s>    <N_3p>    <N_4f>          dZ\n'
header = header .. '=================================================================================================================================\n'
footer = '=================================================================================================================================\n'

if H_4f_ligands_hybridization_lmct == 1 then
    Operators = {H_i, Ssqr, Lsqr, Jsqr, Sk, Lk, Jk, Tk, ldots_4f, N_3p, N_4f, N_L1, 'dZ'}
    header = 'Analysis of the initial Hamiltonian:\n'
    header = header .. '===========================================================================================================================================\n'
    header = header .. 'State         <E>     <S^2>     <L^2>     <J^2>      <Sk>      <Lk>      <Jk>      <Tk>     <l.s>    <N_3p>    <N_4f>    <N_L1>          dZ\n'
    header = header .. '===========================================================================================================================================\n'
    footer = '===========================================================================================================================================\n'
end

io.write(header)
for i, Psi in ipairs(Psis_i) do
    io.write(string.format('%5d', i))
    for j, Operator in ipairs(Operators) do
        if j == 1 then
            io.write(string.format('%12.6f', Complex.Re(Psi * Operator * Psi)))
        elseif Operator == 'dZ' then
            io.write(string.format('%12.2E', dZ[i]))
        else
            io.write(string.format('%10.4f', Complex.Re(Psi * Operator * Psi)))
        end
    end
    io.write('\n')
end
io.write(footer)


if next(spectra) == nil then
    return
end

E_gs_i = Psis_i[1] * H_i * Psis_i[1]

if CalculationRestrictions == nil then
    Psis_m = Eigensystem(H_m, IntermediateRestrictions, 1)
else
    Psis_m = Eigensystem(H_m, IntermediateRestrictions, 1, {{'restrictions', CalculationRestrictions}})
end
Psis_m = {Psis_m}
E_gs_m = Psis_m[1] * H_m * Psis_m[1]

Eedge1 = $Eedge1
DeltaE1 = E_gs_m - E_gs_i

Eedge2 = 0.0
DeltaE2 = 0.0

Emin1 = ($Emin1 - Eedge1) + DeltaE1
Emax1 = ($Emax1 - Eedge1) + DeltaE1
NE1 = $NE1
Gamma1 = $Gamma1

Emin2 = ($Emin2 - Eedge2) + DeltaE2
Emax2 = ($Emax2 - Eedge2) + DeltaE2
NE2 = $NE2
Gamma2 = $Gamma2

DenseBorder = $DenseBorder

G = 0


spectrum = 'Isotropic'
if ValueInTable(spectrum, spectra) then
    if CalculationRestrictions == nil then
        G = G + CreateResonantSpectra(H_m, H_f, {Txy_3p_4f, Txz_3p_4f, Tyz_3p_4f, Tx2y2_3p_4f, Tz2_3p_4f}, {Txy_4f_3p, Txz_4f_3p, Tyz_4f_3p, Tx2y2_4f_3p, Tz2_4f_3p}, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'DenseBorder', DenseBorder}})
    else
        G = G + CreateResonantSpectra(H_m, H_f, {Txy_3p_4f, Txz_3p_4f, Tyz_3p_4f, Tx2y2_3p_4f, Tz2_3p_4f}, {Txy_4f_3p, Txz_4f_3p, Tyz_4f_3p, Tx2y2_4f_3p, Tz2_4f_3p}, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'restrictions1', CalculationRestrictions}, {'restrictions2', CalculationRestrictions}, {'DenseBorder', DenseBorder}})
    end

    Giso = 0
    shift = 0
    for i = 1, #Psis_i do
        for j = 1, 3 * 3 do
            Indexes = {}
            for k = 1, NE1 + 1 do
                table.insert(Indexes, k + shift)
            end
            Giso = Giso + Spectra.Element(G, Indexes) * dZ[i]
            shift = shift + NE1 + 1
        end
    end
    Giso = -1 * Giso

    io.write('Here starts Giso spectrum:')
    print(Giso)
    io.write('Here ends Giso spectrum')
end

spectrum = 'Circular Dichroism'
if ValueInTable(spectrum, spectra) then
    if CalculationRestrictions == nil then
        Grv_ = G + CreateResonantSpectra(H_m, H_f, Tr_3p_4f, Tv_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'DenseBorder', DenseBorder}})
        Grh_ = G + CreateResonantSpectra(H_m, H_f, Tr_3p_4f, Th_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'DenseBorder', DenseBorder}})
        Glv_ = G + CreateResonantSpectra(H_m, H_f, Tl_3p_4f, Tv_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'DenseBorder', DenseBorder}})
        Glh_ = G + CreateResonantSpectra(H_m, H_f, Tl_3p_4f, Th_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'DenseBorder', DenseBorder}})
    else
        Grv_ = G + CreateResonantSpectra(H_m, H_f, Tr_3p_4f, Tv_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'restrictions1', CalculationRestrictions}, {'restrictions2', CalculationRestrictions}, {'DenseBorder', DenseBorder}})
        Grh_ = G + CreateResonantSpectra(H_m, H_f, Tr_3p_4f, Th_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'restrictions1', CalculationRestrictions}, {'restrictions2', CalculationRestrictions}, {'DenseBorder', DenseBorder}})
        Glv_ = G + CreateResonantSpectra(H_m, H_f, Tl_3p_4f, Tv_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'restrictions1', CalculationRestrictions}, {'restrictions2', CalculationRestrictions}, {'DenseBorder', DenseBorder}})
        Glh_ = G + CreateResonantSpectra(H_m, H_f, Tl_3p_4f, Th_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'restrictions1', CalculationRestrictions}, {'restrictions2', CalculationRestrictions}, {'DenseBorder', DenseBorder}})
    end

    Grv = 0
    shift = 0
    for i = 1, #Psis_i do
        for j = 1, 1 do
            Indexes = {}
            for k = 1, NE1 + 1 do
                table.insert(Indexes, k + shift)
            end
            Grv = Grv + Spectra.Element(Grv_, Indexes) * dZ[i]
            shift = shift + NE1 + 1
        end
    end
    Grv = -1 * Grv

    Grh = 0
    shift = 0
    for i = 1, #Psis_i do
        for j = 1, 1 do
            Indexes = {}
            for k = 1, NE1 + 1 do
                table.insert(Indexes, k + shift)
            end
            Grh = Grh + Spectra.Element(Grh_, Indexes) * dZ[i]
            shift = shift + NE1 + 1
        end
    end
    Grh = -1 * Grh

    Glv = 0
    shift = 0
    for i = 1, #Psis_i do
        for j = 1, 1 do
            Indexes = {}
            for k = 1, NE1 + 1 do
                table.insert(Indexes, k + shift)
            end
            Glv = Glv + Spectra.Element(Glv_, Indexes) * dZ[i]
            shift = shift + NE1 + 1
        end
    end
    Glv = -1 * Glv

    Glh = 0
    shift = 0
    for i = 1, #Psis_i do
        for j = 1, 1 do
            Indexes = {}
            for k = 1, NE1 + 1 do
                table.insert(Indexes, k + shift)
            end
            Glh = Glh + Spectra.Element(Glh_, Indexes) * dZ[i]
            shift = shift + NE1 + 1
        end
    end
    Glh = -1 * Glh

    io.write('Here starts Grv spectrum:')
    print(Grv)
    io.write('Here ends Grv spectrum')

    io.write('Here starts Grh spectrum:')
    print(Grh)
    io.write('Here ends Grh spectrum')

    io.write('Here starts Glv spectrum:')
    print(Glv)
    io.write('Here ends Glv spectrum')

    io.write('Here starts Glh spectrum:')
    print(Glh)
    io.write('Here ends Glh spectrum')
end

spectrum = 'Linear Dichroism'
if ValueInTable(spectrum, spectra) then
    if CalculationRestrictions == nil then
        Gvv_ = G + CreateResonantSpectra(H_m, H_f, Tv_3p_4f, Tv_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'DenseBorder', DenseBorder}})
        Gvh_ = G + CreateResonantSpectra(H_m, H_f, Tv_3p_4f, Th_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'DenseBorder', DenseBorder}})
        Ghv_ = G + CreateResonantSpectra(H_m, H_f, Th_3p_4f, Tv_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'DenseBorder', DenseBorder}})
        Ghh_ = G + CreateResonantSpectra(H_m, H_f, Th_3p_4f, Th_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'DenseBorder', DenseBorder}})
    else
        Gvv_ = G + CreateResonantSpectra(H_m, H_f, Tv_3p_4f, Tv_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'restrictions1', CalculationRestrictions}, {'restrictions2', CalculationRestrictions}, {'DenseBorder', DenseBorder}})
        Gvh_ = G + CreateResonantSpectra(H_m, H_f, Tv_3p_4f, Th_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'restrictions1', CalculationRestrictions}, {'restrictions2', CalculationRestrictions}, {'DenseBorder', DenseBorder}})
        Ghv_ = G + CreateResonantSpectra(H_m, H_f, Th_3p_4f, Tv_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'restrictions1', CalculationRestrictions}, {'restrictions2', CalculationRestrictions}, {'DenseBorder', DenseBorder}})
        Ghh_ = G + CreateResonantSpectra(H_m, H_f, Th_3p_4f, Th_4f_3p, Psis_i, {{'Emin1', Emin1}, {'Emax1', Emax1}, {'NE1', NE1}, {'Gamma1', Gamma1}, {'Emin2', Emin2}, {'Emax2', Emax2}, {'NE2', NE2}, {'Gamma2', Gamma2}, {'restrictions1', CalculationRestrictions}, {'restrictions2', CalculationRestrictions}, {'DenseBorder', DenseBorder}})
    end

    Gvv = 0
    shift = 0
    for i = 1, #Psis_i do
        for j = 1, 1 do
            Indexes = {}
            for k = 1, NE1 + 1 do
                table.insert(Indexes, k + shift)
            end
            Gvv = Gvv + Spectra.Element(Gvv_, Indexes) * dZ[i]
            shift = shift + NE1 + 1
        end
    end
    Gvv = -1 * Gvv

    Gvh = 0
    shift = 0
    for i = 1, #Psis_i do
        for j = 1, 1 do
            Indexes = {}
            for k = 1, NE1 + 1 do
                table.insert(Indexes, k + shift)
            end
            Gvh = Gvh + Spectra.Element(Gvh_, Indexes) * dZ[i]
            shift = shift + NE1 + 1
        end
    end
    Gvh = -1 * Gvh

    Ghv = 0
    shift = 0
    for i = 1, #Psis_i do
        for j = 1, 1 do
            Indexes = {}
            for k = 1, NE1 + 1 do
                table.insert(Indexes, k + shift)
            end
            Ghv = Ghv + Spectra.Element(Ghv_, Indexes) * dZ[i]
            shift = shift + NE1 + 1
        end
    end
    Ghv = -1 * Ghv

    Ghh = 0
    shift = 0
    for i = 1, #Psis_i do
        for j = 1, 1 do
            Indexes = {}
            for k = 1, NE1 + 1 do
                table.insert(Indexes, k + shift)
            end
            Ghh = Ghh + Spectra.Element(Ghh_, Indexes) * dZ[i]
            shift = shift + NE1 + 1
        end
    end
    Ghh = -1 * Ghh

    io.write('Here starts Gvv spectrum:')
    print(Gvv)
    io.write('Here ends Gvv spectrum')

    io.write('Here starts Gvh spectrum:')
    print(Gvh)
    io.write('Here ends Gvh spectrum')

    io.write('Here starts Ghv spectrum:')
    print(Ghv)
    io.write('Here ends Ghv spectrum')

    io.write('Here starts Ghh spectrum:')
    print(Ghh)
    io.write('Here ends Ghh spectrum')
end