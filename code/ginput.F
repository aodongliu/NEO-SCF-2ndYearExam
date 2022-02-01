*Deck GInput
      Subroutine GInput(V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Copyright (c) 1988,1990,1992,1993,1995,1998,2003,2009,2016,2021,
C                Gaussian, Inc.  All Rights Reserved.
C
C     This is part of the private, development version of
C     the Gaussian(R) DV system of programs. It is based on
C     the Gaussian(R) 16 system (copyright 2016, Gaussian, Inc.),
C     the Gaussian(R) 09 system (copyright 2009, Gaussian, Inc.),
C     the Gaussian(R) 03 system (copyright 2003, Gaussian, Inc.),
C     the Gaussian(R) 98 system (copyright 1998, Gaussian, Inc.),
C     the Gaussian(R) 94 system (copyright 1995, Gaussian, Inc.),
C     the Gaussian 92(TM) system (copyright 1992, Gaussian, Inc.),
C     the Gaussian 90(TM) system (copyright 1990, Gaussian, Inc.),
C     the Gaussian 88(TM) system (copyright 1988, Gaussian, Inc.),
C     the Gaussian 86(TM) system (copyright 1986, Carnegie Mellon
C     University), and the Gaussian 82(TM) system (copyright 1983,
C     Carnegie Mellon University). Gaussian is a federally registered
C     trademark of Gaussian, Inc.
C
C     This software contains proprietary and confidential information,
C     including trade secrets, belonging to Gaussian, Inc.
C
C     This software is provided under written license and may be
C     used, copied, transmitted, or stored only in accord with that
C     written license.
C
C     The following legend is applicable only to US Government contracts
C     under DFARS:
C
C                        RESTRICTED RIGHTS LEGEND
C
C     Use, duplication or disclosure by the US Government is subject to
C     restrictions as set forth in subparagraph (c)(1)(ii) of the Rights
C     in Technical Data and Computer Software clause at DFARS
C     252.227-7013.
C
C     Gaussian, Inc., 340 Quinnipiac St., Bldg. 40, Wallingford CT 06492
C
C     The following legend is applicable only to US Government contracts
C     under FAR:
C
C                        RESTRICTED RIGHTS LEGEND
C
C     Use, reproduction and disclosure by the US Government is subject
C     to restrictions as set forth in subparagraph (c) of the Commercial
C     Computer Software - Restricted Rights clause at FAR 52.227-19.
C
C     Gaussian, Inc., 340 Quinnipiac St., Bldg. 40, Wallingford CT 06492
C
C     LINK 301
C
C1OV3
C     GAUSSIAN SYSTEM INTEGRAL PACKAGE.
C     OVERLAY THREE CONSISTS OF THE NECESSARY PROGRAMS TO EVALUATE
C     THE ONE- AND TWO-ELECTRON INTEGRALS REQUIRED FOR AN SCF
C     CALCULATION.  See Link1 for a list of the overlay 3 links.
C
C     NOTE THAT THERE IS SOME OVERLAP IN FUNCTIONALITY BETWEEN THE
C     TWO-ELECTRON INTEGRAL PROGRAMS LISTED ABOVE.  IN REGULAR
C     CALCULATIONS, ONE WOULD NORMALLY USE LINKS 301, 302, 303,
C     311 AND 314 TO OBTAIN ALL THE NECESSARY INTEGRALS FOR BASIS
C     FUNCTIONS UP (AND INCLUDING) F-FUNCTIONS.  IF D and F-FUNCTIONS
C     ARE NOT PRESENT IN THE ROUTE, L314 CAN BE OMITTED.
C2
C
C2IOP(5)
C2IOP(6)
C     IOP(5) ... TYPE OF BASIS SET.  The same numbers are used for all
C                basis sets, whether intended for use in expanding AOs
C                (IOp(5)) or in expanding the density (IOp(82)).
C    -1  Same as 0.
C     0  MINIMAL  STO-2G TO STO-6G
C     1  EXTENDED  4-31G,5-31G,6-31G
C     2  MINIMAL STO-NG   (VALENCE FUNCTIONS ONLY)
C     3  EXTENDED LP-N1G (VALENCE BASIS FOR CORELESS HARTREE-FOCK
C        PSEUDOPOTENTIALS)
C     4  EXTENDED 6-311G (UMP2 FROZEN CORE OPTIMIZED) BASIS
C        for first row, MacLean-Chandler (12s,9p)-->(631111,52111)
C        for second row.
C        USE IOP(8) TO SELECT 5D/6D.
C     5  SPLIT VALENCE N-21G (OR NN-21G) BASIS FOR FIRST OR
C        SECOND ROW ATOMS.  (VARIOUS IMPLEMENTATIONS MAY OMIT
C        SECOND ROW ATOMS.)  SEE IOP(6) FOR DETERMINATION OF THE
C        NUMBER OF GAUSSIANS IN THE INNER SHELL.
C     6  LANL ECP basis sets.  IOp(6) selects options.
C     7  GENERAL ... SEE ROUTINE GenBas FOR INPUT INSTRUCTIONS.
C     8  Dunning/Caltech basis sets.  Type selected by IOp(6).
C     9  Stevens/Basch/Krauss/Jasien/Cundari ECP basis
C        sets for H-Lu.  Type selected by IOp(6) for H-Ar.
C        Literature citations in CEPPot.
C    10  CBS basis #1 -- 6-31+g(d,p) on H, He
C                        6-311+G(2df) on Li - Ne
C                        6-311+g(3d2f) on Na - Ar
C    11  CBS basis #2 -- 6-31G, use daggers if any polarization
C    12  CBS basis #3 -- 6-311++G(2df,2p) on H - Ne
C                        6-311++g(3d2f) on Na - Ar
C    13  CBS basis #4 -- 6-31+G(d,p) on H - Si
C                        6-31+G(df,p) on P, S, Cl
C    14  CBS basis #5 -- Large APNO basis set
C    15  CBS basis #6 -- Core correlation basis set
C    16  Dunning cc basis sets, type selected by IOp(6)
C        (=0-4 for V{D,T,Q,5,6}Z) and augmented if IOp(7)=10.
C        IOp(6)=5 for MTsmall basis set.
C    17  Stuttgart/Dresden ECP basis sets.  IOp(6) specifies type.
C        Literature citations in SDDPot.
C    18  Ahlrichs SV basis sets.
C    19  Ahlrichs TZV basis sets.
C    20  MIDI! basis sets.
C    21  EPR-II basis sets.
C    22  EPR-III basis sets.
C    23  UGBS basis set.
C    24  G3large basis set.
C    25  G3MP2large basis set.
C    26  Coreless:  Li,Be 2SDF, B-Ne 2MWB, rest LANL1MB.
C    27  DGauss basis sets, selected by IOp(6)
C    28  Auto-generated, useful only for density basis sets.
C    29  Spherical atomic densities:  a single highly contracted
C        s-Gaussian for each atom.  Only useful for fitting sets.
C    30  One s-gaussian per atom; dummy basis used for MM.
C    31  G3largeXP basis set.
C    32  G3MP2largeXP basis set.
C    33  G3 basis 1 - "6-31G(d)" basis set.
C    34  G3 basis 2 - "6-31+G(d)" basis set.
C    35  G3 basis 3 - "6-31G(2df,d)" basis set.
C    36  G4 QZ HF basis.
C    37  G4 5Z HF basis.
C    38  G4MP2 TZ HF basis.
C    39  G4MP2 QZ HF basis.
C    40  Weigand Coulomb fitting set.
C    41  Ahlrichs SVP Coulomb fitting basis.
C    42  Ahlrichs TZVP Coulomb fitting basis.
C    43  Ahlrichs/Weigand def2-SV basis.
C    44  Ahlrichs/Weigand def2-TZV basis.
C    45  Ahlrichs/Weigand QZV basis.
C    46  Fitting set matched to AO basis, or error if there is none.
C        Converted here to matched value.
C    47  Fitting set matched to AO basis, or /Auto if there is none.
C    48  Uncontracted even-tempered basis on protons only.
C
C     IOP(6) ... NUMBER OF GAUSSIAN FUNCTIONS
C     N  STO-NG,N-31G,LP-N1G,STO-NG-VALENCE, N-21G.
C        NOTE IF IOP(5)=3 AND IOP(6)=8 ;LP-31G FOR LI,BE,B,NA,MG,AL
C                                       LP-41G FOR OTHER ROW1 AND TWO
C                                       ATOMS.
C     DEFAULT OPTIONS ... IOP(6)=0
C     IF IOP(5)=0  N=3  STO-3G
C     IF IOP(5)=1  N=4  4-31G
C     IF IOP(5)=2  N=3  STO-3G  (VALENCE)
C     IF IOP(5)=3  N=3
C     IF IOP(5)=5  N=3
C     WHEN IOP(5)=7 (GENERAL BASES), THIS OPTION IS USED TO CONTROL
C     WHERE THE BASIS IS TAKEN FROM:
C     0 ... Read general basis from the input stream.
C     1 ... Read the general basis from the rw-files and merge with
C           the coordinates in /Mol/.
C     2 ... Read the general basis from the checkpoint file.
C     3 ... Same as 1, for density basis (generated here from 1)
C     4 ... Same as 2, for density basis (generated here from 2)
C    1x ... Read from the alternate file and remove functions/ECPs
C           for inactive atoms.  Used for counterpoise calculations,
C           where one wants to modify the basis differently during
C           different steps.
C    2x ... Read from the other alternate file, saved before the basis
C           is massaged, uncontracted, etc.
C     This option is useful when doing general basis geometry
C     optimizations or properties using a wavefunction on the
C     checkpoint file.    If non-standard ECPs are in use, they are
C     read along with the basis set information.
C
C     When IOp(5)=6 (LANL basis and potentials) this selects the type:
C     0 ... LANL1 ECP, MBS.
C     1 ... LANL1 ECP, DZ.
C     2 ... LANL2 ECP (where available, otherwise LANL1), MBS.
C     3 ... LANL2 ECP (where available, otherwise LANL1), DZ.
C
C     When IOp(5)=8 (Dunning bases) this option selects the type:
C     0 ... Dunning full double-zeta.
C     1 ... Dunning valence double-zeta.
C     2 ... WAG basis (Dunning VDZ on first row, SHC ECP on
C           second row).  See Rappe, Smedley, and Goddard,
C           J. Phys. Chem. 85, 1662 (1981) and J. Phys. Chem.
C           85, 3546 (1981).
C
C     When IOp(5)=9 (CEP basis) this option selects the type (H-Ar only):
C     0 ... CEP-4G.
C     1 ... CEP-31G.
C     2 ... CEP-121G.
C
C     When IOp(5)=17 (Stuttgart/Dresden ECP bases) this option selects
C     the type according to:
C     6 ... SDDAll: SDD for Z > 2
C     7 ... SDD for Z > 18 with SEG basis for Lanthanides & Actinides,
C           D95 or 6-31G and no ECP otherwise
C     8 ... SDDOld: same as SDD with old Lanthanide & Actinide basis
C
C     When IOp(5)=26 (Coreless basis) this selects the choice of basis
C     (the same ECPs are used regardless):
C     0 ... Default (3)
C     1 ... Primitives which match the ECPs.
C     2 ... Functions from extended Huckel theory.
C     3 ... VSTO-4G basis for 1st row, along with LP-31G potential.
C   N>3 ... Huckel basis for method N-1
C
C     When IOp(5)=27 (DGauss basis sets):
C     1 ... DGDZVP
C     2 ... DZVP2
C     3 ... DGTZVP
C     4 ... DGA1 (fitting basis)
C     5 ... DGA2 (fitting basis)
C
C     When IOp(5)=48 (Even tempered uncontracted basis):
C     MMMMMMMEE ... Tightest exponent is (MMMMMMM/1000000)*10**EE,
C                   or (MMMMMMM/1000000)*10**(-EE) if negative.
C
C2IOP(7)
C     IOP(7) ... DIFFUSE AND POLARIZATION FUNCTIONS.
C          0 ... NONE.
C          1 ... D-FUNCTIONS ON HEAVY ATOMS (2ND ROW ONLY FOR
C                3-21G).
C          2 ... 2 D-FUNCTIONS ON HEAVY ATOMS (SCALED UP AND DOWN
C                BY A FACTOR OF 2 FROM THE STANDARD SINGLE-D VALUE).
C          3 ... ONE SET OF D-FUNCTIONS AND ONE SET OF F-FUNCTIONS
C                ON HEAVY ATOMS.
C                (indicates an extra tight 2df with ccp basis sets.
C          4 ... TWO SETS OF D-FUNCTIONS AND ONE SET OF F-FUNCTIONS
C                ON HEAVY ATOMS.
C          5 ... Three sets of d functions.
C          6 ... Three sets of d functions and one set of f-functions.
C          7 ... Three sets of d functions and two sets of f-functions.
C          8 ... CBS-Q d(f),d,p polarization basis
C          9 ... Tight d for VnZ+1 (W1 theory)
C         10 ... A SET OF DIFFUSE SP FUNCTIONS ON HEAVY ATOMS.
C         20 ... Augment non-hydrogens only (cc basis sets only).
C         30 ... maug-:  Main group(SP), TM(SP).
C         40 ... H(SP), Main group(SP), TM(SP)
C         50 ... Jul- aug:  up to LVal on non-H,He.
C         60 ... Jun- aug:  up to LVal-1 on non-H,He.
C         70 ... May- aug:  up to LVal-2 on non-H,He.
C         80 ... Apr- aug:  up to LVal-3 on non-H,He.
C        100 ... P-FUNCTIONS ON HYDROGENS; interpret first digit as
C                pol level for ugbs.
C        200 ... 2 SETS OF P-FUNCTIONS ON HYDROGENS.
C        300 ... ONE SET OF P-FUNCTIONS AND ONE SET OF D-FUNCTIONS
C                ON HYDROGENS.
C        400 ... TWO SETS OF P-FUNCTIONS AND ONE SET OF D-FUNCTIONS
C                ON HYDROGENS.
C        500 ... Three sets of p-functions.
C        600 ... Three sets of p-functions and one set of d-functions.
C        700 ... (2d,d,p) -- 2d on 2nd and later atoms, 1d on 1st row atoms.
C       1000 ... A diffuse s function on hydrogens (Pople-style basis sets).
C       1000 ... Truhlar-style calendar basis sets (inconsistent s and p
C                diffuse functions).
C       N000 ... Number of times to augment (cc-pvxz basis sets).
C      M0000 ... Maximum L for diffuse functions is L(valence)-M.
C
C     When IOp(5)=48 (Even tempered uncontracted basis):
C     SSSSSNNNLL ... NNN exponents, angular momentum 0 through LL,
C                    scaling by 1/Sqrt(SSSSS/1000).
C
C2IOP(8)
C     IOP(8) ... SELECTION OF PURE/CARTESIAN FUNCTIONS.
C     0 ... SELECTION DETERMINED BY THE BASIS
C           N-31G         ... 6D/7F
C           N-311G        ... 5D/7F
C           N-21G*        ... 5D
C           STO-NG*       ... 5D
C           LP-N1G*       ... 5D
C           LP-N1G**      ... 5D
C           GENERAL BASIS ... 5D/7F
C     1 ... FORCE 5D.
C     2 ... FORCE 6D.
C    10 ... FORCE 7F.
C    20 ... FORCE 10F.
C
C2IOp(9)
C     IOp(9) ... Where 308 should store dipole velocity integrals.
C          0 ... Usual place (572).
C         -1 ... Write over the dipole length integrals (518).
C          N ... Store in RWF N.
C
C2IOp(10)
C     IOp(10) ... Modification of internally stored bases (default 12000):
C           0 ... None.
C           1 ... Read in general basis data in addition to setting up a
C                 standard basis.
C          10 ... Massage the data in Common /B/ and Common /Mol/.
C          20 ... Massage the data in Common /B/ and Common /Mol/, but don't
C                 change ian if nuc charge changed.
C         100 ... Add ghost atoms to /B/ so that every shell is on a separate
C                 center.
C        1000 ... Split S=P AO basis shells into separate S and P shells.
C        2000 ... Do not split S=P AO shells.
C       10000 ... Split S=P=D=... AO shells into S=P, D, F, ...
C       20000 ... Do not split AO S=P=D... shells.
C      100000 ... Uncontract the AO basis and removes duplicate primitives.
C      200000 ... Uncontract the density basis and removes duplicate primitives.
C      300000 ... Uncontract both basis sets and removes duplicate primitives.
C      400000 ... Same as 1 but don't remove the duplicates.
C      500000 ... Same as 2 but don't remove the duplicates.
C      600000 ... Same as 3 but don't remove the duplicates from the AO basis.
C      700000 ... Same as 3 but don't remove the duplicates from the density basis
C      800000 ... Same as 3 but don't remove the duplicates from both basis.
C     1000000 ... Modification 1 for Fermi-contact spin-spin coupling
C     2000000 ... Modification 2 for Fermi-contact spin-spin coupling
C
C3DEFAULTS
C     STO-NG STANDARD SCALE-FACTORS.  For VSTO-nG, the values for H-Ar
C     can be determined by Slater's rules:
C     H=1.2,He=1.7,Li-Ne=0.325*(IA-1),Na-Ar=(0.65*I-4.95)/3
C
C     ATOM   1S   2SP   3SP
C
C     H     1.24
C     HE    1.69
C     LI    2.69  0.80
C     BE    3.68  1.15
C     B     4.68  1.50
C     C     5.67  1.72
C     N     6.67  1.95
C     O     7.66  2.25
C     F     8.65  2.55
C     NE    9.64  2.88
C     NA   10.61  3.48  1.75
C     MG   11.59  3.90  1.70
C     AL   12.56  4.36  1.70
C     SI   13.53  4.83  1.75
C     P    14.50  5.31  1.90
C     S    15.47  5.79  2.05
C     CL   16.43  6.26  2.10
C     A    17.40  6.74  2.33
C
C     INNER SHELLS ARE BEST ATOM VALUES ... J.CHEM.PHYS. 38, 2686 (1963)
C     OUTER SHELL HAS BEEN SELECTED ON THE BASIS OF NUMEROUS
C     OPTIMIZATION STUDIES ON VARIED SMALL MOLECULES
C
C     N-31G (ALSO N-31G* AND N-31G**) STANDARD SCALE-FACTORS
C
C     HYDROGEN
C
C            1S   1S*
C
C     H     1.20  1.15
C
C     FIRST ROW ATOMS
C
C     ATOM   1S   2SP   2SP*
C
C     B     1.00  1.03  1.12
C     C     1.00  1.00  1.04
C     N     1.00  0.99  0.98
C     O     1.00  0.99  0.98
C     F     1.00  1.00  1.00
C
C     SECOND ROW ATOMS
C
C     ATOM   1S   2SP   3SP   3SP*
C
C     P     1.00  1.00  0.98  1.02
C     S     1.00  1.00  0.98  1.01
C     CL    1.00  1.00  1.00  1.01
C
C     LP-N1G  SCALE=1.0 FOR LI-AR (INNER AND OUTER)
C
C
C     STANDARD POLARIZATION EXPONENTS FOR N-31G* AND N-31G** BASES
C
C     ATOM    VALUE
C     ----    -----
C      H       1.1
C      LI      0.2
C      BE      0.4
C      B       0.6
C      C-NE    0.8
C
C
C     STANDARD POLARIZATION EXPONENTS FOR STO-NG* BASIS.
C
C     ATOM    VALUE
C     ----    -----
C     NA, MG   0.09
C     AL-CL    0.39
C
C2IOP(11)
C     IOP(11) ... CONTROL OF TWO-ELECTRON INTEGRAL STORAGE FORMAT.
C     0 ... REGULAR INTEGRAL FORMAT IS USED.
C     1 ... RAFFENETTI '1' INTEGRAL FORMAT IS USED.  CAN ONLY
C           BE USED WITH THE CLOSED SHELL SCF.
C     2 ... RAFFENETTI '2' INTEGRAL FORMAT.  SUITABLE FOR USE WITH
C           THE OPEN SHELL (UHF) SCF.
C     3 ... RAFFENETTI '3' INTEGRAL FORMAT.  SUITABLE FOR USE WITH
C           OPEN SHELL RHF SCF AND THE POST-SCF PROCEDURES, but not
C           yet accepted by them.
C     9 ... USE ILSW TO DECIDE BETWEEN RAFFENETTI 1 AND 2.
C
C2IOp(12)
C     IOP(12) ... Flag for semi-empirical runs, to account for sparkles,
C                 translation vectors and d functions properly:
C           1 ... CNDO
C           2 ... INDO
C           3 ... ZINDO/1
C           4 ... ZINDO/S
C           5 ... MINDO3
C           6 ... MNDO
C           7 ... AM1
C           8 ... PM3
C           9 ... DFTB
C          10 ... PM6
C          11 ... PDDG
C
C2IOp(13)
C
C     IOp(13) ... Nuclear center whose Fermi contact terms are
C                 to be added to the core hamiltonian.  The
C                 magnitude is specified by IOp(15).
C
C2IOP(14)
C
C     IOP(14) ... Addition of electrostatic integrals to core
C                 hamiltonian.
C           0 ... No.
C         -1x ... SCRF calculation -- multiply moments by fudge factor
C                 for charged species.
C          -7 ... Same as 0.
C          -6 ... Read coefficients of field, starting with electric
C                 field, up through 34 elements (hexadecapoles) in
C                 free format, blank terminated.
C          -5 ... Read components of electric field only from /Gen/ on
C                 checkpoint file.
C          -4 ... Read components of moments off rwf 521 on chk file.
C          -3 ... Read components of electric field only from /Gen/.
C          -2 ... Read components of moments off rwf 521, or set to zero
C                 if the file is not present.
C          -1 ... Yes, read 12 cards with x,y,z components of electric
C                 field, followed by xx,yy,zz,xy,xz,yz electric field
C                 gradient, xxx,yyy,zzz,xyy,xxy,xxz,xzz,yzz,yyz,xyz
C                 field second derivatives, and xxxx,yyyy,zzzz,xxxy,
C                 xxxz,yyyx,yyyz,zzzx,zzzy,xxyy,xxzz,yyzz,xxyz,yyxz,
C                 zzxy field third derivatives in format (3D20.10).
C                 (These correspond to dipole, quadrupole, octopole,
C                 and hexadecapole perturbations).
C        1-34 ... Just component number n in the above order with
C                 magnitude given by IOp(15).
C
C      The nuclear repulsion energy is also modified appropriately, and
C      the electric field is stored in Gen(2-4).
C
C2IOp(15)
C
C     IOp(15) ... Magnitude of electric field.
C           N ... N * 0.0001.
C
C2IOP(16)
C     IOP(16) ... Pseudopotential option
C           0 ... Default.  ECPs if defined with the basis set.
C           1 ... Yes, read if general basis.
C           2 ... No.
C          00 ... Default (10).
C          10 ... Read ECPs for QM atoms.
C          20 ... Read ECPs for EE charge centers only.
C          30 ... Read two input sections, for QM then EE charge centers.
C         000 ... Default (100).
C         100 ... Spin-orbit ECP coefficients are used as-is, appropriate
C                 for published Stuttgart potentials.
C         200 ... Spin-orbit ECP coefficients are scaled by 2/(2l+1),
C                 appropriate for CRENBL potentials.
C
C2IOP(17)
C     IOP(17) ... SPECIFICATION OF PSEUDOPOTENTIALS
C    -2 ... Same as 0.
C    -1 ... Read potential in old format.
C     0 ... Default, based on IOp(5).
C     1 ... USE INTERNALLY STORED 'CORELESS HARTREE-FOCK'
C     2 ... Goddard/Smedley SECE/SHC potentials.
C     3 ... Stevens/Basch/Krauss CEP potentials.
C     4 ... LANL1 potentials.
C     5 ... LANL2 potentials.
C     6-7   unused
C     8 ... READ IN FROM CARDS (SEE PINPUT FOR DETAILS)
C     9 ... Dresden/Stuttgart potentials - SDD combination
C    10 ... Dresden/Stuttgart potentials - SDD for Z > 18, D95V, no ECP
C                                          otherwise.
C    11 ... Dresden/Stuttgart potentials - SDF
C    12 ... Dresden/Stuttgart potentials - SHF
C    13 ... Dresden/Stuttgart potentials - MDF
C    14 ... Dresden/Stuttgart potentials - MHF (first set)
C    15 ... Dresden/Stuttgart potentials - MHF (second set)
C    16 ... Dresden/Stuttgart potentials - MWB (first set)
C    17 ... Dresden/Stuttgart potentials - MWB (second set)
C    18 ... Dresden/Stuttgart potentials - MWB (third set)
C    19 ... Pseudopotentials for all coreless basis.
C    20 ... Alternative potentials for coreless basis.
C    21 ... Psuedopotentials for the def2SV, def2TZV, and QZV basis sets.
C
C2IOP(18)
C     IOP(18) ... Printing of pseudopotentials
C     0 ... Print only when input is from cards or if GFPrint was
C           specified.
C     1 ... Print
C     2 ... Don't print
C
C2IOP(19)
C     IOP(19) ... SPECIFICATION OF SUBSTITUTION POTENTIAL TYPE
C     0 ... DONT/T USE ANY SUBSTITUTION POTENTIALS
C     N ... REPLACE THE STANDARD POTENTIAL OF THIS RUN (EG.CHF),
C           WITH A SUBSTITUTION POTENTIAL OF TYPE N WHEREVER SUCH
C           A SUBSTITUTION POTENTIAL EXISTS.
C
C2IOp(20)
C
C     IOp(20) ... Size of buffers for integral file.
C     (no longer used, fixed at 131072 wp words).
C
C2IOp(21)
C
C     IOp(21) ... Size of buffers for integral derivative file.
C     (no longer used, fixed at 131072 wp words).
C
C2IOP(22)
C     IOP(22) ... CONTROL OF THE PRE-CUTOFF IN THE TWO-ELECTRON
C                 D-INTEGRAL PROGRAM.  Used only in L312.
C     0 ... NO PRE-CUTOFF.
C     1 ... PRE-CUTOFFS DESIGNED FOR THE 6-31G* BASIS.
C
C2IOp(23)
C     IOp(23) ... Disable use of certain basis functions.
C           0 ... Use all basis functions.
C           1 ... Read in a list of basis function numbers in
C                 Format (10I5), terminated by a blank line,
C                 and set their dialgonal core Hamiltonian
C                 elements to +100.0.
C
C2IOP(24)
C     IOP(24) ... Printing of gaussian function table.
C           0 ... Default (don't print).
C           1 ... Print old-fashioned table.
C          10 ... Print as GenBas input.
C         100 ... Print in more readable format.
C        1000 ... Print shell coordinates.
C       00000 ... Print AO basis using default primitive normalization;
C       10000 ... Print AO basis using coefficients of raw primitives
C       20000 ... Print AO basis using coefficients of AO normalized primitives
C       30000 ... Print AO basis using coefficients of J normalized primitives
C      000000 ... Print density basis using default primitive normalization;
C      100000 ... Print density using coefficients of raw primitives
C      200000 ... Print density using coefficients of AO normalized primitives
C      300000 ... Print density using coefficients of J normalized primitives
C     1000000 ... Print current basis set for all elements as gen input,
C                 reading the name of the file to write from input
C     2000000 ... Print current basis set for all elements as gen input,
C                 reading the name of the file to write from input but
C                 omitting program revision string in the header line.
C
C2IOP(25)
C     IOP(25) ... NUMBER OF LAST TWO-ELECTRON INTEGRAL LINK.
C     -2 ... Use integrals from a previous job ... read /IBF/ from
C            the checkpoint file.
C     -1 ... We are re-using integrals produced earlier in the
C            current calculation ... use the /IBF/ already on the RWF.
C      0 ... WE ARE NOT USING TWO-ELECTRON INTEGRALS.
C      1 ... Direct SCF.
C     >0 ... LINK NUMBER.
C
C2IOP(26)
C    IOP(26) ... ACCURACY OPTION.
C          0 ... DEFAULT.  INTEGRALS ARE COMPUTED TO 10**-10 ACCURACY.
C          1 ... TEST.  DO ALL INTEGRALS AS WELL AS POSSIBLE in L311.
C          2 ... STO-3G.  USE OLD very inaccurate CUTOFFS IN LINK 311.
C         10 ... TEST.  DO ALL INTEGRALS AS WELL AS POSSIBLE in L314.
C         20 ... Sleazy.  Use looser cutoffs in L314.
C
C2IOP(27)
C     IOP(27) ... Computing and storing of small two-electron integrals.
C           0 ... Discard integrals with magnitude less than 10**-12.
C           N ... Discard integrals with magnitude less than 10**-N.
C
C2IOP(28)
C     IOP(28) ... Special SP code control.
C           0 ... Default, use IsAlg.
C           1 ... All integrals with d's -- L311 does nothing.
C           2 ... SP integrals in link 311, d and higher elsewhere.
C           3 ... All integrals done in l314 using Prism.
C
C2IOP(29)
C     IOP(29) ... Accuracy in L302:
C           0 ... Default (10**-13).
C           N ... 10**-N.
C
C2IOP(30)
C     IOP(30) ... CONTROL OF TWO-ELECTRON INTEGRAL SYMMETRY.
C     0 ... TWO-ELECTRON INTEGRAL SYMMETRY IS TURNED OFF.
C     1 ... TWO-ELECTRON INTEGRAL SYMMETRY IS TURNED ON.  NOTE,
C           HOWEVER, THE SET2E WILL INTERROGATE ILSW TO SEE IF
C           THE SYMMETRY RW-FILES EXIST.  IF THEY DON'T, SYMMETRY
C           HAS BEEN TURNED OFF elseWHERE, AND SET2E WILL ALSO
C           TURN IT OFF HERE.
C2IOP(31)
C     IOP(31) ... USE OF SYMMETRY IN COMPUTING GRADIENT (Obsolete).
C
C2IOP(32)
C     IOP(32) ... Whether to check the eigenvalues of the overlap matrix:
C           0 ... Default (205).
C           1 ... Yes.
C           2 ... No.
C           3 ... Yes, and reduce expansion space if linear dependence
C                 is found (NYI).
C           4 ... Yes, and use Schmidt orthogonalization to reduce
C                 expansion space.
C           5 ... Yes, using SVD to reduce expansion space.
C           6 ... Set up SAOs as with 5 but using diagonalization
C                 instead of SVD.
C           9 ... Set up a unit matrix for the transformation.
C         100 ... Try to make the new set of vectors as much like the
C                 previous set, if any.
C         200 ... Do SVD ignoring the previous orthonormal set, if any.
C        1000 ... Use schmidt orthogonalized to match to previous o.n. set.
C        2000 ... Use symmetric orthogonalization with Jacobi
C                 diagonalization to match to previous o.n. set.
C       10000 ... Check orthonormality of generated set in RAOMat.
C       20000 ... Do not check orthonormality of generated set in RAOMat.
C
C2IOP(33)
C     IOP(33) ... INTEGRAL PACKAGE PRINTING.
C     0 ... NO INTEGRALS ARE PRINTED.
C     1 ... PRINT ONE-ELECTRON INTEGRALS.
C     3 ... PRINT TWO-ELECTRON INTEGRALS IN STANDARD FORMAT.
C     4 ... PRINT TWO-ELECTRON INTEGRALS IN DEBUG FORMAT.
C     5 ... COMBINATION OF 1 AND 3.
C     6 ... COMBINATION OF 1 AND 4.
C
C2IOP(34)
C     IOP(34) ... DUMP OPTION.
C     0 ... NO DUMP.
C     1 ... CONTROL WORDS PRINTED (AS USUAL).
C     2 ... ADDITIONALLY, COMMON/B/ IS DUMPED AT THE BEGINNING
C           OF EACH INTEGRAL LINK.
C     3 ... ADDITIONALLY, THE INTEGRALS ARE PRINTED (STANDARD FORMAT).
C
C2IOP(35)
C
C     IOP(35) ... LAST INTEGRAL DERIVATIVE LINK (No longer used in
C                 overlay 3).
C     0 ... WHATEVER LINK STARTS WRITING THE INTEGRAL DERIVATIVE
C           FILE SHOULD ALSO CLOSE IT.
C     N ... IS THE NUMBER OF THE LAST TWO-ELECTRON INTEGRAL DERIVATIVE
C           PROGRAM.
C
C2IOp(36)
C
C     IOp(36) ... Matrices to compute in L303 and L308:
C          -1 ... None
C           0 ... Default (dipole).
C           1 ... Dipole.
C           2 ... Quadrupole.
C           3 ... Octopole.
C           4 ... Hexadecapole.
C          00 ... Default (same as 20).
C          10 ... Do not compute absolute overlaps.
C          20 ... Compute absolute overlap over contracted functions.
C          30 ... Compute absolute overlap over both contracted and
C                 over primitive functions.
C         000 ... Default, same as 100.
C         100 ... L308 should compute (del r + r del) in addition
C                 to Del and r x Del.
C         200 ... L308 should just Del and r x Del.
C
C2IOP(37)
C     IOp(37) ... Sorting of ao 2e integrals in L315;
C                 type of integrals in L320.
C           0 ... Default (No Sort).
C           1 ... Yes, if l315 in the route.
C           2 ... No.
C          00 ... Default (10)
C          10 ... L320 produces integrals.
C          20 ... L320 produces integral derivatives.
C          30 ... L320 produces integrals and integral derivatives.
C
C2IOP(38)
C     IOp(38) ... Algorithm for 1e integrals:
C           0 ... Default in 302, same as 1.
C           1 ... PRISM.
C           2 ... Rys.
C          00 ... Default in 308, same as 1.
C          10 ... PRISM.
C          20 ... Explicit spdf code.
C
C2IOP(39)
C     IOp(39) ... Initialization of force and force constant rwfs.
C           0 ... Initialize.
C           1 ... Leave alone.
C
C2IOp(40)
C   IOp(40) ... Neglect of integrals:
C         0 ... keep all integrals.
C         1 ... neglect four center integrals.
C         2 ... neglect three center two-electron integrals as well.
C         3 ... neglect 2e integrals with diatomic differential overlap.
C        10 ... neglect three center one-electron integrals.
C        20 ... neglect 1e integrals with diatomic differential overlap.
C        30 ... Do only overlap and not other 1e integrals.
C
C2IOp(41)
C
C   IOp(41) ... Various semi-empirical methods.
C         0 ... No NDDO
C         1 ... NDDO
C        00 ... Default use of NDDO beta parameters (arithmetic mean for
C               indo parameters, geometric mean for NDDO/1 or read-in
C               parameters).
C        10 ... Arithmetic mean in NDDO.
C        20 ... Geometric mean in NDDO.
C       000 ... Default parameters (same as 5).
C       100 ... Read parameters for atomic numbers 1-18 in the order
C               Scale (D20.12), followed by ((HDiag(J,I),J=1,3),I=1,18)
C               (Format 3D20.12), followed by ((Beta(J,I),J=1,3),I=1,18)
C       200 ... Read parameters from rwf.
C       300 ... Read parameters from chk.
C       400 ... Original INDO/2 Beta and HDiag Parameters.
C       500 ... GNDDO/1 parametrization.
C      0000 ... Use STO-3G scale factors.
C      1000 ... Use Slater's rules scale factors.
C     00000 ... Default (unit overlap matrix).
C     10000 ... Use the unit matrix for the overlap.
C     20000 ... Use the real overlap matrix.
C    100000 ... Do CNDO/2.
C    200000 ... Do INDO/2.
C    300000 ... Do ZINDO/1 (NYI).
C    400000 ... Do ZINDO/S.
C    500000 ... Do MINDO/3 (NYI).
C    600000 ... Do MNDO.
C    700000 ... Do AM1.
C    800000 ... Do PM3.
C    900000 ... Do PM3MM.
C   1000000 ... Do Harris functional through L511.
C   1100000 ... Do Harris functional scaling atomic densities for current
C               charge and multiplicity.
C   1200000 ... Harris XC but regular Coulomb iteration.
C   1300000 ... Harris (XC and atomic densities) through regular code.
C   1400000 ... Regular SCF with separate K, for testing.
C   1500000 ... J as usual but NDDO for K.
C   1600000 ... used internally as part of 15.
C   1700000 ... DFT-SCTB with tabulated parameters.
C   1800000 ... DFT-SCTB with analytic expressions.
C   1900000 ... EHT-SC.
C   2000000 ... Set 2e terms to zero.
C   2100000 ... Harris XC and DFTB-style charge iteration.
C   2200000 ... Harris XC and improved DFTB-style charge iteration.
C   2300000 ... PM6PFD with overlap.
C   2400000 ... PM6PFD with overlap and Harris XC.
C   2500000 ... PM6PFD with overlap and approximate XC.
C   2600000 ... NDDO with Mayer Bond Order correlation corrections.
C   27-38   ... Reserved for other methods with 2e integrals.
C   3900000 ... PM6
C   4000000 ... PMDDG
C   41      ... PM6E
C   42      ... PM7
C   43      ... PM6 with T transformed to OAO.
C   44      ... PM7TS.
C   45      ... PM7MOPAC.
C   46-98   ... Assumed to be ZDO methods.
C   9900000 ... External program.
C  100      ... Amoeba.
C  101-     ... Assumed to be MM methods.
C
C2IOp(42)
C
C   IOp(42) ... How to form NDDO core hamiltonian in L317:
C         0 ... Default (same as 1).
C         1 ... Read the integrals sequentially.
C         2 ... Load all the integrals into memory.
C
C2IOp(43)
C   IOp(43) ... Handling of background charge distribution.
C           00 ... same as 11;
C            1 ... consider external charges;
C            2 ... do not consider external charges;
C           10 ... consider self-consistent solvent charges;
C           20 ... do not consider self-consistent solvent charges;
C
C2IOp(44)
C
C     IOp(44) ... integral rejection using L318.
C           0 ... keep all integrals.
C           1 ... neglect four center transformed integrals.
C           2 ... neglect four center and 3 center (ab|ac) integrals.
C           3 ... neglect four center and three center (0,0) integrals.
C           4 ... NDDO approximation -- no (ab|xx) and no <a|X|b>
C           5 ... NDDO on 2e and V ints only -- T and S unchanged.
C           6 ... Do not transform 2e integrals, only 1e.
C
C2IOp(45)
C
C     IOp(45) ... transformation matrix in L318.
C           0 ... use S**-1/2.
C           1 ... just orthogonalize functions on the same center.
C           2 ... Use unit matrix (for debugging).
C
C     IOp(45) ... Order of multipoles in SCRF for L303.
C
C2IOp(46)
C
C     IOp(46) ... Whether to abort the job if badbas detects an error:
C           0 ... Default (yes).
C           1 ... yes.
C           2 ... no.
C
C2IOp(47)
C
C     IOp(47) ... Flags for use in PRISM and CalDFT throughout the program.
C          -2 ... Force use of only the MD paths for all calculations.
C          -1 ... Force use of only the OS path for all calculations.
C           Bit flags:
C           0 ... If bit 0 is set (use AllowP array) then read in
C                 a list of allowed paths.
C           1 ... Use expanded matrix logic for PBC exact exchange.
C           2 ... Reverse choice of whether to precompute distance matrix
C                 during numerical quadrature.
C           3 ... Skip consistency checks for XC quadrature
C           4 ... Do not do extra work to use cutoffs better, currently
C                 only affects CalDFT.
C           5 ... Reverse normal choice of diagonal/canonical sampling
C                 in Prism and PrmRaf.  The default is diagonal only
C                 on vector machines.
C           6 ... Trace input and output using Linda/subprocess.
C           7 ... Force single matrix code in CPKS.
C           8 ... Force all near field in FMM.
C           9 ... Turn off dynamic allocation of parallel work in
C                 CalDSu, CoulSu, and FMMEnt.
C          10 ... Force square loops, currently only in PrismC.
C          11 ... Turn off dynamic work allocation among Linda workers.
C          12 ... Reverse normal choice of Scat20 vs. replicated
C                 Fock matrices.  Default is to use replicated
C                 matrices only on Fujitsu and NEC.
C          13 ... Turn on Schwartz screening only in FoFCou, turning
C                 off heuristic screening.
C          14 ... Force separate evaluation of J and K terms.
C          15 ... Forbid use of gather/scatter digestion even for small
C                 numbers of density matrices.
C          16 ... Insist on gather/scatter digestion even for large
C                 numbers of matrices.  Does not affect FoFRaf, which
C                 only does inner loops over matrices.
C          17 ... Forbid use of Schwartz screening in FoFCou.
C          18 ... Don't compute on Linda master.
C          19 ... Do nuclear contribution in FoFCou even for non-PBC
C          20 ... Do not use special Coulomb algorithm in FoFCou.
C          21 ... Force dynamic parallel work logic even for single
C                 processor tasks.
C          22 ... Turn off use of Sqrt(P) in density-based cutoffs.
C          23 ... Use tabulated numerical values for atomic densities
C                 instead of Gaussian expansions.
C          24 ... Do allocation for parallel 2e integrals but run sequentially.
C          25 ... Do allocation for parallel XC but run sequentially.
C          26 ... Make all atoms large in XC quadrature.
C          27 ... Make all shells large in XC quadrature.
C          28 ... Do not symmetry reduce grid points on unique atoms.
C          29 ... Turn on use of precomputed XC weights.
C          30 ... Make Linda workers run sequentially.
C          31 ... Reserved for flag for calls to OneElI, etc. in parallel regions.
C          32 ... Use Euler-2 radial grid instead of Mura grid.
C          33 ... Set thresholds for CPU/GPU split to all GPU, for testing.
C
C2IOp(48)
C     IOp(48)   ... Options for FMM:
C    RRLLNNTTWW ... RR:  Range (default 2)
C                   LL:  LMax (default from tolerance)
C                   NN:  Number of levels (default 8)
C                   TT:  Tolerance (default 18)
C                   WW:  IWS (default 2).
C
C2IOp(49)
C     IOp(49)   ... More bitwise options for FMM and 2e integrals  The bits are:
C             0 ... Indicates whether FMM can be used by FoFCou.
C             1 ... Uncontract all shell pairs.
C             2 ... Apply symmetry to derivative distributions (NYI).
C             3 ... Do not save as many multipole expansions as possible in
C                   memory.
C             4 ... Turn on FMM print.
C             5 ... Convert to sparse storage under FoFCou for testing.
C             6 ... Split primitives for better boxification.
C             7 ... Default UseUAB/Use 256.
C             8 ... UseUAB, if 128 set.
C             9 ... Turn off parallelism in FMM (does not use parallel logic).
C            10 ... Set up for parallel FMM but run loops sequentially.
C            11 ... Do not default to FMM.
C            12 ... Force FMM on.
C            13 ... Set by PsmSet to indicate whether the NAtoms test for defaulting
C                   FMM was passed.
C            14 ... Turn on parallelism in FMM during CPHF.  Default is off because
C                   one job fails on one machine, which may still be a bug.
C            15 ... Force use of old box-box screening.
C            16 ... Do not Include 1/R or Erf(R)/R in box-box screening
C            17 ... Force use of non-cubic logic.
C            18 ... Turn off box-box screening.
C            19 ... Skip FF exchange.
C         20-22 ... Pure function control:
C                   0 ... Default, same as 1
C                   1 ... Convert densities, etc. to Cartesian.
C                   2 ... Transform 2e integrals to pure before digestion.
C                   3 ... Generate 2e integrals over real spherical
C                         harmonics.
C                   4 ... Generate 2e integrals over complex spherical
C                         harmonics.
C                   5 ... Generate 2e integrals over spinors.
C                   6 ... Generate 2e integrals over large and small
C                         components.
C
C2IOp(51)
C     IOp(51) ... Parameters for FMM box length (MMMMMNNNN):
C       MMMMM ... Box length when doing Coulomb will be MMMMM/1000 Bohr.
C                 The default is 2.5 Bohr.
C        NNNN ... Box length when doing Exchange will be NNNN/1000 Bohr.
C                 The default is 0.75 Bohr.
C     If doing both Coulomb and exchange at the same time, the maximum of
C     the two values is used.
C
C2IOp(52)
C     IOp(52) ... Turn off normal evaluation of ECP integrals.
C           0 ... Default:  if needed, ECP integrals are evaulated in L302.
C           1 ... Old routines will be used, so L302 does not do ECP ints.
C
C2IOp(53)
C     IOp(53) ... Accuracy in ECP integral evaluation:
C           0 ... Default.
C          -1 ... No Cutoffs
C           N ... 10**-N.
C
C2IOp(54)
C     IOp(54) ... Type of core density to use with ECPs:
C          -1 ... None
C           0 ... Default (None)
C           1 ... Non-relativistic
C           2 ... Relativistic
C
C2IOp(55)
C     IOp(55) ... Use of sparse storage:
C      N<-100 ... Yes, cutoff 5 x 10 ** (N+100)
C   -100<N<-4 ... Cutoff 10**(N+5) for testing new code.
C          -4 ... Reserved (used for nosparse in parsing).
C          -3 ... Yes, intermediate accuracy (10**-6)
C          -2 ... Yes, crude accuracy (10**-6)
C          -1 ... Yes, default accuracy (10**-8).
C           0 ... No
C           N ... Yes, cutoff 10**(-N)
C
C2IOp(56)
C     IOp(56) ... Cutoff for intermediate matrices during sparse
C                 operations:
C           0 ... 100 times smaller than storage cutoff.
C           N ... 10**(-N).
C
C2IOp(57)
C     IOp(57) ... No. of core electrons for Stuttgart/Dresden ECP's.
C
C2IOp(58)
C     IOp(58) ... Cholesky control options.
C
C2IOp(59)
C     IOp(59) ... Threshold for throwing away eigenvectors of S:
C       -NMLL ... M.LL*10**-N.
C           0 ... Default (10**-6)
C           N ... 10**-N.
C
C2IOp(60)
C     IOp(60) ... Control of orthogonalization and simplification of
C                 generalized contraction basis sets.
C          -1 ... Turn off orthogonalization and simplification
C           0 ... Default (2).
C           1 ... Orthogonalize and remove primitives with 0 coefficients
C                 (exact transformation).
C           2 ... Orthogonalize and remove primitives with 0 or small
C                 coefficients.
C           N ... Orthogonalize and remove primitives with coefficients
C                 less than 10**(-N)
C
C2IOp(61)
C     IOp(61) ... Sparse Semiempirical Hamiltonian Cutoffs in L302:
C          XX ... F(Mu,Lambda) atom--atom cutoff criterion (angstroms)
C                 Mu, Lambda are basis functions on different atoms.
C                 (defaults to 15 angstroms).
C        XX00 ... F(Mu,Nu) atom--atom cutoff criterion (angstroms)
C                 Mu, Nu are basis functions on the same atom.
C                 (defaults to no F(Mu,Nu) cutoff).
C
C2IOp(62)
C     IOp(62) ... Maximum allowed error in S over orthogonalized basis functions:
C          -2 ... Take value from chk file.
C          -1 ... Default.
C           0 ... Default (10**-9), unless overridden by saved values.
C           N ... 10**-(N).
C
C2IOp(63)
C     IOp(63) ... Debug option to test point charge FMM.
C           0 ... No.
C           1 ... Yes.
C           2 ... Yes, read parameters.
C          10 ... Also do forces.
C
C2IOp(64)
C     IOp(64) ... Set value for ILSW derivative flag.  Only active if IOp(39)=0.
C          -2 ... Set to zero
C          -1 ... Set to -1.
C           0 ... Leave alone.
C           N ... set to N.
C
C2IOp(65)
C     IOp(65) ... Number of k-points:
C          -1 ... Just Gamma point.
C           N ... About N points.
C          -N ... Old logic for NRecip=N.
C
CIOp(66)
C     IOp(66) ... Over-ride setting of NThInc in lineary dependence cutoff:
C          -2 ... Recover value from chk file.
C          -1 ... 0
C           0 ... 0 unless over-ridden by saved values.
C           N ... Scale value from IOp(59) by 2^N.
C
CIOp(67)
C     IOp(67) ... Electric-field dependent functions:
C           0 ... Default (on if already present in basis read from
C                 rwf or chk, otherwise off).
C           1 ... No.
C           2 ... Yes, with standard values.
C           3 ... Yes, with read-in values.
C
CIOp(70)
C     IOp(70) ... SCRF flag.
C           0 ... Default (1)
C           1 ... Use defaults.
C           2 ... Read setting from checkpoint.
C           3 ... Read setting from the input stream.
C           4 ... Read setting from checkpoint and modify them by reading
C                 from the input stream.
C           5 ... Read from rwf.
C        0100 ... Flag for macroiterations.
C        1000 ... SCI-PCM.
C        2000 ... D-PCM.
C        2100 ... C-PCM.
C        2200 ... IEF-PCM.
C        2300 ... IVC-PCM
C        3000 ... Not currently used.
C        4000 ... Onsager.
C       10000 ... Generate COSMOTHERMO output.
C       20000 ... Do COSMO style CPCM:  Klamt radii, iterative (implies g03defaults)
C       30000 ... Do SMD parametrization of non-electrostatic terms.
C      x00000 ... Flag for PCM family options:
C                 1 = include cavity-field effects;
C                 2 = setting for accurate DeltaG of solvation;
C                 3 = setting to reproduce G03 behaviour;
C     1000000 ... Flag to skip PCMInp as L124 already did it or we're doing
C                 flavor X of ONIOM-PCM;
C     2000000 ... Write the PCM charges on the checkpoint file.
C     3000000 ... Read the PCM charges from the checkpoint file.
C     4000000 ... Read and write the PCM charges from and to the checkpoint file.
C     5000000 ... Write the non-equilibrium PCM charges on the checkpoint file.
C     6000000 ... Read the non-equilibrium PCM charges from the checkpoint file.
C     7000000 ... Write the CC non-equilibrium PCM charges on the checkpoint file.
C     8000000 ... Read the CC non-equilibrium PCM charges from the checkpoint file.
C     9000000 ... Write the cLR non-equilibrium PCM charges on the checkpoint file.
C    00000000 ... Default, same as 30000000;
C    10000000 ... Do the PCM electrostatic cavity;
C    20000000 ... Do the PCM non-electrostatic cavity;
C    30000000 ... Do both the PCM electrostatic and non-electrostatic cavities;
C    40000000 ... Do neither the PCM electrostatic nor non-electrostatic cavities;
C
CIOp(71)
C     IOp(71) ... IDeriv level flag (for SCRF setup).
C
CIOp(72)
C     IOp(72) ... Solvent type flag (for SCRF setup).
C
C2IOp(74)
C     IOp(74) ... Type of exchange and correlation potentials:
C         -92 ... RSX-QIDH
C         -91 ... RSX-0DH
C         -90 ... wB97M(2) (nonvariational, only accessible via IOp)
C         -89 ... wB97M-V
C         -88 ... B97M-V
C         -87 ... wB97X-V
C         -86 ... wB97X-D3
C         -85 ... VV10 (Vydrov, Van Voorhis, JCP 133 244103 2010)
C         -84 ... revM11
C         -83 ... revM06   (Wang,Verma,Jin,Truhlar,He,PNAS,2018,115,10257)
C         -82 ... revM06-L (Wang,Jin,Yu,Truhlar,He,PNAS,2017,114,8487)
C         -81 ... M11plus (revised M11 + R35 correlation).
C         -80 ... GAMxc (Hu, Truhlar).
C         -79 ... PBE-QIDH.
C         -78 ... PBE0-DH.
C         -77 ... DSD-PBEP86 (double hybrid, DFT-D3).
C         -76 ... PW6B95-D3.
C         -75 ... PW6B95.
C         -74 ... M08-HX.
C         -73 ... MN15.
C         -72 ... MN15-L.
C         -71 ... LC-wHPBE.
C         -70 ... MN12-SX.
C         -69 ... N12-SX.
C         -68 ... MN12-L.
C         -67 ... N12.
C         -66 ... M11L.
C         -65 ... SOGGA11X.
C         -64 ... M11.
C         -63 ... SOGGA11.
C         -62 ... HISSaPBE.
C         -61 ... HISSbPBE.
C         -60 ... B2PLYP-D3 (double hybrid, DFT-D3).
C         -59 ... B97-D (DFT-D3).
C         -58 ... wB97X-D
C         -57 ... wB97X
C         -56 ... wB97
C         -55 ... M06-2X
C         -54 ... M06
C         -53 ... M06-L
C         -52 ... M06-HF
C         -51 ... HSEH1PBE.
C         -50 ... mPW2PLYP-D (double hybrid).
C         -49 ... B2PLYP-D (double hybrid).
C         -48 ... mPW2PLYP (double hybrid).
C         -47 ... B2PLYP (double hybrid).
C         -46 ... PAPF-D.
C         -45 ... PAPF.
C         -44 ... APF-D.
C         -43 ... APF.
C         -42 ... B97-D.
C         -41 ... LC-wPBE.
C         -40 ... CAM-B3LYP.
C         -39 ... OAPF.
C         -38 ... M052X.
C         -37 ... M05.
C         -36 ... HSE1PBE.
C         -35 ... TPSSh
C         -34 ... BMK.
C         -33 ... X3LYP
C         -32 ... t-HCTH hybrid
C         -31 ... t-HCTH
C         -30 ... OmPW3PBE.
C         -29 ... OmPW1PBE.
C         -28 ... OmPW1LYP.
C         -27 ... OmPW1PW91.
C         -26 ... PBEH1PBE.
C         -25 ... HSE2PBE.
C         -24 ... O3LYP.
C         -23 ... HCTH407.
C         -22 ... HCTH147.
C         -21 ... B97-2.
C         -20 ... B97-1.
C         -19 ... HCTH93.
C         -18 ... B98.
C         -17 ... B1B95.
C         -16 ... BA3PBE.
C         -15 ... BA1PBE.
C         -14 ... PBE3PBE.
C         -13 ... PBE1PBE
C         -12 ... mPW3PBE.
C         -11 ... mPW1PBE.
C         -10 ... mPW1LYP.
C          -9 ... LG1LYP.
C          -8 ... B1LYP.
C          -7 ... mPW91PW91.
C          -6 ... Becke3 with Perdew 91 correlation.
C          -5 ... Becke3 using VWN/LYP for correlation.
C          -4 ... Becke 3 with Perdew 86 correlation.
C          -3 ... Becke "Half and Half" with LYP/VWN correlation.
C          -2 ... Becke "Half and Half":  0.5 HF + 0.5 LSD
C          -1 ... Do only coulomb part; skip exchange-correlation.
C          00 ... Default, same as 100.
C          01 ... Vosko-Wilk-Nusair method 5 correlation.
C          02 ... Lee-Yang-Parr correlation.
C          03 ... Perdew 81 correlation.
C          04 ... Perdew 81 + Perdew 86 correlation.
C          05 ... VWN 80 (LSD) correlation
C          06 ... VWN 80 (LSD) + Perdew 86 correlation
C          07 ... [unused]
C          08 ... PW91
C          09 ... PBE
C          10 ... VSXC
C          11 ... Bc96
C          18 ... VWN5+P86
C          19 ... LYP+VWN5 for scaling
C          20 ... KCIS correlation
C          21 ... Becke-Roussel correlation (NYI)
C          22 ... PKZB correlation
C          23 ... TPSSc
C          24 ... t-HCTH (JCP 116, 9559 (2002))
C          25 ... t-HCTH hybrid (JCP 116, 9559 (2002))
C          26 ... BMK (Boese and Martin, JCP 121, 3405 (2004))
C          27 ... M05 (Zhao,Schultz,Truhlar, JCP 123 (2005) 161103)
C          28 ... M05-2X (Zhao,Schultz,Truhlar, JCTC 2006 in press)
C          29 ... OAPF (Austin, Petersson, Frisch, ...)
C          30 ... B97-D (Grimme, JCC 2006, 27, 1787)
C          31 ... APF (Austin, Petersson, Frisch, ...)
C          32 ... PAPF (Austin, Petersson, Frisch, ...)
C          33 ... M06-HF (Zhao,Truhlar, JPC A 2006, 110, 13126)
C          34 ... M06-L (Zhao,Truhlar, JCP 2006, 125, 194101)
C          35 ... M06 (Zhao,Truhlar, Theo Chem Acc 2008, 120, 215)
C          36 ... M06-2X (Zhao,Truhlar, Theo Chem Acc 2008, 120, 215)
C          37 ... wB97 (J.-D. Chai, M. Head-Gordon, JCP 128, 084106 (2008))
C          38 ... wB97X (J.-D. Chai, M. Head-Gordon, JCP 128, 084106 (2008))
C          39 ... wB97X-D (J.-D. Chai, M. Head-Gordon, PCCP 10, 6615 (2008))
C          40 ... revTPSSc
C          41 ... SOGGA11 (Peverati, Zhao, Truhlar, JPCL 2, 1991 (2011))
C          42 ... M11 (Peverati, Truhlar, JPCL 2, 2810 (2011))
C          43 ... SOGGA11-X (Peverati, Truhlar, JCP 135, 191102 (2011))
C          44 ... M11-L (Peverati, Truhlar, JPCL 3, 117 (2012))
C          45 ... N12 (Peverati, Truhlar, JCTC 8, 2310 (2012))
C          46 ... MN12-L (Peverati, Truhlar, PCCP DOI: 10.1030/c2cp42025b)
C          47 ... N12-SX (Peverati, Truhlar, PCCP submitted)
C          48 ... MN12-SX (Peverati, Truhlar, PCCP submitted)
C          49 ... CVDFT correlation
C          50 ... CCDFT correlation
C          51 ... MN15-L correlation.
C          52 ... MN15 correlation.
C          53 ... M08-HX correlation.
C          54 ... BC95 correlation in PW6B96.
C          55 ... GAMc correlation (Yu, Truhlar).
C          56 ... M11 correlation with parametrized R35 term (using XCParU).
C          57 ... LYP+R35 correlation (using XCParU).
C          58 ... [reserved for more development work using XCParU].
C          59 ... [reserved for more development work using XCParU].
C          60 ... M11plus correlation.
C          61 ... revM06-L correlation.
C          62 ... revM06 correlation.
C          63 ... revM11 correlation.
C        1000 ... Hartree-Fock exchange.
C        2000 ... Hartree-Fock-Slater exchange (Alpha = 2/3)
C        3000 ... X-alpha exchange (alpha= 0.7)
C        4000 ... Becke 1988 exchange
C        5000 ... LG exchange (deprecated)
C        6000 ... PW91 exchange
C        7000 ... Gill 96 exchange (deprecated)
C        8000 ... PW86 exchange
C        9000 ... mPW exchange
C       10000 ... PBE exchange
C       11000 ... [reserved to map 300]
C       12000 ... VSXC exchange
C       13000 ... revPW86 exchange
C       14000 ... B98 (JCP 108,9624(1998) eq.2c ) exchange
C       15000 ... HCTH (JCP 109,6264 (1998) exchange
C       16000 ... B97-1 (CPL 316,160(2000)) exchange
C       17000 ... B97-2 (JCP 115,9233(2001)) exchange
C       18000 ... HCTH147 exchange
C       19000 ... HCTH407 exchange
C       20000 ... OPTX exchange
C       21000 ... OPTX exchange as in O3LYP
C       22000 ... XVa exchange (NYI)
C       23000 ... Becke-Roussel '88 exchange
C       24000 ... PKZB exchange
C       25000 ... TPSSX exchange
C       26000 ... HSE03 (JCP 118,8207(2003)) exchange
C       27000 ... PBEHole (JCP 109,3313(1998)) exchange
C       28000 ... Old mPW exchange (local scaling in non-local term)
C       29000 ... t-HCTH (JCP 116, 9559 (2002))
C       30000 ... t-HCTH hybrid (JCP 116, 9559 (2002))
C       31000 ... X (0.765*B88+0.235*PW91) (PNAS 101(2004) 2673)
C       32000 ... BMK (Boese and Martin, JCP 121, 3405 (2004))
C       33000 ... M05 (Zhao,Schultz,Truhlar, JCP 123 (2005) 161103)
C       34000 ... M05-2X (Zhao,Schultz,Truhlar, JCTC 2006 in press)
C       35000 ... OAPF (Austin, Petersson, Frisch, ...)
C       36000 ... B97-D (Grimme, JCC 2006, 27, 1787)
C       37000 ... APF (Austin, Petersson, Frisch, ...)
C       38000 ... PAPF (Austin, Petersson, Frisch, ...)
C       39000 ... HSE + Henderson
C       40000 ... M06-HF (Zhao,Truhlar, JPC A 2006, 110, 13126)
C       41000 ... M06-L (Zhao,Truhlar, JCP 2006, 125, 194101)
C       42000 ... M06 (Zhao,Truhlar, Theo Chem Acc 2008, 120, 215)
C       43000 ... M06-2X (Zhao,Truhlar, Theo Chem Acc 2008, 120, 215)
C       44000 ... wB97 (J.-D. Chai, M. Head-Gordon, JCP 128, 084106 (2008))
C       45000 ... wB97X (J.-D. Chai, M. Head-Gordon, JCP 128, 084106 (2008))
C       46000 ... wB97X-D (J.-D. Chai, M. Head-Gordon, PCCP 10, 6615 (2008))
C       47000 ... HISS (Henderson,Izmaylov,Scuseria,Savin, JCP 127, 22103 (2007))
C       48000 ... revTPSSX
C       49000 ... SOGGA11 (Peverati, Zhao, Truhlar, JPCL 2, 1991 (2011))
C       50000 ... M11 (Peverati, Truhlar, JPCL 2, 2810 (2011))
C       51000 ... SOGGA11-X (Peverati, Truhlar, JCP 135, 191102 (2011))
C       52000 ... M11-L (Peverati, Truhlar, JPCL 3, 117 (2012))
C       53000 ... N12 (Peverati, Truhlar, JCTC 8, 2310 (2012))
C       54000 ... MN12-L (Peverati, Truhlar, PCCP DOI: 10.1030/c2cp42025b)
C       55000 ... N12-SX (Peverati, Truhlar, PCCP submitted)
C       56000 ... MN12-SX (Peverati, Truhlar, PCCP submitted)
C       57000 ... [reserved to produce B values for XDM]
C       58000 ... [reserved to run HF + XDM]
C       59000 ... PiLDA (R3.5 functional).
C       60000 ... LDA via R3.5 logic for testing.
C       61000 ... PiLOS (R3.5 functional).
C       62000 ... PiPBE (R3.5 functional).
C       63000 ... PiPBELOS (R3.5 functional).
C       64000 ... LDA similarity metric based on R3.5 logic.
C       65000 ... wPiLDA (R3.5 functional).
C       66000 ... M11plus exchange.
C       67000 ... revM06-L exchange.
C       68000 ... revM06 exchange.
C       69000 ... revM11 exchange.
C       70000 ... CVDFT exchange.
C       71000 ... MN15-L exchange.
C       72000 ... MN15 exchange.
C       73000 ... M08-HX exchange.
C       74000 ... PW6 exchange in PW6B96.
C       75000 ... GAMx exchange (Yu, Truhlar).
C       76000 ... M11 exchange with parametrized R3.5 term (using XCParU).
C       77000 ... [reserved for more development work using XCParU].
C
C     So 1000 is Hartree-Fock, 2000 is Hartree-Fock-Slater,
C     2005 is Local Spin Density, and 4002 is BLYP.
C
C     1xxxxxx ... Do Hirao's long range correction (JCP 115(2001) 3540).
C     2xxxxxx ... Do Harris XC with full J.
C     3xxxxxx ... Do Harris with the specified functional.
C     4xxxxxx ... Do Harris XC with DFTB-style J.
C     5xxxxxx ... Do Harris XC with improved DFTB-style J.
C
C2IOp(75)
C     IOp(75) ... Number of radial and angular points in numerical
C                 integration for DFT:
C           0 ... Default (-5).
C           1 ... SG1 pruned grid.
C           2 ... Even sleazier grid than SG1 used for CPHF.
C           3 ... Pruned (75,194) which is not good for much.
C           4 ... FineGrid.
C          -4 ... FineGrid unless uncontracting, then 199302.
C           5 ... UltraFine.
C          -5 ... UltraFine unless uncontracting, then 199590.
C           7 ... SuperFine.
C          -7 ... SuperFine unless uncontracting, then 299974.
C      IIIJJJ ... III radial points, JJJ angular points.
C     -IIIJJJ ... III radial points, and a spherical product angular
C                 grid with JJJ theta points and 2*JJJ phi points.
C
C2IOp(76)
C     IOp(76) ... Mixing of HF and DFT.  Negative values correspond
C                 to standard combinations of HF exchange, local and
C                 non-local exchange, and local and non-local
C                 correlation.
C         -36 ... PBE-QIDH coefficients.
C         -35 ... PBE0-DH coefficients.
C         -34 ... DSD-PBEP86 coefficients.
C         -33 ... PW6B95 and PW6B95-D3 coefficients.
C         -32 ... M08-HX coefficients.
C         -31 ... MN15 coefficients.
C         -30 ... SOGGA11-X coefficients.
C         -29 ... HSEH1, N12-SX and MN12-SX coefficients.
C         -28 ... M06-2X coefficients.
C         -27 ... M06, wB97, wB97X, wB97X-D, HISS-B, HISS-A, M11,
C                 M11plus and LC-wHPBE coefficients.
C         -26 ... M06-HF coefficients.
C         -25 ... mPW2PLYP coefficients.
C         -24 ... B2PLYP coefficients.
C         -23 ... APF coefficients.
C         -22 ... Unused.
C         -21 ... LC-wPBE coefficients.
C         -20 ... CAM-B3LYP coefficients.
C         -19 ... OAPF coefficients.
C         -18 ... M05-2X coefficients.
C         -17 ... TPSSh coefficients.
C         -16 ... BMK coefficients.
C         -15 ... X3LYP coefficients.
C         -14 ... tHCTH coefficients.
C         -13 ... B1B95/M05 coefficients.
C         -12 ... HSE1PBE, HSE2PBE coefficients.
C         -11 ... Unused
C         -10 ... O3LYP coefficients.
C          -9 ... B97-2 coefficients.
C          -8 ... B97-1 coefficients.
C          -7 ... HCTH coefficints.
C          -6 ... B98 coefficients.
C          -5 ... mPW91PW91 coefficients.
C          -4 ... Becke 3 coefficients:  aLSD + (1-a)HF + b(dBx) +
C                 VWN + c(LYP-VWN), with a=0.8 b=0.72 c=0.81
C                 Note that Becke actually used Perdew correlation
C                 rather than LYP.
C          -3 ... Becke "Half and Half" 0.5 HF + 0.5 Xc + Corr
C          -2 ... Coefficients of 0 and 0 (no exchange).
C          -1 ... Coefficients of 0.0 and 1.0 for DFT and HF,
C                 respectively.
C           0 ... Default: pure HF, DFT or mixed in accord with IOp(76)
C  MMMMMNNNNN ... Mixture of MMMMM/10000 DFT exchange and NNNNN/10000
C                 HF exchange.  The DFT exchange factor multiplies any
C                 implied by IOp(74) or set by IOp(77).
C
C2IOp(77)
C     IOp(77) ... Mixing of local and non-local exchange:
C          -1 ... 0 for both.
C           0 ... Default (coefficients of 1 and zero or as determined
C                 by IOp(76)
C  MMMMMNNNNN ... MMMMM/10000 non-local plus NNNNN/10000 local.  Sign is
C                 applied to the local term.
C
C                 For the HSE03 functional, these coefficients scale the
C                 short range (MMMMM) and long range (NNNNN) terms.
C
C2IOp(78)
C     IOp(78) ... Mixing of local and non-local correlation:
C          -1 ... 0 for both.
C           0 ... Default (coefficients of 1 and zero as determined
C                 by IOp(76)
C  MMMMMNNNNN ... MMMMM/10000 non-local plus NNNNN/10000 local.  Sign is
C                 applied to the local term.
C
C     In L510, 1 to set up for CAS-MP2 or 2 to do spin-orbit calculation.
C
C2IOp(79)
C IOp(79) ... Range cutoff in Becke weights.
C       0 ... Default (SS weights)
C      -1 ... Use SS weights.
C      -2 ... Use Becke weights with default cutoff of 30 au.
C      -3 ... Use Savin weights.
C   -M<-3 ... Use SS weights with XCal = M/1000.
C       N ... Use Becke weights with cutoff N Bohr.
C
C2IOp(80)
C     IOp(80) ... Range for microbatching in DFT.  Negative to turn off
C                 screening of basis functions and grid points.
C                 1000000000 turns of microbatching logic.
C
C2IOp(81)
C     IOp(81) ... 6d integrated terms to include in functional:
C           0 ... Default (none, or as prescribed by IOp(74)).
C           1 ... VV10 dispersion (JCP 133 244103 2010)
C           2 ... VV10 dispersion as in wB97X-V
C           3 ... VV10 dispersion as in wB97M(2)
C           4 ... rVV10 dispersion (PRB 87 041108 2013) (NYI)
C         999 ... Turn off N^6 terms which would be on by default.
C        N000 ... Set screening algorithm to N.  Default is 2, 9 meands
C                 no screening
C       M0000 ... Set screening threshold to 10^-14+M
C      L00000 ... How to handle N^6 terms (default 3)
C                 1 ... In single CalDSu call with regular XC.
C                 2 ... In separate CalDSu call using mixed grids.
C                 3 ... In separate CalDSu call using the reduced
C                       grid in both integrations.
C                 4 ... Use NumIn6.
C                 5 ... Use NumIn6 with old algorthm.
C     K000000 ... Use grid K steps down from global default for N^6 terms.
C                 Default is 9, 9=> same as general grid.
C                 (Defaults are handled by DecN6F).
C    J0000000 ... Bit 0:  add the silly constant to the total energy.
C
C2IOp(82)
C     IOp(82) ... Fitting density basis set for Coulomb in DFT.
C          -1 ... None.
C           0 ... Default (-1).
C           N ... Same numbering of basis sets as for AO basis, including
C                 7=General basis.  See comments for IOp(5) and IOp(6)
C                 28=Generate automatically from AO basis.
C
C2IOp(83)
C     IOp(83) ... Equivalent of IOp(6) for density basis.
C                 For auto-generated basis sets:
C          MN ... -1 keep all generated functions.  Otherwise, an AO
C                 shell with angular momentum LAO generates a DBF shell
C                 with angular momenta 0 up to LDB, where if LVal is the
C                 highest valence (occupied) LAO then if LAO<=LVal,
C                 LDB=2*LAO, while if LAO>LVal LDB = LAO + Max(LVal,1) + M
C                 If N>0 then LDB is limited to N-1, i.e., all angular
C                 momenta of N or higher are discarded.
C
C2IOp(84)
C     IOp(84) ... Equivalent of IOp(7) for density basis.
C                 For auto-generated basis sets:
C           0 ... Default (4022)
C           1 ... Use all products of AOs.
C           2 ... Use only AO primitives squared in fititing basis.
C          10 ... Do not split shells
C          20 ... Split F and higher shells away from S=P=D.
C         N00 ... Use 1.5 + N/4 as the test for similar exponents during
C                 auto-generation of fitting sets.
C        1000 ... Use old (G03) algorithm.
C        2000 ... Use new algorithm.
C        3000 ... Use algorithm 3.
C        4000 ... New iterative merging of shells, monotonic L.
C
C2IOp(85)
C     IOp(85) ... Pure vs. Cartesian functions in density basis.
C           0 ... Default (pure).
C           1 ... Pure.
C           2 ... Cartesian.
C
C2IOp(86)
C     IOp(86) ... Discard basis functions based on angular momentum:
C           0 ... No.
C           N ... Discard basis functions with angular momentum >= N.
C
C2IOp(87)
C     IOp(87) ... Discard density basis functions based on angular momentum:
C           0 ... No.
C           N ... Discard density basis functions with angular momentum >= N.
C
C2IOp(88)
C     IOp(88) ... Modification of internally stored density and NEO basis sets:
C           0 ... None.
C           1 ... Read in general basis data in addition to setting up a
C                 standard density basis.
C          10 ... Massage the data in Common /B/ and Common /Mol/.
C         100 ... Add ghost atoms to /B/ so that every shell is on a separate
C                 center.  This is also done if requested in IOp(10).
C        1000 ... Split S=P density basis shells into separate S and P shells.
C        2000 ... Do not split S=P density shells.
C       10000 ... Split S=P=D=... density shells into S=P, D, F, ...
C       20000 ... Do not split density S=P=D... shells.
C      100000 ... Read in general basis data in addition to setting up a
C                 standard NEO basis.
C
C2IOp(89)
C     IOp(89) ... Set up for density fitting.
C           0 ... Default (102 if a fitting set has been included and
C                 pure DFT is being used, 1 otherwise).
C           1 ... Do not use density fits.
C           2 ... Use fits, forming Z = modified A^-1.
C           3 ... Use fits, solving iterative with stored A.
C           4 ... Use fits, solving iterative with direct products,
C                 with A formed to generate preconditioning.
C           5 ... Iterative, no formation of A.
C           6 ... Form A' over neutral distributions via multiplies by A.
C           7 ... Form A' over neutral distributions via direct products.
C         1xx ... Form inverse matrix once.
C         2xx ... Solve iteratively with no preconditioning
C         3xx ... Solve iteratively with diagonal preconditioning.
C         4xx ... Solve iteratively with symmetric block-diagonal
C                 preconditioning.
C         5xx ... Solve iteratively with non-symmetric block-diagonal
C                 preconditioning.
C         6xx ... Solve non-iterative using precomputed A'^-1.
C       1xxxx ... Put all functions into a single block in forming
C                 the preconditioning matrix.
C      1xxxxx ... Form the full preconditioning matrix (not block-diagonal)
C     0xxxxxx ... Default, same as 1xxxxxx
C     1xxxxxx ... Don't set up fitting if exact exchange is in use.
C     2xxxxxx ... Set up fitting regardless and do one fit with the
C                 converged SCF density.
C     3xxxxxx ... Set up fitting regardless and use for Coulomb during
C                 iterations even if exact exchange is used (NYI).
C    10000000 ... Fit using Coulomb operator (default).
C    20000000 ... Fit using overlaps.
C
C2IOp(90)
C     IOp(90) ... Thresholds for density fitting
C        MMNN ... 10**(-MM) on iterative solution, default MM=09.
C                 10**(-NN) on generalized inverse, default NN=06.
C
C2IOp(91)
C     IOp(91) ... Scalar relativistic core Hamiltonian:
C           0 ... Default (1)
C           1 ... Non-relativistic.
C           2 ... RESC.
C           3 ... Douglass-Kroll-Hess 0th order.
C           4 ... Douglass-Kroll-Hess 2nd order.
C           5 ... DKH 4th order, including SO terms.
C           6 ... Complete X2C, including SO terms.
C          00 ... Default (10).
C          10 ... Do Boettinger scaling of 1e SO to approximate effect
C                 of 2e terms.
C          20 ... Do not rescale SO terms.
C         100 ... Multiply SO terms by 10 for debugging.
C         N00 ... Multiply SO terms by 10 * 10^(N-1) for debugging.
C        1000 ... Multiply SO terms by half.
C        2000 ... Multiply SO terms by two.
C        3000 ... Multiply SO terms by -two.
C        4000 ... Multiply SO terms by zero.
C
C2IOp(92)
C     IOp(92) ... Whether read-in basis sets are in terms of normalized
C                 primitives?
C           0 ... Default (3232).
C           1 ... AO coefficients are for raw primitives;
C           2 ... AOs have overlap normalization;
C           3 ... AOs have Coulomb normalization.
C          10 ... DBF coefficients are for raw primitives;
C          20 ... DBFs have overlap normalization;
C          30 ... DBFs have Coulomb normalization;
C         100 ... Do not normalize AOs contraction coefficients.
C         200 ... Use overlap normalization for AOs contraction coefficients.
C         300 ... Use Coulomb normalization for AOs contraction coefficients.
C        1000 ... Do not normalize DBFs contraction coefficients.
C        2000 ... Use overlap normalization for DBFs contraction coefficients.
C        3000 ... Use Coulomb normalization for DBFs contraction coefficients.
C
C2IOp(93)
C     IOp(93) ... Nuclear charge distribution:
C           0 ... Default (1, unless scalar relativistic)
C           1 ... Point nuclei
C           2 ... Single s Gaussians using formula of Quiney et. al
C           3 ... Very tight single s Gaussians, for debugging.
C           4 ... Same as 2 but exponents are 100x smaller, for debugging.
C           5 ... Filatov-Cremer radii
C         10x ... Include nuclear charge distributions in DBF set.
C        Mxxx ... Use method M to handle nuclear charges during density fitting.
C       00000 ... Default (1).
C       10000 ... Use nuclear density in core Hamiltonian if present.
C       20000 ... Do not use nuclear density in core Hamiltonian even if present.
C      100000 ... Load array of nuclear radii with the values used.
C
C2IOp(94)
C     IOp(94) ... Range of PBC cells in Bohr.
C           0 ... default (100).
C           N ... N Bohr.
C          -M ... Multiply usual range by M.
C
C2IOp(95)
C     IOp(95) ... Minimum number of PBC cells.
C          -N ... At least N cells in each direction.
C           0 ... Based on range estimate (IOp(94)).
C           N ... At least N cells total.
C
C2IOp(96)
C     IOp(96) ... Number of PBC cells for DFT:
C           0 ... As many as look significant.
C           N ... At least N.
C
C2IOp(97)
C     IOp(97) ... Number of PBC cells for exact exchange:
C           0 ... As many as look significant.
C           N ... At least N.
C
C2IOp(98)
C     IOp(98) ... Maximum number of density matrices in PBC.
C           0 ... Default, based on number of cells having overlap with cell 0.
C           N ... No more than N matrices.
C
C2IOp(99)
C     IOp(99) ... Whether to set up precomputed quadrature grid in L302:
C           0 ... Default (4 if doing DFT, -1 otherwise).
C          -1 ... No
C           1 ... Yes, storing only grid parameters.
C           2 ... Yes, storing grid parameters and weights.
C           3 ... Yes, storing grid parameters, weights, and point coordinates.
C           4 ... Yes, storing only dimensions.
C
C2IOp(100)
C     IOp(100) ... Minimum Number of PBC cells for PBC-MP2
C            0 ... Same as for HF exchange.
C            N ... N.
C
C2IOp(101)
C     IOp(101) ... Maximum range of cells
C           -N ... No more than N in each direction
C            0 ... No limit.
C            N ... No more than N total.
C
C2IOp(102)
C     IOp(102) ... Number of density fittings solutions to save
C                  from previous SCF iterations.  Default is 6
C                  (using 5 previous solutions plus the current
C                  right-hand side to generate the initial guess).
C                  Negative to use projected equations rather than
C                  least-squares.
C
C2IOp(103)
C     IOp(103) ... Maximum number of vectors allowed in expansion
C                  space during iterative density fitting.  Default
C                  is Max(NDBF/2,1000).
C
C2IOp(104)
C     IOp(104) ... Maximum number of iterations during iterative
C                  density fitting.  Default is Max(1000,NDBF+100).
C
C2IOp(105)
C     IOp(105) ... Re-use of PBC cell data.
C            0 ... Default (re-use if present).
C            1 ... Reuse.
C            2 ... Do not reuse.
C            3 ... Read from chk file.
C
C2IOp(106)
C     IOp(106) ... Over-ride default number of atoms threshold for
C                  turning on FMM (for debugging).  This number is
C                  scaled up appropriately if symmetry is in use,
C                  to compensate for the loss of some symmetry with FMM.
C            0 ... Default (60)
C            N ... N atoms for the C1 case.
C
C2IOp(107)
C     IOp(107) ... Omega for short/long range Hartree-Fock exchange
C            0 ... Standard HF exchange
C   MMMMMNNNNN ... Short range HF exchange with NNNNN/10000 and long range
C                  exchange with MMMMM/10000.
C
C2IOp(108)
C     IOp(108) ... Omega for short/long range DFT exchange
C            0 ... Standard DFT exchange or default from functional
C   MMMMMNNNNN ... Short range DFT exchange with NNNNN/10000 and long range
C                  DFT exchange with MMMMM/10000.
C
C2IOp(109)
C     IOp(109) ... Omega for short/long range DFT correlation
C            0 ... Standard DFT correlation or default from functional
C   MMMMMNNNNN ... Short range DFT correlation with NNNNN/10000 and long range
C                  DFT correlation with MMMMM/10000.
C
C2IOp(110)
C     IOp(110) ... Threshold in precomputed XC quadrature grid,
C                  over-riding default or value in IOp(27).
C            0 ... As implied by IOp(27).
C            N ... 10^-N.
C
C2IOp(111)
C     IOp(111) ... Extra PBC printing:
C            1 ... Print table of cells.
C
C2IOp(112)
C     IOp(112) ... Huckel/DFTB parameters:
C            0 ... Default (13)
C            3 ... Hoffman parameters.
C            4 ... Pykko parameters.
C            5 ... Huckel initial guess parameters.
C           00 ... Default (10 for Huckel, 20 for DFTB).
C           10 ... Use standard parameters.
C           20 ... Read parameters to over-ride the standard ones.
C           30 ... Read parameters from rwf file 738.  Implies basis
C                  set also read from the rwf.
C           40 ... Read parameters from chk file 738.  Implies basis
C                  set also read from the chk.
C       1000xx ... Old algorithm (read parameters in L302).
C
C2IOp(113)
C     IOp(113) ... Generate SABF data
C           00 ... Default (121).
C            1 ... Generate AO basis function SABF data if symmetry is on.
C            2 ... Make AO SABF data C1 regardless.
C           10 ... Generate density basis function SABF data if symmetry is on.
C           20 ... Make density basis SABF data C1 regardless.
C          100 ... Generate NEO basis function SABF data if symmetry is on.
C          200 ... Make NEO SABF data C1 regardless.
C
C2IOp(114)
C     IOp(114) ... Factor for number of significant basis functions
C                  allocation in XC quadrature allocation.
C            0 ... Default:  use amount computed by LdMGrd.
C            N ... Scale values by N/10.
C
C2IOp(115)
C     IOp(115) ... Factor for number of significant atoms allocation in
C                  XC quadrature allocation.
C            0 ... Default:  use amount computed by LdMGrd.
C            N ... Scale values by N/10.
C
C2IOp(116)
C     IOp(116) ... Type of SCF:
C           -2 ... Take from the chk file.
C           -1 ... Ignore ILSW and determine on the fly.
C            0 ... Take from ILSW.
C            1 ... Real RHF.
C            2 ... Real UHF.
C            3 ... Complex RHF
C            4 ... Complex UHF
C            5 ... Complex, but use ILSW to decide whether RHF/UHF.
C            7 ... GHF using real basis functions
C           11 ... Complex RHF, complex spherical harmonic basis.
C           12 ... Complex UHF, complex spherical harmonic basis.
C           15 ... GHF, complex spinorbital basis (NYI).
C           19 ... GHF, spinor basis (NYI).
C           23 ... DF, spinor basis (NYI).
C          101 ... Real ROHF.
C          201 ... Unrestricted if derivatives are being done but RO
C                  single points; used for RO-compound methods.
C
C2IOp(117)
C     IOp(117) ... Handling Spin-orbit ECPs.
C            0 ... Default; include them if present and doing GHF.
C            1 ... Always compute SO terms.
C            2 ... Never compute SO terms.
C
C2IOp(118)
C     IOp(118) ... Extra memory for integral evaluation.
C            0 ... None.
C            N ... Add N words to the estimated memory requirements for
C                  direct integral evaluation, in all links.
C           -N ... Always use N words regardless of basis or order of
C                  derivatives.
C
C2IOp(119)
C     IOp(119) ... Coefficients of short/long range Hartree-Fock exchange.
C            0 ... Standard HF exchange.
C   MMMMMNNNNN ... NNNNN/10000 short range and MMMMM/10000 long range exchange.
C                  The signs can be changed by IOp(130) (see below).
C
C2IOp(120)
C     IOp(120) ... Coefficients of short/long range DFT exchange.
C            0 ... Standard DFT exchange or default from functional.
C   MMMMMNNNNN ... NNNNN/10000 short range and MMMMM/10000 long range.
C                  The signs can be changed by IOp(131) (see below).
C
C2IOp(121)
C     IOp(121) ... Coefficients of short/long range DFT correlation.
C            0 ... Standard DFT correlation or default from functional.
C   MMMMMNNNNN ... NNNNN/10000 short range and MMMMM/10000 long range.
C                  The signs can be changed by IOp(132) (see below).
C
C2IOp(122)
C     IOp(122) ... QM/QM embedding (L302):
C            2 ... Evaluate QM:QM full density embedding contribution.
C
C2IOp(123)
C     IOp(123) ... Phase convention for complex orbitals:
C            0 ... Normal; largest coefficient set to 1.
C            1 ... Largest coefficient set to i in each orbital.
C            2 ... Largest coefficient set to i in first orbital, i^2
C                  in second, etc.
C            3 ... Largest coefficient set to phase 60 degrees.
C            4 ... Largest coefficient set to phase 60 degrees, then
C                  120, etc.
C
C2IOp(124)
C     IOp(124) ... Empirical dispersion term.
C            0 ... Default (same as 2)
C            1 ... Add it regardless;
C            2 ... Add it for the DFT functionals for which it has been
C                  defined and parametrized and for which a specific name
C                  has been defined in Link1.
C            3 ... Add it for the DFT functionals for which it has been
C                  defined and parametrized;
C            4 ... Do not add it regardless;
C           10 ... Force dispersion type 1 (APF-D).
C           20 ... Force dispersion type 2 (Grimme B97-D).
C           30 ... Force dispersion type 3 (Grimme DFT-D3).
C           40 ... Force dispersion type 4 (Grimme DFT-D3(BJ)).
C           50 ... Force dispersion type 5 (Grimme D3, PM7 version).
C          000 ... Whether to change Grimme dispersion based on functional.
C                  Defaulted based on lowest digit.
C          100 ... Do the change.
C          200 ... Do not do the change.
C        NNxxx ... Use Grimme parameters for hybrid functional NN
C                  (see IOp(74)).
C      MMMMxxx ... Use Grimme parameters for pure functional MMMM
C                  (see IOp(74)).
C     10000000 ... Kill the job when atomic parameters are unavailable.
C     20000000 ... Continue the calculation even if some of the atomic
C                  parameters are unavailable.
C
C2IOp(125)
C     IOp(125) ... Scaling of AA/BB and AB components of E(2).
C           -3 ... 0 for AB
C           -2 ... 0 for AA/BB
C           -1 ... 0 for both.
C            0 ... Default (1 for both).
C   MMMMMNNNNN ... MMMMM/10000 for AA/BB, NNNNN/10000 for AB.
C                  A value over 50000 for either translates to -(NNNNN-50000).
C
C2IOp(126)
C     IOp(126) ... Omega for short/long range 1/r operator in E(2,AA)
C                  and E(2,BB) evaluation.
C            0 ... Standard 1/r operator.
C            N ... Short range 1/r operator with N/10000
C   MMMMMNNNNN ... Short range 1/r operator with NNNNN/10000 and long range
C                  1/r operator with MMMMM/10000.
C
C2IOp(127)
C     IOp(127) ... Omega for short/long range 1/r operator in E(2,AB)
C                  evaluation.
C            0 ... Standard 1/r operator.
C   MMMMMNNNNN ... Short range 1/r operator with NNNNN/10000 and long range
C                  1/r operator with MMMMM/10000.
C
C2IOp(128)
C     IOp(128) ... Coefficients of short/long range combination of 1/r
C                  operator in E(2,AA) and E(2,BB) evaluation.
C            0 ... Standard 1/r operator.
C   MMMMMNNNNN ... NNNNN/10000 short range and MMMMM/10000 long range.
C                  The signs can be changed by IOp(133) (see below).
C
C2IOp(129)
C     IOp(129) ... Coefficients of short/long range combination of 1/r
C                  operator in E(2,AB) evaluation.
C            0 ... Standard 1/r operator.
C   MMMMMNNNNN ... NNNNN/10000 short range and MMMMM/10000 long range.
C                  The signs can be changed by IOp(134) (see below).
C
C2IOp(130)
C     IOp(130) ... Coefficient of full range of HF exchange.
C            0 ... Standard full range HF exchange.
C           -1 ... 0 full range coefficient.
C        NNNNN ... NNNNN/10000 full range coefficient.
C       100000 ... Use the negative of the short range coeff as set by IOp(119).
C       200000 ... Set the short range coeff to zero.
C      1000000 ... Use the negative of the long range coeff as set by IOp(119).
C      2000000 ... Set the long range coeff to zero.
C     10000000 ... Use the negative of the mid range coeff as set by IOp(138).
C     20000000 ... Set the mid range coeff to zero.
C
C2IOp(131)
C     IOp(131) ... Coefficient of full range of DFT exchange.
C            0 ... Standard full range DFT exchange.
C           -1 ... 0 full range coefficient.
C        NNNNN ... NNNNN/10000 full range coefficient.
C       100000 ... Use the negative of the short range coeff as set by IOp(120).
C       200000 ... Set the short range coeff to zero.
C      1000000 ... Use the negative of the long range coeff as set by IOp(120).
C      2000000 ... Set the long range coeff to zero.
C
C2IOp(132)
C     IOp(132) ... Coefficient of full range of DFT correlation.
C            0 ... Standard full range DFT correlation.
C           -1 ... 0 full range coefficient.
C        NNNNN ... NNNNN/10000 full range coefficient.
C       100000 ... Use the negative of the short range coeff as set by IOp(121).
C       200000 ... Set the short range coeff to zero.
C      1000000 ... Use the negative of the long range coeff as set by IOp(121).
C      2000000 ... Set the long range coeff to zero.
C     10000000 ... Use the negative of the mid range coeff as set by IOp(138).
C     20000000 ... Set the mid range coeff to zero.
C
C2IOp(133)
C     IOp(133) ... Coefficient of full range of 1/r operator in E(2,AA)
C                  and E(2,BB) evaluation.
C            0 ... Standard full range 1/r operator.
C           -1 ... 0 full range coefficient.
C        NNNNN ... NNNNN/10000 full range coefficient.
C       100000 ... Use the negative of the short range coeff as set by IOp(128).
C      1000000 ... Use the negative of the long range coeff as set by IOp(128).
C
C2IOp(134)
C     IOp(134) ... Coefficient of full range of 1/r operator in E(2,AB)
C                  evaluation.
C            0 ... Standard full range 1/r operator.
C           -1 ... 0 full range coefficient.
C        NNNNN ... NNNNN/10000 full range coefficient.
C       100000 ... Use the negative of the short range coeff as set by IOp(129).
C      1000000 ... Use the negative of the long range coeff as set by IOp(129).
C
C2IOp(135)
C     IOp(135) ... Setup for semi-empirical.
C            0 ... Default (1 for AM1/PMn full-matrix, 2 for sparse and other
C                  methods).
C            1 ... New code.
C            2 ... Old code.
C           Nx ... Flags for AM1Par (default 2020)
C           10 ... Generate standard parameters.
C           20 ... Read parameters from rwf.
C           30 ... Read parameters from chk.
C           40 ... Read parameters from chk if present; otherwise generate.
C           50 ... Do not produce any standard parameters.
C          100 ... Read additional parameters from the input stream.
C          200 ... Read additional parameters from the input stream using
C                  MOPAC format and units.
C          300 ... Read additional parameters in both formats, Gaussian
C                  internal format first.
C         1000 ... Save parameters on rwf.
C         2000 ... Do not save parameters on rwf.
C
C2IOp(136)
C     IOp(136) ... Printing of semi-empirical parameters:
C       0 ... Default (2 unless IPrint>=2 or parameters read in).
C       1 ... Print parameters for elements used in this calculation.
C       2 ... Do not print parameters.
C       3 ... Print parameters for all elements.
C      00 ... Default (10)
C      10 ... Print parameters in human-readable form.
C      20 ... Print parameters in input format.
C      30 ... Print parameters in both formats.
C     000 ... Default (100)
C     100 ... Print only non-zero parameters
C     200 ... Print all parameters including zero parameters
C
C2IOp(137)
C     IOp(137) ... Control of FMM for nuclear repulsion.
C            0 ... Default:  Use for 5K or more atoms.
C            N ... Use for N or more atoms.
C           -1 ... Always use FMM.
C           -2 ... Never use FMM; necessary when doing external point
C                  charges if one coincides with a (ghost) nucleus.
C
C2IOp(138)
C     IOp(138) ... Mid-range coefficients for split-range functionals:
C   MMMMMNNNNN ... NNNNN/10000 HF and MMMMM/10000 XC.
C
C2IOp(139)
C     IOp(139) ... Check for XC functional derivatives required by the
C                  current calculation:
C           -1 ... Don't check.
C            0 ... Default (2).
C            N ... Check for functional and up to N-1th derivatives.
C
C2IOp(140)
C     IOp(140) ... Override PCM solution method.
C            0 ... Leave unchanged.
C            1 ... Force inversion.
C            2 ... Force iterative.
C            3 ... Force simultaneous in L502.
C
C2IOp(141)
C     IOp(141) ... Override PCM FoFCou accuracy parameter.
C            0 ... Leave unchanged.
C            N ... 10**-N.
C
C2IOp(142)
C     IOp(142) ... Convergence for iterative PCM solution.
C            0 ... Default, 10**-6
C            N ... 10**-N.
C
C2IOp(143)
C     IOp(143) ... Iteration limit for PCM solution.
C            0 ... Default (400)
C            N ... N.
C
C2IOp(144)
C     IOp(144) ... Threshold for discarding small surface elements.
C            0 ... Default (1.d-12).
C            N ... 10**-N.
C
C2IOp(145)
C     IOp(145) ... Generate and print atom/shell blocking information.
C            0 ... Default (2).
C            1 ... Yes.
C            2 ... No.
C           0x ... Default (70).
C           1x ... Block by atoms.
C           2x ... Block by shells.
C           3x ... Block by single elements.
C           4x ... One big block.
C           5x ... Block by functional group, using the connectivity.
C           6x ... Block shells until the desired minimum block size is reached,
C                  but keeping shells on the same atom in the same block.
C           7x ... Like 6, but permute shells based on connectivity.
C        MMMxx ... Minimum block size is MMM
C
C2IOp(158)
C     IOp(158) ... Over-ride defaults for PCM:
C           00 ... Normal default for model (26).
C            1 ... DPCM
C            2 ... CPCM
C            3 ... Isotropic non-symmetric IEFPCM
C            4 ... Anisotropic non-symmetric IEFPCM
C            5 ... Ionic non-symmetric IEFPCM
C            6 ... Isotropic symmetric IEFPCM.
C           10 ... Add spheres, default ofac=0.8 rmin=0.5
C           20 ... Do not add spheres
C
C2IOp(159)
C     IOp(159) ... Over-ride defaults for atomic densities:
C            0 ... Normal defaults.
C         NN>0 ... Default is NN.
C        -NN<0 ... Force type NN regardless of input.
C
C2IOp(160)
C     IOp(160) ... Operation of l316:
C            0 ... Default (1121).
C            1 ... Print out 2e integrals.
C            2 ... Do not print out 2e integrals.
C           10 ... Write fortran unformatted matrix element file,
C                  using the default name ("Gau-#####.EUF", where #####
C                  is the PID) in the scratch directory.
C           20 ... Do not write matrix element file.
C           30 ... Write the matrix element file, reading the file name
C                  from an input section (with terminating blank line).
C          100 ... Include only active nuclei in the molecule data on the
C                  file
C          200 ... Include all centers in the molecule data on the file.
C         1000 ... Use full size integers for labels of packed matrices.
C         2000 ... Use Integer*4 for labels of packed matrices; ignored on
C                  machines which do not support I*4.
C        10000 ... Use the same size integer labels for 4d matrices
C                  (2e integrals) as for other data.
C        20000 ... Use Integer*2 labels for 4d matrices; ignored on machines
C                  which do not support 16-bit integers.
C    100000000 ... Store binary data with no record marks, appropriate
C                  to reading in c/c++/perl/python.
C
C2IOp(161)
C     IOp(161) ... Saving/Restoring L302 results for SCF=Restart:
C            0 ... Default (22)
C            1 ... Save the XC dimensioning and orthonormal vectors
C                  on the chk file as well as the rwf.
C            2 ... Do not store on the chk file
C           10 ... Restore the information from the chk file if present.
C           20 ... Do not restore the information.
C
C2IOp(162)
C     IOp(162) ... Update dimensioning information in fragment rwf:
C            0 ... Default (122).
C            1 ... Yes.
C            2 ... No.
C           10 ... Generate orthonormal sets for each block in L302.
C           20 ... Do not generate orthonormal sets for each block.
C          1xx ... Set fragment type 1 (blocking ignoring overlap
C                  between fragments).
C          2xx ... Set fragment type 2 (GI model).
C          3xx ... Set fragment type 3 (ST model).
C
C2IOp(163)
C     IOp(163) ... Basis for Rung 3.5 functionals.
C            0 ... Default (-2).
C   LLLLMMMMNN ... NN functions, starting at exponent MMMM/100,
C                  scaling LLLL/100.  Negative to use old testing
C                  code.
C           -1 ... Standard tight value (0200500017)
C           -2 ... Standard medium value (0250500013)
C           -3 ... Same as -2 but forbid extra thresholds for m11plus.
C
C2IOp(164)
C     IOp(164) ... Set up BEBO parameters?
C            0 ... Default (No).
C            1 ... Yes with default parameters.
C            2 ... No.
C            3 ... Yes with parameters read from input.
C            4 ... Yes with parameters read from chk file.
C            5 ... Yes with parameters from chk file if present.
C            6 ... Yes with parameters read from rwf file.
C           00 ... Default parameter set (10).
C           10 ... -0.05H for every pair.
C
C2IOp(165)
C     IOp(165) ... Generate and test d/dx V = d/dx S^-1/2 for testing.
C            0 ... No.
C            1 ... Yes
C           Nx ... Use step-size 10^-N in numerical differentiation,
C                  default 4.
C          Mxx ... Use threashold 10^-M for linear dependence test,
C                  default 6.
C
C2IOp(166)
C     IOp(166) ... PCM point density:
C            0 ... Default (5 pts/A^2 for default quadrature).
C            N ... N pts/A^2.
C           -N ... TSAre=N, forces old quadrature.
C
C2IOp(167)
C     IOp(167) ... Size of core (for general basis input):
C            0 ... Default for internal basis sets, minimal if GBS.
C            N ... N-zeta core.
C
C2IOp(168)
C2IOp(169)
C     IOp(168),IOp(169) ... Bitmap of allowed prism paths if non-zero,
C            0-24 in word 168, 25-49 in 169, or 168=-1/-2 for OSOnly/MDOnly.
C
C2IOp(170)
C     IOp(170) ... Set-up XDM dispersion.
C            0 ... No. Default(3107110)
C            1 ... Yes, include XDM in the SCF energy self-consistently.
C            2 ... Yes, evaluate XDM using the converged SCF density.
C            3 ... No, but leave the rwf on disk if present.
C           10 ... Use XDM/PostG atomic volumes.
C           20 ... Use atomic volumes from the atomic densities.
C           30 ... Use atomic volumes from the old atomic densities.
C          100 ... Use XDM/PostG atomic polarizabilities.
C          200 ... Use atomic polarizabilities from the APFD dispersion.
C        NN000 ... Bitmap for which terms to include (0,1,2=C6,C8,C10)
C       100000 ... Use the usual denominators for the R^-N terms.
C       200000 ... Use the SAM form of the denominators as in APFD.
C      1000000 ... Kill the job when atomic parameters are unavailable.
C      2000000 ... Continue the calculation even if some of the atomic
C                  parameters are unavailable.
C     10000000 ... Evaluate the analytical gradients of the XDM atomic
C                  quantities and use straightforward code to assemble the gradient.
C     20000000 ... Same as 1 but use intermediated gradient code.
C     30000000 ... Evaluate all integrated gradient terms in CalDFT.
C
C2IOp(171)
C     IOp(171) ... Switching functions in XDM.
C            0 ... Default (05070744)
C            N ... Use switching function N for hole overstepping.
C           N0 ... Use switching function N for Hirshfeld regularization.
C         NN00 ... Range of NN/100 bohrs for hole overstepping switching.
C     MMNN0000 ... Range from 10^-MM to 10^-NN for Hirshfeld regularization.
C
C2IOp(172)
C     IOp(172) ... Damping/switching function for APF empirical dispersion.
C            0 ... Default (-5, see details in subroutines R6DAPF)
C
C2IOp(173)
C     IOp(173) ... Range for APF switching function.
C            0 ... Default (50)
C         NNNN ... A value of NNNN/1000 of the hard cutoff.
C
C2IOp(174)
C     IOp(174) ... S6 scale factor in Grimme's D2/D3/D3BJ dispersion.
C            0 ... Default (see subroutine R6DS6)
C           -1 ... Set S6 to 0.
C     NNNNNNNN ... A value of NNNNNNNN/1000000.
C
C2IOp(175)
C     IOp(175) ... S8 scale factor in Grimme's D2/D3/D3BJ dispersion.
C            0 ... Default (see subroutine R6DS8)
C           -1 ... Set S8 to 0.
C     NNNNNNNN ... A value of NNNNNNNN/1000000.
C
C2IOp(176)
C     IOp(176) ... SR6 scale factor in Grimme's D2/D3/D3BJ dispersion.
C            0 ... Default (see subroutine R6DSR6)
C           -1 ... Set SR6 to 0.
C     NNNNNNNN ... A value of NNNNNNNN/1000000.
C
C2IOp(177)
C     IOp(177) ... A1 parameter in Becke-Johnson damping for D3BJ and XDM.
C            0 ... Default (see subroutine R6DABJ/XDMABJ)
C           -1 ... Set A1 to 0.
C     NNNNNNNN ... A value of NNNNNN/1000000.
C
C2IOp(178)
C     IOp(178) ... A2 parameter in Becke-Johnson damping for D3BJ and XDM.
C            0 ... Default (see subroutine R6DABJ/XDMABJ)
C           -1 ... Set A2 to 0.
C     NNNNNNNN ... A value of NNNNNN/1000000 Ang.
C
C2IOp(179)
C     IOp(179) ... Few more control options for Rung 3.5 functionals auxiliary
C                  basis set expansion.
C            0 ... Default (161)
C            1 ... Use the Coulomb metric.
C            2 ... Use the overlap metric.
C           M0 ... Tolerance for the inversion of the auxiliary basis
C                  overlap, 10^-M.
C          100 ... Use normalized auxiliary functions.
C          200 ... Use unnormalized auxiliary functions.
C
C2IOp(180)
C     IOp(180) ... Exponent for Rung 3.5 correlation basis function.
C           0  ... Default (see subroutine R35Inp)
C      NNNNNN  ... An exponent of NNNNNN/10000.
C
C2IOp(181)
C     IOP(181) ... Choice of EM perturbation terms and gauge for external
C                  electric pulses
C           -3 ... Read from chk.
C           -2 ... Read from rwf.
C           -1 ... Same as zero.
C            0 ... Dafault, No external pulses
C            1 ... Electric-dipole term
C           1X ... Electric-quadrupole terms
C          1XX ... Magnetic-dipole term
C         0XXX ... Length gauge
C         1XXX ... Velocity gauge
C
C2IOp(182)
C     IOp(182) ... Scaling of AA/BB and AB components of D(2) in CIS(D).
C           -3 ... 0 for AB
C           -2 ... 0 for AA/BB
C           -1 ... 0 for both.
C            0 ... Default (1 for both).
C   MMMMMNNNNN ... MMMMM/10000 for AA/BB, NNNNN/10000 for AB.
C                  A value over 50000 for either translates to -(NNNNN-50000).
C
C2IOp(183)
C     IOp(183) ... Damping factor for w in CIS(D).
C           -1 ... 0.
C            0 ... Default (10000).
C        NNNNN ... NNNNN/10000
C
C2IOp(184)
C2IOp(185)
C     Control of absorbing boundary conditions.
C
C2IOp(186)
C     IOp(186) ... Scaling factor for the non-local XC dispersion term.
C           -1 ... 0.
C            0 ... Default (100000).
C        NNNNN ... NNNNN/100000
C
C2IOp(187)
C     IOp(187) ... Nuclear (NEO) basis, corresponding to IOp(5).
C           -1 ... No NEOs.
C            0 ... Default (-1).
C            7 ... Read general basis.
C           48 ... Uncontracted even tempered basis.
C      0000000 ... NEO hydrogen mass, default 8.
C      1000000 ... Use precise proton mass.
C      2000000 ... Use precise deuteron mass.
C      3000000 ... Use precise triton mass.
C      8000000 ... Use precise mass specified in the atom specification
C                  (/Mol/); all NEO atoms must have the same mass.
C      9000000 ... Use SHS group proton mass. instead of exact value.
C
C2IOp(188)
C     IOp(188) ... Nuclear (NEO) basis, corresponding to IOp(6) if
C                  general basis.
C
C2IOp(189)
C     IOp(189) ... Nuclear (NEO) basis corresponding to IOp(7).
C
C2IOp(190)
C     IOp(190) ... Pure vs. Cartesian functions in NEO basis:
C           0 ... Default (pure).
C           1 ... Pure.
C           2 ... Cartesian.
C?
#include "commonmol.inc"
#include "commonb2.inc"
#include "repall.inc"
#include "commonlab.inc"
#include "commonlp2.inc"
      Integer MaxIOp
      Parameter (MaxIOp=DEFMAXIOP)
      Common /IOp/    IOp(MaxIOp)
      Common/IO/In,IOut,IPunch
      Common/PhyCon/PhyCon(30)
      Common/MUnit/IUnit(20)
      Common/SymInf/IDum1,IDum2,JTrans(3,8),TRot(3,3),TrVec(3)
      Logical RdBas, RdChk, NDDO, PrtNDO, ChkIBF, BadBas, TstBas, AM1,
     $  RdExtB, OptGen, DelSml, DoSplt, SEHam, DFT, DFTBHm, HarTHm,
     $  RdExDB, DoMaDB, UncoAO, UncoDB, AppNDB, DoElCv, DoNElC, ZINDO,
     $  AllowP(50), FMM, ChkHyp, PCM, Strict, UseChk, NewSE, AM1Ham,
     $  UpdFrg, R35XC, RDupAO, RDupDB, LoadRN, XSC, OldXSC, HavBas,
     $  RdExNE, DoPMM
      Integer FilNum, FMFlag, FMFlg1, FMFl1I, ECPAcc, ChToWP
      Parameter (IRWGen=501,IRwBUM=506,IRwECP=512,IRwSt=520,IOFInf=521,
     $  IONDDO=566,IRwFX=584,IRwFFX=585,IRwS1=588,IRwH1=589,IRwF1=590,
     $  IODens=603,IOChrg=605,IOEDrv=619,IRwMag=695,IRwBSv=716,
     $  IRwDSv=718,IOFDPr=721,IODFDP=722,IRwDBF=724,IONuDB=732,
     $  IRwSCS=761,IRwPPT=787,IRwBEB=791,IRwNAC=810,IRwPls=819,
     $  LRwGen=1000,MaxEl=200,MaxLM=4,MOmega=DEFMOMEGA,NOmega=DEFNOMEGA,
     $  MxXCNm=DEFMAXXCNAME,MaxCrd=1024,MinPFr=1,NFld=3,
     $  MaxMM=((MaxLM+1)*(MaxLM+2)*(MaxLM+3))/6,NScaDF=DEFNSCADF)
      Character*(MxXCNm) ExName, CoName, N6Name, NameS*1, NS(1)*1
      Dimension V(*),Gen(LRwGen),CelVec(3,3),ScaDFX(NScaDF),XMask(2),
     $  EField(0:MaxMM-1,NFld),XX(1),JJ(1),Omega(MOmega,NOmega),
     $  SOScal(1),ICard(MaxCrd),LMaxEl(MaxEl),ZNucEl(MaxEl)
      Save Zero,XX,JJ,IJunk,RJunk,Pt5,One,Two,SOScl0
      Data Zero/0.0d0/,XX/0.0d0/,JJ/0/,IJunk/0/,RJunk/0.0d0/,Pt5/0.5d0/,
     $  One/1.0d0/,Two/2.0d0/,Ten/10.0d0/,SOScl0/-1.23456d10/
 1010 Format(' The combination of multiplicity',I2,
     $  ' and',I6,' electrons is impossible.')
 1180 Format(1X,I5,' basis functions, ',I5,' primitive gaussians, ',
     $I5,' cartesian basis functions',
     $/,1X,I5,' alpha electrons    ',I5,' beta electrons',
     $/,1X,5X,' nuclear repulsion energy',F21.10,' Hartrees.')
 1900 Format(' Dump of /Label/Label:')
 1903 Format(' Dump of /Label/ITtile:')
 1907 Format(' Dump of /Label/IRtCrd:')
 1910 Format(' There are',I4,' occupied orbitals but only',I4,
     $       ' basis functions!')
 1920 Format(' Nuclear repulsion after scaling=',F20.10,' Hartrees.')
 1930 Format(' Nuclear repulsion after ',A,' charges =',F20.10,
     $  ' Hartrees.')
 1980 Format(' Background charge distribution read from rwf.')
 1990 Format(' Self-consistent solvent charge distribution read from',
     $  ' rwf.')
 2030 Format(' Setting ILSW derivative flag to',I10,'.')
 2110 Format(' AtZEff=',5F15.7)
 2130 Format(' Nuclear repulsion from inactive atom pairs=',F20.10,
     $       ' Hartrees.')
 2140 Format(' PCM non-electrostatic energy = ',F20.10,' Hartrees.')
 2150 Format(' Nuclear repulsion after PCM non-electrostatic terms =',
     $  F20.10,' Hartrees.')
 2160 Format(' Unrecognized IOp(74)=',I20,'.')
 2170 Format(' Pair potential = ',F20.10,' Hartrees.')
 2180 Format(' Nuclear repulsion after pair potential =',F20.10,
     $  ' Hartrees.')
 2190 Format(' Using SOScal=',1PD16.8)
C
      Call Drum(V,MDV)
      Call GetNOp(NOp1,NOp2)
      Call InitB(IOp(34),1,MM,MMDF)
      Call FileIO(2,-IRwGen,LRwGen,Gen,0)
      IBasis = IOp(5)
      If(IBasis.eq.-1) IBasis = 0
      NGIC = IOp(6)
      NDDO = Mod(IOp(41),10).ne.0
      NDMean = (Mod(IOp(41),100)/10) - 1
      NDPar = Mod(IOp(41),1000)/100
      If(NDPar.eq.0) NDPar = 5
      ISlatr = Mod(IOp(41),10000)/1000
      IHSE = Mod(IOp(41),100000)/10000 - 1
      If(IHSE.eq.-1) IHSE = 0
      IXNDO = Mod(IOp(41),100000000)/100000
      XSC = IXNDO.ge.17.and.IXNDO.le.19
      ZINDO = IXNDO.ge.3.and.IXNDO.le.4
      IHMeth = IXNDO + 1
      IHukP = Mod(IOp(112),1000000)
      OldXSC = (Mod(IOp(112),10000000)/1000000).eq.1
      IRdPar = Mod(IHukP,100)/10
C     DFTB=Chk requires GenChk.
      If(XSC.and..not.OldXSC.and.IBasis.eq.2.and.
     $  (IRdPar.eq.3.or.IRdPar.eq.4)) then
        IBasis = 7
        If(IRdPar.eq.3) then
          NGIC = 1
        else
          NGIC = 2
          endIf
        endIf
      LocB0 = Mod(NGIC,10)
      LocB1 = Mod(NGIC,100)/10
      RdBas = IBasis.eq.7.and.LocB0.gt.0
      RdChk = IBasis.eq.7.and.LocB0.eq.2
      IPolar = IOp(7)
      If(IPolar.eq.-1) then
        IPolar = 0
      else if(IBasis.eq.14) then
        If(IOp(7).lt.0) IPolar = 10
        If(ICharg.lt.0.and.Mod(IPolar,10).eq.0) IPolar = IPolar + 1
         endIf
C
C     Determine 5D/6D and 7F/10F.  This will come from the basis if
C     it is being read from the rwf or chk file, or from IOp's if
C     a new basis is being generated.
C
      MZCore = IOp(167)
      If(RdChk) then
        Call RdWrB(3,-IUnit(9))
        IPureD = I5DB1
        IPureF = I7FB1
        If(MZCore.eq.0) then
          Call XLSW(-2,28,1,MZCore,V)
          MZCore = MZCore + 1
          endIf
      else if(RdBas) then
        If(LocB1.eq.2) then
          Call RdWrB(3,IRwBUM)
        else if(LocB1.eq.1) then
          Call RdWrB(3,IRwBSv)
        else
          Call RdWrB(3,0)
          endIf
        IPureD = I5DB1
        IPureF = I7FB1
        If(MZCore.eq.0) then
          Call ILSW(2,28,MZCore)
          MZCore = MZCore + 1
          endIf
      else
        IPureD = Mod(IOp(8),10) - 1
        IPureF = Mod(IOp(8),100)/10 - 1
        endIf
      If(IPureD.lt.-1.or.IPureD.gt.1.or.IPureF.lt.-1.or.IPureF.gt.1)
     $  Call GauErr('Illegal IOp(8) in GInput.')
      RdExtB = Mod(IOp(10),10).eq.1
      IDoMas = Mod(IOp(10),100)/10
      DoSplt = (Mod(IOp(10),1000)/100).eq.1
      IMxEq0 = Mod(IOp(10),10000)/1000
      IMxEq1 = Mod(IOp(10),100000)/10000
      If(IMxEq0.eq.1) then
        LMxEqB = 0
      else if(IMxEq1.eq.1) then
        LMxEqB = 1
      else
        LMxEqB = -1
        endIf
      IUnCon = Mod(IOp(10),1000000)/100000
      UnCoAO = IUnCon.eq.1.or.IUnCon.eq.3.or.IUnCon.eq.4.or.
     $  (IUnCon.ge.6.and.IUnCon.le.8)
      UnCoDB = IUnCon.eq.2.or.IUnCon.eq.3.or.IUnCon.eq.5.or.
     $  (IUnCon.ge.6.and.IUnCon.le.8)
      RDupAO = UnCoAO.and.IUnCon.ne.4.and.IUnCon.ne.6.and.IUnCon.ne.8
      RDupDB = UnCoDB.and.IUnCon.ne.5.and.IUnCon.ne.7.and.IUnCon.ne.8
      IFCMod = Mod(IOp(10),10000000)/1000000
      RdExDB = Mod(IOp(88),10).eq.1
      DoMaDB = (Mod(IOp(88),100)/10).eq.1
      IMxEq0 = Mod(IOp(88),10000)/1000
      IMxEq1 = Mod(IOp(88),100000)/10000
      If(IMxEq0.eq.1) then
        LMxEqD = 0
      else if(IMxEq1.eq.1) then
        LMxEqD = 1
      else
        LMxEqD = -1
        endIf
      ISETyp = IOp(12)
      IXYZ = IOp(14)
      If(IXYZ.eq.-7) IXYZ = 0
      RMag = GFloat(IOp(15)) / GFloat(10000)
      JPseud = Mod(IOp(16),10)
      JECPCh = Mod(IOp(16),100)/10
      ISECPS =  Mod(IOp(16),1000)/100
      JReadP = IOp(17)
      If(IBasis.eq.30) then
        JECPCh = 1
        JPseud = 0
        JReadP = 0
      else if(JECPCh.eq.0) then
        JECPCh = 1
        endIf
C     If(JPseud.eq.1.and.JReadP.eq.8.and..not.
C    $  (IBasis.eq.7.and.NGIC.ne.2))
C    $  Call GauErr('Pseudo=Read along with default ECPs.')
      If(JReadP.eq.-2) JReadP = 0
      If(JPseud.ne.2.and.JReadP.eq.0) then
        If(IBasis.eq.2) then
          JReadP = 0
        else if(IBasis.eq.3) then
          JReadP = 1
        else if(IBasis.eq.6.and.NGIC.le.1) then
          JReadP = 4
        else if(IBasis.eq.6) then
          JReadP = 5
        else if(IBasis.eq.7) then
          If(NGIC.gt.0) then
            JReadP = 8
          else
            JPseud = 2
            endIf
        else if(IBasis.eq.8) then
          If(NGIC.eq.2) then
            JReadP = 2
          else
            JPseud = 2
            endIf
        else if(IBasis.eq.9) then
          JReadP = 3
        else if(IBasis.eq.17) then
          If(NGIC.eq.8) then
            JReadP = 15
          else
            JReadP = 8 + NGIC
            endIf
        else if(IBasis.eq.26) then
          If(NGIC.eq.1) then
            JReadP = 19
          else if(NGIC.eq.2.and.JReadP.eq.0) then
            JPseud = 2
          else
            JReadP = 20
            endIf
        else if(IBasis.ge.43.and.IBasis.le.45) then
          JReadP = 21
        else
          JPseud = 2
          endIf
        endIf
      JPrinP = IOp(18)
      IBOut = IOp(24)
      IBPr = Mod(IBOut,10)
      IPtGIn = Mod(IBOut,100)/10
      IPtGTb = Mod(IBOut,1000)/100
      IPtShC = Mod(IBOut,10000)/1000
      IPtNAO = Mod(IBOut,100000)/10000
      IPtNDF = Mod(IBOut,1000000)/100000
      IPtAll = Mod(IBOut,10000000)/1000000
      If(IPtNAO.eq.0) IPtNAO = 2
      If(IBPr.eq.1.or.IPtGTb.eq.1) JPrinP = 1
      LLink = IOp(25)
      IAccFD = IOp(27)
      IPrint = IOp(33)
      IDump = IOp(34)
      IniFX = IOp(39)
      I74 = IOp(74)
      I74H = Mod(I74,10000000)/1000000
      If(I74H.ge.2) then
        If(IHMeth.eq.1) then
          If(I74H.eq.2) then
            IHMeth = 13
          else if(I74H.eq.3) then
            IHMeth = 14
          else if(I74H.eq.4) then
            IHMeth = 22
          else if(I74H.eq.5) then
            IHMeth = 23
          else
            Write(IOut,2160) I74
            Call Lnk1E(0)
            endIf
          endIf
        I74 = I74 - I74H*1000000
        endIf
      AM1 = AM1Ham(IHMeth)
      INewSE = Mod(IOp(135),10)
      ICntSE = IOp(135)/10
      I2 = Mod(ICntSE,1000)/100
      If(Mod(ICntSE,10).eq.0) ICntSE = ICntSE + 1
      If(I2.eq.0) ICntSE = ICntSE + 100
      NewSE = INewSE.le.1.and.AM1
      IPrPar = IOp(136)
      NucFMM = IOp(137)
      If(IBasis.eq.30) then
        If(ISETyp.eq.0) ISETyp = -2
      else if(IXNDO.ne.0.and.SEHam(IHMeth)) then
        ISlatr = 1
        If(ISETyp.eq.0) then
          If(IXNDO.le.9) then
            ISETyp = Min(IXNDO,8)
          else if(IXNDO.eq.17.or.IXNDO.eq.18) then
            ISETyp = 9
          else if(IXNDO.eq.40) then
            ISETyp = 11
          else if(AM1Ham(IHMeth)) then
            ISETyp = 10
            endIf
          endIf
        If(IBasis.eq.0) IBasis = 2
        If(NGIC.eq.0) then
          If(ZINDO.or.NewSE) then
            NGIC = 6
          else
            NGIC = 3
            endIf
          endIf
        endIf
      If(NGIC.eq.-1) NGIC = 0
      If(ISETyp.eq.-1) ISETyp = 0
      NGIC0 = Mod(NGIC,10)
      Call ILSW(2,62,IBkgCh)
      If(Mod(IBkgCh,10).eq.1.and.Mod(IOp(43),10).ne.2) then
        IRdChg = 31
      else
        IRdChg = 3
        endIf
      If((Mod(IBkgCh,100)/10).eq.1.and.(Mod(IOp(43),100)/10).ne.2) then
        IRdSCS = 31
      else
        IRdSCS = 3
        endIf
      ChkHyp = IBasis.ne.30
      TstBas = IOp(46).eq.1.or.(IOp(46).eq.0.and.ChkHyp)
      IPFlag = IOp(47)
      FMFlag = IOp(48)
      FMFl1I = IOp(49)
      NFxFlg = IOp(51)
      ECPAcc = IOp(53)
      ICorTp = IOp(54)
      ISparse = IOp(55)
      ISpInt = IOp(56)
      ICore = IOp(57)
      ICCP = IOp(60)
      If(ICCP.eq.0) then
        If(IBasis.eq.30.or.(IBasis.eq.7.and.NGIC0.ne.0.and..not.RdExtB))
     $  then
          ICCP = -1
        else
          ICCP = 2
          endIf
        endIf
      RDupAO = RDupAO.and.ICCP.ne.-1
      RDupDB = RDupDB.and.ICCP.ne.-1
      OptGen = ICCP.ge.1
      DelSml = ICCP.gt.2
      Strict = ICCP.eq.1
      If(ICCP.gt.2) then
        Cut = GFloat(10)**(-ICCP)
      else
        Cut = Zero
        endIf
      If(ICCP.lt.0) ICCP = 0
      ICkFMM = IOp(63)
      IDerIL = IOp(64)
      IEFBas = Mod(IOp(67),10)
      ISCRF  = IOp(70)
      IDeriv = IAbs(IOp(71))
      IDenBs = IOp(82)
      NGICDB = IOp(83)
      If(IDenBs.eq.0) then
        If(IBasis.eq.7.and.NGIC.ne.0) then
          IDenBs = 7
          NGICDB = NGIC + 2
        else
          IDenBs = -1
          endIf
      else if(IDenBs.eq.46.or.IDenBs.eq.47) then
        Call MapDBS(IBasis,NGIC,IPolar,IDenBs,NGICDB)
        If(IDenBs.eq.-1)
     $    Call GauErr('No stored fitting set matched to this AO basis.')
        endIf
      If(IPtNDF.eq.0) then
        If(IDenBs.eq.29) then
          IPtNDF = 1
        else
          IPtNDF = 3
          endIf
        endIf
      IPolDB = IOp(84)
      IPurDB = Max(IOp(85)-1,0)
      LMaxBF = IOp(86) - 1
      LMaxDF = IOp(87) - 1
      IDKHam = Mod(IOp(91),100)
      IDKHam = Max(IDKHam-1,0)
      ISclSO = Mod(IOp(91),1000)/100
      ISclS1 = Mod(IOp(91),10000)/1000
      IAONrm = Mod(IOp(92),10)
      If(IAONrm.eq.0) IAONrm = 2
      IDBNrm = Mod(IOp(92),100)/10
      If(IDBNrm.eq.0) IDBNrm = 3
      IAONCC = Mod(IOp(92),1000)/100
      If(IAONCC.eq.0) IAONCC = 2
      IDBNCC = Mod(IOp(92),10000)/1000
      If(IDBNCC.eq.0) IDBNCC = 3
      INucDB = Mod(IOp(93),100)
      If(INucDB.eq.0) then
        If(IDKHam.gt.0) then
          INucDB = 2
        else
          INucDB = 1
          endIf
        endIf
      IApNDB = Mod(IOp(93),1000) / 100
      AppNDB = IApNDB.eq.1
      LoadRN = (Mod(IOp(93),1000000)/100000).eq.1
      ISAO = Mod(IOp(113),10)
      If(ISAO.eq.0) ISAO = 1
      ISDBF = Mod(IOp(113),100)/10
      If(ISDBF.eq.0) ISDBF = 2
      ISNEO = Mod(IOp(113),1000)/100
      If(ISNEO.eq.0) ISDBF = 1
      IOpClI = IOp(116)
      IPhase = IOp(123)
      Call ILSW(1,61,IPhase)
      IPrtBl = Mod(IOp(145),10)
      IBkTyp = Mod(IOp(145),100)/10
      If(IPrtBl.eq.0) IPrtBl = 2
      MinBkS = IOp(145)/100
      IDfPCM = IOp(158)
      ITyADD = IOp(159)
      IUpdFr = Mod(IOp(162),10)
      UpdFrg = IUpdFr.eq.1
      IFrgT1 = Mod(IOp(162),100)/10
      IFrgTp = Mod(IOp(162),1000)/100
      IBEBO = IOp(164)
      IEMPul = IOp(181)
      NEOBas = Mod(IOp(187),1000000)
      If(NEOBas.eq.0.or.XSC.or.IBasis.eq.30) NEOBas = -1
      INEOM = Mod(IOp(187),10000000)/1000000
      If(INEOM.eq.0) INEOM = 8
      NGICNE = IOp(188)
      If(NEOBas.eq.7.and.NGICNE.eq.0.and.IBasis.eq.7) NGICNE = NGIC
      IPolNE = IOp(189)
      IPurNE = Max(IOp(190)-1,0)
      RdExNE = (Mod(IOp(88),1000000)/100000).eq.1
C
C     Conver NDPar to internal form (see SetNDO).
      If(NDPar.eq.1.and.IBasis.eq.7.and.NGIC.eq.1) NDPar = -2
      If(NDPar.eq.1.and.IBasis.eq.7.and.NGIC.eq.2) NDPar = -3
      If(NDPar.eq.2) NDPar = -2
      If(NDPar.eq.3) NDPar = -3
      If(NDPar.eq.1) NDPar = -1
      If(NDPar.ge.4) NDPar = NDPar - 4
      If(NewSE) Call AM1Par(In,IOut,IPrint,IHMeth,ICntSE,IPrPar,NAtoms,
     $  IAn,IAtTyp,V,MDV)
      ISEAll = 1
      Call GetSEP(IOut,IPrint,IHMeth,0,0,.True.,ISEAll,MaxASE,NSEAtI,
     $  NSEAtR,NCorPA,MaxCPA,NSEElI,NSEElR,MxElBE,jRSESc,jISEP,jRSEP,
     $  jSECP,jIElP,jRElP,jBELP,jBEBO,jNext,V,MDV)
      Call RdChg(IRdChg,0,NumChg,NumDBF,LChg,LDBF,jNext-1,jChg,jDBF,IV,
     $  V,MDV)
      Call RdChg(IRdSCS,IRwSCS,NSCSCh,NSCSDB,LSCSCh,LSCSDB,IV-1,jSCSCh,
     $  jSCSDB,jV,V,MDV)
      jNext = jV
      Call GetCon(0,NAtoms,0,JJ,JJ,jV,MxBond,IConct,INBond,IIBond,
     $  IRBTyp,IV,V,MDV)
      Call ISet(-1,MaxEl,LMaxEl)
      Call ASet(MaxEl,-One,ZNucEl)
      MDV1 = MDV - IV + 1
      If(XSC.and..not.OldXSC) then
        IBB = IV
        NBB = 0
        MaxBB = MDV1*IntPWP(0)
        Call RdRout(In,IOut,220,ICard,MaxCrd,V(IBB),MaxBB,NBB,Junk)
        IV = IV + ChToWP(NBB)
        MDV1 = MDV - IV + 1
        Call GetHP(IOut,IPrint,.False.,IXNDO,IHukP,NAtoms,NBasis,IAn,
     $    IAtTyp,AtmChg,C,V(IBB),NBB,MaxEl,LMaxEl,ZNucEl,V(IV),MDV1)
        endIf
C
C     Defaults for pure vs. Cartesian functions.  All basis sets use pure
C     f and higher by default, except the SBK basis sets used with the
C     Stuttgart/Dresden ECPs, which Dolg advises should be used with 10f,
C     since that is how they were developed.
C
      If(IPureD.eq.-1) then
        If(IBasis.eq.1.or.IBasis.eq.5.or.IBasis.eq.8.or.IBasis.eq.11.or.
     $    IBasis.eq.9.or.IBasis.eq.17.or.
     $    IBasis.eq.13.or.IBasis.eq.15.or.
     $    IBasis.eq.33.or.IBasis.eq.34.or.IBasis.eq.35) then
          IPureD = 1
        else
          IPureD = 0
          endIf
        endIf
      If(IPureF.eq.-1) then
        If(IBasis.eq.9.or.IBasis.eq.17) then
          IPureF = 1
        else
          IPureF = 0
          endIf
        endIf
      NumD = 5 + IPureD
      NumF = 7 + 3*IPureF
      Call ILSW(1,2,IPureD)
      Call ILSW(1,16,IPureF)
      Call ILSW(1,28,0)
      If(MZCore.eq.0) then
        If(IBasis.eq.8.and.NGIC.eq.0) then
          MZCore = 2
        else
          MZCore = 1
          endIf
        endIf
      Call ILSW(1,28,MZCore-1)
      If(IPrint.gt.0) Call CorPrS(IOut,'Coordinates in L301:',0,NAtoms,
     $  IAn,.True.,IAtTyp,C,PhyCon(1))
      Call MolBas(In,IOut,IPrint,IDump,1,IPtAll,NewSE,IBasis,NGIC,
     $  IPolar,IPureD,IPureF,RdExtB,.False.,.False.,IAONrm,INucDB,
     $  LMxEqB,ICore,ISlatr,ISETyp,ICharg,Multip,NE,NAtoms,IAn,IAtTyp,
     $  AtmChg,C,RNFroz,AtParM(1,MAtNRd),IHMeth,IPFlag,AllowP,MaxASE,
     $  NSEAtI,V(jISEP),NSEAtR,V(jRSEP),MaxEl,LMaxEl,MM,MMDV,V(IV),MDV1)
      If(IDoMas.ne.0.or.UnCoAO.or.DoSplt) then
        Call RdWrB(5,0)
        I5DB2 = IPureD
        I7FB2 = IPureF
        Call Renor1(IOut,(IAONrm-1),0,1,1,NShelB,ShelNB,ShelTB,ShelCB,
     $    ShelAB,ShlADB,JAnB,EXXB,C1B,C2B,C3B,C4B)
        Call RdWrB(-1,IRwBUM)
        endIf
      Call GetPBC(NPDir,CelVec)
C
C     Load pseudo-potentials (ECPs) if requested.  This must be done
C     before checking the charge and multiplicity, since the number of
C     electrons replaced by ECPs must be known.
C
      Call RdWrECP(8,0,NAtoms)
      If(JPseud.ne.2.and.JReadP.ne.0) then
        If(IBasis.eq.7) then
          NGICP = NGIC
        else
          NGICP = 0
          endIf
        Call PInput(IOut,IPrint,JReadP,JPrinP,NGICP,JECPCh,ISECPS,
     $    NAtoms,IAn,IAtTyp)
      else if(ITqry(IRwECP).gt.0) then
        Call FileIO(5,IRwECP,0,V,0)
        endIf
C
C     Massage after possibly loading RNFroz with ECPs.
C
      If(IDoMas.ne.0) then
        NAtOld = NAtoms
        Call RWCon(3,0,NAtoms,0,MxBond,JJ,JJ,XX)
        Call Massge(IOut,IDump,IDoMas,.True.,MM,MMDF,NPDir)
        If(NAtoms.ne.NAtOld) then
          If(ITqry(IRwECP).gt.0) Call RdWrECP(1,IRwECP,NAtoms)
          If(MxBond.gt.0) then
            Call AllCon('GInput',0,MxBond,NAtoms,IV,ICon,INBond,IIBond,
     $        IRBTyp,jNext,V,MDV)
            Call DumCon(NAtOld,NAtoms,MxBond,V(INBond),V(IIBond),
     $        V(IRBTyp))
            endIf
          endIf
        endIf
      If(IFCMod.gt.0) then
        LMaxFC = 2
        IAtMx = IV
        IEnd = IAtMx + (LMaxFC+1)*NAtoms - 1
        Call TstCor(IEnd,MDV,'GInput-FCModB')
        Call FCModB(IOut,IPrint,IFCMod,NAtoms,IAn,IAtTyp,C,LPSkip,MM,
     $    MMDF,LMaxFC,V(IAtMx))
        endIf
      I5DB1 = IPureD
      I7FB1 = IPureF
      IPtFlg = 100*IPtNAO + 10
      Call ILSW(2,26,NoSym)
      Call PostPB(In,IOut,IPrint,IDump,1,RDupAO,LMxEqB,LMaxBF,OptGen,
     $  Strict,DelSml,Cut,.True.,.True.,.False.,UnCoAO,DoSplt,IBPr,
     $  IEFBas,IOFInf,IRwGen,NAtoms,IAn,IAtTyp,AtmChg,C,NBasis,NPrims,
     $  IXYZ,MaxLM,NFld,RMag,EField,IAONrm,1,IPtGIn,IPtGTb,IPtFlg,
     $  IPtShC,NoSym,NOpAll,NOp2,ISAO,Gen,V(IV),MDV1)
      Call PutLbl(IPureD,IPureF,1,0,-2)
      If(IDump.ne.0) then
        Write(IOut,1900)
        LenP = LenStr(5,Label,0)
        Call StrOut(IOut,Label,LenP,1)
        Write(IOut,1903)
        LenP = LenStr(5,ITitle,0)
        Call StrOut(IOut,ITitle,LenP,1)
        Write(IOut,1907)
        LenP = LenStr(5,IRtCrd,0)
        Call StrOut(IOut,IRtCrd,LenP,1)
        If(IDump.ge.2) Call BDump(2)
        endIf
      If(JPseud.ne.2.and.IPtGIn.eq.1) Call PrnPot(IOut,IPrint,2,MxAtEC,
     $  MxTECP,1,NAtoms,IAn,NLP,CLP,ZLP,KFirst,KLast,LMax,LPSkip,RNFroz,
     $  0)
C
C     Compute NAE, NBE, NE, and the nuclear repulsion energy.
C     Check the charge and multiplicity.
C
      If(IBasis.eq.2.or.IBasis.eq.30.or.ISETyp.eq.-2.or.ISETyp.gt.0)
     $  Call VSTONF(IOut,NewSE,ISETyp,NAtoms-NPtMol,IAn,IAtTyp,AtZNuc,
     $  MaxASE,NSEAtI,V(jISEP),MaxEl,ZNucEl,RNFroz)
      Call LDAtmC(IOut,IPrint,NAtoms,NPtMol,ICharg,IAn,IAtTyp,RNFroz,
     $  AtChMM,AtZNuc,NE,AtmChg)
      If(Mod(NE,2).eq.Mod(Multip,2)) then
        Write(IOut,1010) Multip, NE
        Call Lnk1E(0)
        endIf
      If(Multip.ge.0) then
        NAE = (NE+Multip-1)/2
        NBE = (NE-Multip+1)/2
      else
        NAE = (NE+Multip+1)/2
        NBE = (NE-Multip-1)/2
        endIf
C
C     Initialize the PE environment, if present
C
      Call GtCPar(IOut,IPrint,0,NAtoms,IAn,IAtTyp,MaxAtm,IAtPrM,AtParM,
     $  IV,NIAtP,NRAtP,NDAtP,jIAtP,jRAtP,jV,V,MDV)
      Call MMPIni(IOut,IPrint,NAtoms,NDAtP,IAtTyp,C,V(jIAtP),V(jRAtP),
     $  V(IConct),jV,DoPMM,LPMAll,IPMAll,IV,V,MDV)
      MDV1 = MDV - IV + 1
C
C     Set up an NEO basis.
C
      If(HavBas(IOut,IPrint,3,IBasis,NGIC,NEOBas,NGICNE)) then
        Call MolBas(In,IOut,IPrint,IDump,3,IPtAll,NewSE,NEOBas,NGICNE,
     $    IPolNE,IPurNE,IPurNE,RdExNE,.False.,.False.,IAONrm,INucDB,
     $    LMxEqB,ICore,ISlatr,ISETyp,ICharg,Multip,NE,NAtoms,IAn,IAtTyp,
     $    AtmChg,C,RNFroz,AtParM(1,MAtNRd),IHMeth,IPFlag,AllowP,MaxASE,
     $    NSEAtI,V(jISEP),NSEAtR,V(jRSEP),MaxEl,LMaxEl,MM,MMDV,V(IV),
     $    MDV1)
C       If(DoMaDB) Call Massge(IOut,IDump,1,.False.,MM,MMDF,NPDir)
        IPtFlg = 100*IPtNAO + 30
        Call PostPB(In,IOut,IPrint,IDump,3,RDupDB,LMxEqB,LMaxBF,.False.,
     $    Strict,DelSml,Cut,.False.,.False.,.False.,.False.,DoSplt,0,
     $    IEFBas,IOFInf,IRwGen,NAtoms,IAn,IAtTyp,AtmChg,C,NNEOBs,NNEOPr,
     $    IXYZ,MaxLM,NFld,RMag,EField,IAONrm,1,IPtGIn,IPtGTb,IPtFlg,0,
     $    NoSym,NOpAll,NOp2,ISDBF,Gen,V(IV),MDV1)
      else
        NEOBas = -1
        endIf
      Call SetOpC(IOpClI,Multip,IDeriv,NEOBas)
      Call NEOMol(IOut,IPrint,NEOBas,.True.,NAtoms,IAtTyp)
      Call NEOMas(IOut,IPrint,NEOBas,INEOM,NAtoms,IAn,IAtTyp,AtmWgt)
      If(NEOBas.gt.0.and.ITqry(IRwECP).gt.0) then
        Call FixECP(JECPCh,NAtoms,IAtTyp)
        Call RdWrECP(1,IRwECP,NAtoms)
        endIf
C
C     Handle semi-empirical parameters and then compute classical
C     (nuclear repulsion and poossibly SE corrections).
C
      If(XSC.and..not.OldXSC) Call GetHP(IOut,IPrint,.True.,IXNDO,IHukP,
     $  NAtoms,NBasis,IAn,IAtTyp,AtmChg,C,V(IBB),NBB,MaxEl,LMaxEl,
     $  ZNucEl,V(IV),MDV1)
      IUseCN = ICNUse(0,NAtoms,IAtTyp,MicOpt)
      If(NewSE) then
        Call AM1Nuc(IOut,IPrint,IUseCN,IHMeth,0,1,NAtoms,IAn,IAtTyp,
     $    AtmChg,C,MaxASE,V(ISEAll),V(jRSESc),NSEAtI,V(jISEP),NSEAtR,
     $    V(jRSEP),NCorPA,MaxCPA,V(jSECP),NSEElI,V(jIElP),NSEElR,
     $    V(jRElP),MxBond,V(INBond),V(IIBond),V(IRBTyp),ENR,XX,XX,V(IV),
     $    MDV1)
      else if(DoPMM) then
C       For PMM, do the nuclear repulsion with the static multipole
C       contributions.
        Call MMPNuc(IOut,IPrint,.False.,NAtoms,AtmChg,C,NDAtP,V(jRAtP),
     $    LPMAll,V(IPMAll),ENR,EField(1,1),XX,V(IV),MDV1)
        Call SavPMM(IOut,IPrint,0,LPMAll,V(IPMAll))
      else if(.not.XSC) then
        Call RepNN(IOut,IPrint,.True.,0,NucFMM,.False.,NAtoms,IAtTyp,
     $    AtmChg,C,ENR,0,XX,V(IV),MDV1)
        endIf
      Gen(41) = ENR
      If(IUseCN.eq.1.and..not.NewSE.and..not.DFTBHm(IHMeth)) then
        Call RepNN(IOut,IPrint,.True.,-1,NucFMM,.False.,NAtoms,IAtTyp,
     $    AtmChg,C,ENRI,0,XX,V(IV),MDV1)
      else
        ENRI = Zero
        endIf
      Gen(52) = ENRI
      NAt3 = 3 * NAtoms
      If(ICkFMM.ne.0) then
        IFRef = IV
        IVF = IFRef + NAt3
        IdVF = IVF + NAtoms
        IV2 = IdVF + NAt3
        MDV2 = MDV - IV + 1
        Call TstCor(IV2,MDV,'GInput-CkPFMM')
        Call CkPFMM(IOut,IPrint,IDump,ICkFMM,1,-1,Zero,NAtoms,AtmChg,C,
     $    V(IFRef),V(IVF),V(IdVF),V(IV2),MDV2)
        endIf
C
C     Set default ZEffective values for spin-orbit, now that we know which
C     atoms have ECPs.
C
      Call LodZEf(NAtoms,IAn,AtmChg,AtZEff)
      If(IPrint.gt.0) Write(IOut,2110) (AtZEff(I),I=1,NAtoms)
C
C     Compute number of cartesian functions.
C
      Call GetMBN(1,1,NShell,ShellT,ShellC,NBas6D)
C
C     Check that there are enough basis functions.
      If(NAE.gt.NBasis) then
        Write(IOut,1910) NAE, NBasis
        If(TstBas) Call Lnk1E(0)
        endIf
C
C     Initialize nuclear and field derivative rwfs, generalized density
C     rwf and stability test flag.
C
      NAt3P = NAt3 + 3*NPDir
      If(IniFX.eq.0) then
        Call TossMO(0,.True.)
        Call ILSW(1,5,0)
        Call ILSW(1,20,0)
        Call ILSW(1,56,0)
        Call ILSW(1,15,0)
        Call ILSW(2,51,NMtPBC)
        NMtPBC = Max(NMtPBC,1)
        Call ILSW(1,51,NMtPBC)
        Call ILSW(1,71,0)
        IDSet = IDerIL
        If(IDSet.eq.-2) IDSet = 0
        If(IDerIL.ne.0) then
          Write(IOut,2030) IDSet
          Call ILSW(1,24,IDSet)
          endIf
        Gen(62) = 0
        LenEDv = 1 + 3 + 3*NAt3P + 6 + 6*NAt3P + 10 + 3*NAt3P + 3*NAt3P
     $    + 9 + 9*NAt3P + 18 + 18*NAt3P + 6 + 6*NAt3P
        Call TstCor(IV+LenEDv,MDV,'GInput-InitEDv')
        Call IniBc2(IRwFx,2*NAt3P,V(IV),MDV1)
        LFFx = ITqry(IRwFFX)
        If(LFFX.gt.0) Call FileIO(5,IRwFFX,0,V(IV),0)
        Call IniBc2(IOEDrv,LenEDv,V(IV),MDV1)
        LenMag = 28 + 31*(NAt3P/3)
        Call IniBc2(IRwMag,LenMag,V(IV),MDV1)
        Call IniBc2(IRwNAC,NAt3P,V(IV),MDV1)
        Call AClear(10,V(IV))
        Call FileIO(1,-IRwSt,InToWP(10),V(IV),0)
        Call ConDDF(IODens,1)
        Call BNPack(0,-1,V(IV))
        Call FileIO(1,-IODens,1,V(IV),0)
        Call BNPack(0,0,V(IV))
        Call ConDDF(IOChrg,1)
        Call FileIO(1,-IOChrg,1,V(IV),0)
        Call ILSW(1,6,0)
        If(ITqry(IOFDPr).gt.0) then
          Call FileIO(2,-IOFDPr,1,XMask,0)
          Call BNUPak(XMask,MskFDP,NFrqRd)
          Call BNPack(0,NFrqRd,XMask)
          Call FileIO(1,-IOFDPr,1,XMask,0)
          endIf
        If(ITqry(IODFDP).gt.0) then
          Call FileIO(2,-IODFDP,2,XMask,0)
          Call BNUPak(XMask,MskFDP,NFrqRd)
          Call BNUPak(XMask(2),NDeriv,ND12)
          Call BNPack(0,NFrqRd,XMask)
          Call BNPack(0,0,XMask(2))
          Call FileIO(1,-IODFDP,2,XMask,0)
          endIf
C       Clear these files because some links decide what to do based on
C       their existence or length.
        If(ITqry(IRwS1).gt.0) Call FileIO(5,IRwS1,0,V,0)
        If(ITqry(IRwH1).gt.0) Call FileIO(5,IRwH1,0,V,0)
        If(ITqry(IRwF1).gt.0) Call FileIO(5,IRwF1,0,V,0)
        endIf
C
C     Set up NDDO parameters.
C
      If(NDDO) then
        IFlags = IV
        IScVEN = 2
        IScVEE = 3
        IScVNN = 4
        IHDiag = 5
        IBeta = IHDiag + NBasis
        IFirst = IBeta + NBasis
        IShlBf = IFirst + InToWP(NBasis)
        LenNDO = IShlBf + InToWP(NBasis) - IV
        Call TstCor(IV+LenNDO,MDV,'GInput-NDDO')
        If(NDPar.eq.-2) Call FileIO(2,-IONDDO,LenNDO,V(IFlags),0)
        If(NDPar.eq.-3)
     $    Call FileIO(2,-FilNum(IONDDO,IUnit(9)),LenNDO,V(IFlags),0)
        If(NDPar.le.-2) Call BNUPak(V(IFlags),INDOF,Junk)
        If(NDPar.le.-2) IHSE = Mod(INDOF,4) / 2
        If(NDPar.le.-2) NDMean = Mod(INDOF,2)
        PrtNDO = IBPr.ne.0.or.NDPar.eq.-1.or.IPrint.gt.0
        Call SetNDO(In,IOut,PrtNDO,NDMean,IHSE,NDPar,NShell,NBasis,IAn,
     $    JAn,ShellT,ShellC,V(IScVEN),V(IScVEE),V(IScVNN),
     $    V(IHDiag),V(IBeta),V(IFirst),V(IShlBf))
        INDOF = 2*IHSE + NDMean
        If(NDPar.ge.-1) Call BNPack(INDOF,0,V(IFlags))
        Call FileIO(1,-IONDDO,LenNDO,V(IFlags),0)
        ENRSca = Gen(41) * V(IScVNN)
        Gen(41) = ENRSca
        endIf
C
      Write(IOut,1180) NBasis, NPrims, NBas6D, NAE, NBE, ENR
      If(NAE.lt.0.or.NBE.lt.0) Call GauErr('Bad numbers of electrons.')
      If(IPSave(0).eq.0.and.ENRI.ne.Zero) Write(IOut,2130) ENRI
      If(NDDO) Write(IOut,1920) ENRSca
C
C     Check dimensioning limits.
C
      ChkIBF = LLink.gt.1
      I1 = IV
      I2 = I1 + InToWP(NAtoms)
      Call TstCor(I2+InToWP(NAtoms)-1,MDV,'GInput-BadBas')
      Call RdWrB(3,0)
      If(BadBas(IOut,IPrint,TstBas,ChkIBF,ChkHyp,DoSplt,V(I1),V(I2)))
     $  Call GauErr('Error in basis set or dimensions detected.')
C
C     Generate/Read control information for PRISM, CalDFT, FMM, and ECPs.
C
      NBsUse = NBasis
      I75 = IOp(75)
      If(I75.eq.0) I75 = -5
      If(IFCMod.ne.0.and.I75.eq.-4) then
        I75 = 199302
      else if(IFCMod.ne.0.and.I75.eq.-5) then
        I75 = 199590
      else if(IFCMod.ne.0.and.I75.eq.-7) then
        I75 = 299974
      else if(I75.le.-4.and.I75.ge.-999) then
        I75 = IAbs(I75)
        endIf
      Call GetDFT(IOut,IPrint,.True.,IHMeth,IOp(79),IOp(74),I75,IOp(76),
     $  IOp(77),IOp(78),IOp(80),IOp(81),IOp(107),IOp(108),IOp(109),
     $  IOp(119),IOp(120),IOp(121),IOp(125),IOp(126),IOp(127),IOp(128),
     $  IOp(129),IOp(130),IOp(131),IOp(132),IOp(133),IOp(134),IOp(124),
     $  IOp(138),IOp(139),IOp(182),IOp(183),IOp(186),IExCor,IXCFnc,
     $  N6Func,IRadAn,DFT,ExName,CoName,N6Name,ScaHFX,ScaDFX,IRanWt,
     $  IRanGd,ICorTp,IEmpDi,Omega)
      Call PsmSet(IOut,IDump,NPDir,IPFlag,FMFlag,FMFl1I,NFxFlg,
     $  ECPAcc,ICorTp,ISparse,ISpInt,IOp(106),IOp(118),IOp(168),
     $  IOp(169),IHMeth,NAtoms,NOpAll,MxAtSO,NEqAll,ScaHFX,IAtTyp)
      Call SetPFl(IOut,IPrint,0,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,
     $  IHMeth,Omega,0,Junk1,Junk2,Junk3,Junk4,V(IV),MDV1)
C     Generate Rung 3.5 aux basis any time IOp(163) is set, e.g. for L608.
      IXCR35 = IExCor
      If(IOp(163).ne.0.and..not.R35XC(IExCor,-1,-1)) IXCR35 = 64000
      If(R35XC(IXCR35,-1,-1)) Call R35Inp(IOut,IPrint,IXCR35,IOp(163),
     $  IOp(179),IOp(180),V(IV),MDV1)
      SOScal(1) = One
      If(ISclSO.gt.0) SOScal(1) = Ten**ISclSO
      If(ISclS1.eq.1) SOScal(1) = SOScal(1)*Pt5
      If(ISclS1.eq.2) SOScal(1) = SOScal(1)*Two
      If(ISclS1.eq.3) SOScal(1) = -SOScal(1)*Two
      If(ISclS1.eq.4) SOScal(1) = Zero
      If(SOScal(1).ne.One.and.(IPSave(0).eq.0.or.IPrint.gt.0))
     $  Write(IOut,2190) SOScal(1)
      If(SOScal(1).eq.Zero) SOScal(1) = SOScl0
      Call XLSW(1,7,2,IJunk,SOScal)
C
C     Initialize integral file data (/IBF/).
      IDer2E = Max((Mod(IOp(37),100)/10)-1,0)
      Call Set2E(IOut,IDump,IOp(25),IOp(11),IOp(30),IDer2E,ScaHFX)
C
C     Set up a density basis for finite-sized nuclei.
C
      If(INucDB.gt.1) then
        Call InitB(0,1,MM,MMDF)
        Call NucDBF(IOut,IPrint,INucDB,LoadRN,NAtoms,IAn,IAtTyp,AtmChg,
     $    AtParM(1,MAtNRd),C,IPFlag,AllowP,IHMeth,PhyCon(1),MM,MMDF,
     $    V(IV),MDV1)
        Call RdWrB(-12,0)
      else if(ITqry(IONuDB).gt.0) then
        Call FileIO(5,IONuDB,0,V,0)
        endIf
C
C     Set up a density basis for spherical atomic charge distributions.
C
      If(HarTHm(IHMeth)) then
        Call InitB(0,1,MM,MMDF)
        ITypAD = 23
        Call ATSDBs(IOut,IPrint,1,ITypAD,0,ICharg,NE,NAtoms,IAn,IAtTyp,
     $    AtmChg,C,.True.,MM,MMDF,V(IV),MDV1)
        Call Renor1(IOut,2,0,0,2,NShell,ShellN,ShellT,ShellC,ShellA,
     $    ShlADF,JAn,EXX,C1,C2,C3,C4)
        Call RdWrB(17,0)
        endIf
C
C     Handle external charges and self-consistent solvent charge distibution.
C
      If(NumChg.gt.0) then
        If(IPSave(0).eq.0) Write(IOut,1980)
        ICtENR = 1
        If(NucFMM.eq.-2) ICtENR = 4
        Call ExtCNR(IOut,IPrint,ICtENR,NAtoms,IAn,AtmChg,C,LChg,NumChg,
     $    V(jChg),IHMeth,FMM,FMFlag,FMFlg1,NFxFlg,IPFlag,AllowP,1,0,
     $    ECNR,XX,XX,V(IV),MDV1)
        Gen(41) = Gen(41) + ECNR
        Write(IOut,1930) 'external point',Gen(41)
        endIf
      Gen(55) = Zero
      If(NSCSCh.gt.0.or.NSCSDB.gt.0) then
        If(IPSave(0).eq.0) Write(IOut,1990)
        Call SCSCNR(IOut,IPrint,1,NumChg,LChg,V(jChg),NAtoms,IAn,IAtTyp,
     $    AtmChg,C,LSCSCh,LSCSDB,NSCSCh,NSCSDB,V(jSCSCh),V(jSCSDB),
     $    IHMeth,FMM,FMFlag,FMFlg1,NFxFlg,IPFlag,AllowP,NBasis,Zero,
     $    SCSNR,XX,V(IV),MDV1)
        Gen(41) = Gen(41) + SCSNR
        Write(IOut,1930) 'self-consistent solvent',Gen(41)
        Gen(55) = SCSNR
        endIf
C
C     Compute and print the VV10 constant energy shift.
      Call VV10bN(IOut,N6Func,NAE,NBE,ScaDFX(5),Gen(41))
C
C     Empirical dispersion term.
C
      Call GtCPar(IOut,IPrint,0,NAtoms,IAn,IAtTyp,MaxAtm,IAtPrM,AtParM,
     $  IV,NIAtP,NRAtP,NDAtP,jIAtP,jRAtP,jV,V,MDV)
      Call R6DIni(IOut,IOp(172),IOp(173),IOp(174),IOp(175),IOp(176),
     $  IOp(177),IOp(178))
      Call TStamp(1,'Bef disp')
      Call R6Disp(IOut,IPrint,'Nuclear repulsion',1,IEmpDi,IHMeth,
     $  IExCor,IXCFnc,ScaHFX,2,NAtoms,NPDir,IAtTyp,IAn,1,C,NDAtP,
     $  V(jIAtP),V(jRAtP),MxBond,V(INBond),V(IIBond),V(IRBTyp),Gen(41),
     $  XX,XX,V(jV),MDV-jV+1)
      Call TStamp(1,'Aft disp')
C
C     XDM initialization.
C
      Call XDMIni(IOut,IPrint,IOp(170),IOp(171),IHMeth,ScaHFX,IExCor,
     $  IXCFnc,IBasis,NGIC,IPolar,NAtoms,IAn,IAtTyp,V(IV),MDV1)
C
C     Now generate or read in the density basis.
C
      If(HavBas(IOut,IPrint,2,IBasis,NGIC,IDenBs,NGICDB)) then
        Call MolBas(In,IOut,IPrint,IDump,2,IPtAll,NewSE,IDenBs,NGICDB,
     $    IPolDB,IPurDB,IPurDB,RdExDB,AppNDB,LoadRN,IDBNrm,INucDB,
     $    LMxEqD,ICore,ISlatr,ISETyp,ICharg,Multip,NE,NAtoms,IAn,IAtTyp,
     $    AtmChg,C,RNFroz,AtParM(1,MAtNRd),IHMeth,IPFlag,AllowP,MaxASE,
     $    NSEAtI,V(jISEP),NSEAtR,V(jRSEP),MaxEl,LMaxEl,MM,MMDV,V(IV),
     $    MDV1)
        If(DoMaDB) Call Massge(IOut,IDump,1,.False.,MM,MMDF,NPDir)
        INrmCC = 2
        If(IDenBs.eq.29) INrmCC = 0
        IPtFlg = 100*IPtNDF + 20
        Call PostPB(In,IOut,IPrint,IDump,2,RDupDB,LMxEqD,LMaxDF,.False.,
     $    Strict,DelSml,Cut,.False.,.False.,AppNDB,UnCoDB,DoSplt,0,
     $    IEFBas,IOFInf,IRwGen,NAtoms,IAn,IAtTyp,AtmChg,C,NDBF,NDBPrm,
     $    IXYZ,MaxLM,NFld,RMag,EField,IDBNrm,INrmCC,IPtGIn,IPtGTb,
     $    IPtFlg,0,NoSym,NOpAll,NOp2,ISNEO,Gen,V(IV),MDV1)
      else
        IDenBs = -1
        endIf
C
C     Set global SCRF flag and initialize PCM. The 6th and 7th digits are
C     to be preserved from being reset by the GetPCM call since they're
C     set by Link1 during the various sub-calculations of a ONIOM-PCM run.
C
      IDS = -4
      NameS = ' '
      NS(1) = ' '
      LNS = 0
      LIS = 0
      LRS = 0
      If(IPrint.ge.4) Call SolvDB(IOut,IDS,NameS,LNS,NS,LIS,JJ,LRS,XX)
      IRdPCM = Mod(ISCRF,10)
      ISCRF4 = Mod(ISCRF,100000)/10000
      ISCRF = ISCRF - IRdPCM - ISCRF4*10000
      ISCRF6 = Mod(ISCRF,10000000)/1000000
      ISCRF7 = Mod(ISCRF,100000000)/10000000
      If(IRdPCM.eq.2.or.IRdPCM.eq.4.or.IRdPCM.eq.5) then
        UseChk = IRdPCM.eq.2.or.IRdPCM.eq.4
        Call GetPCM(IOut,IPrint,UseChk,0,1,ISCRF,RJunk)
        endIf
      If(PCM(ISCRF)) then
        IDervU = IDeriv
        Call ILSW(2,24,IDSet)
        If(IDSet.eq.-1) then
          IDervU = 0
        else if(IDSet.ge.1) then
          IDervU = Min(IDeriv,IDSet)
          endIf
        ISCRF = Mod(ISCRF,1000000)
        IDoInp = ISCRF6
        IDoCav = ISCRF7
        If(IDoCav.eq.0) IDoCav = 3
        DoElCv = IDoCav.eq.1.or.IDoCav.eq.3
        DoNElC = IDoCav.eq.2.or.IDoCav.eq.3
        IPCSol = IOp(140)
        IPCAcc = IOp(141)
        If(IPCAcc.eq.0.and.IAccFD.ne.0) IPCAcc = IAccFD
        IPCQC = IOp(142)
        IPCMxI = IOp(143)
        IPCThr = IOp(144)
        IPCPDn = IOp(166)
        IOptAA = 0
        If(IDoInp.ne.1) Call PCMInp(In,IOut,IPrint,ICharg,NAtoms,ISCRF,
     $    IRdPCM,IDervU,IOp(72),IOptAA,ISCRF4,ISCRF6,DoElCv,IPCSol,
     $    IPCAcc,IPCQC,IPCMxI,IPCThr,IPCPDn,IDfPCM,IAn,IAtTyp,C,V(IV),
     $    MDV1)
        Call PCMIOp('load',' ',IJunk,RJunk)
        NOpPCM = NOpAll
        If(DoElCv) Call ElCav(IOut,IPrint,NoSym,NOpPCM,NPDir,0,0,CelVec,
     $    XX,XX,PhyCon(1),NAtoms,IAn,C,V(IV),MDV1)
        If(DoNElC) then
          Call NECSet(IOut,IPrint,ICharg,NAtoms,IAn,IAtTyp,C,.True.,JJ,
     $      XX,XX,XX,JJ,XX,XX,XX,V(IV),MDV1)
          IDrvNE = 0
          If(IOp(71).lt.0) IDrvNE = IDervU
          Call NElCav(IOut,IPrint,IDrvNE,NoSym,NOpPCM,NPDir,0,0,CelVec,
     $      XX,XX,NAtoms,IAn,IAtTyp,C,.True.,JJ,XX,XX,XX,XX,XX,XX,JJ,XX,
     $      JJ,XX,JJ,XX,JJ,XX,JJ,XX,JJ,XX,V(IV),MDV1)
          GNEl = PCMGNE(IOut,IPrint)
          If(GNEl.ne.Zero) then
            If(IPSave(0).eq.0) Write(IOut,2140) GNEl
            Gen(41) = Gen(41) + GNEl
            Write(IOut,2150) Gen(41)
            endIf
          endIf
C
C       If there are external charges try to check whether they're
C       actually contained in the cavity.
C
        If(NumChg.gt.0) then
          jQChg = IV
          jCChg = jQChg + NumChg
          jEnd  = jCChg + 3*NumChg - 1
          Call TstCor(jEnd,MDV,'GInput-GausCk')
          jV = jEnd + 1
          MDV2 = MDV - jV + 1
          Call CpyChg(LChg,NumChg,V(jChg),V(jQChg),V(jCChg))
          Call GausCk(IOut,IPrint,NumChg,V(jQChg),V(jCChg),IPFlag,
     $      AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,NAtoms,IAn,IAtTyp,
     $      AtmChg,C,V(jV),MDV2)
          endIf
        Call PCMIOp('save',' ',IJunk,RJunk)
      else if(ISCRF.ne.0) then
        Call SCRFIn(IOut,IPrint,ISCRF,IRdPCM,V(IV),MDV1)
        endIf
      IDS = -5
      If(IPrint.ge.4) Call SolvDB(IOut,IDS,NameS,LNS,NS,LIS,JJ,LRS,XX)
      IStart = 1
      IPrt1 = Max(IPrint,1)
C
C     Pair potential.
C
      Call RdPrPt(IOut,IPrint,IRwPPT,0,NAtoms,IStart,NPot,IHead,IITPot,
     $  IATPot,Next,V,MDV)
      If(NPot.gt.0) then
        Call PairPt(.True.,0,NAtoms,1,C,NPot,V(IITPot),V(IATPot),EPPT,
     $    XX,XX)
        Write(IOut,2170) EPPT
        Gen(41) = Gen(41) + EPPT
        Write(IOut,2180) Gen(41)
        endIf
      Call FileIO(1,-IRwGen,LRwGen,Gen,0)
      If(IPrtBl.eq.1) Call BlkLst(IOut,IPrt1,IBkTyp,MinBkS,0,NAtoms,
     $  NBasis,NShell,I5DB1,I7FB1,JAn,ShellT,ShellC,IStart,NBlk,IstBfB,
     $  LstBfn,NTBlk,LenBT,IstBlT,LstBlk,Next,V,MDV)
      Call ILSW(1,69,ITyADD)
C
C     Update dimensioning information for fragments.
      If(UpdFrg) then
        IFrgTU = Max(IFrgTp,1)
        Call GetFrg(IOut,IPrint,0,IFrgTU,.False.,NAtoms,NBasis,IAtFrg,1,
     $    NFrag,IHavMO,NACore,NAVal2,NAVal1,NBCore,NBVal2,NBVal1,IFrStp,
     $    NDNBsC,NumNBC,IIChOn,IMulOn,INBCom,IMpF2M,IMpM2F,IIndFr,
     $    IIFrSt,IIFrAt,IChFrg,V,MDV)
        If(NFrag.eq.0)
     $    Call GauErr('Fragment setup requested but no fragments.')
        IEnd = IChFrg + NFrag - 1
        Call TstCor(IEnd,MDV,'GInput-RdWrFr')
        Call RdWrB(7,0)
        Call RdWrB(-11,0)
        If(IFrgTp.eq.0) then
          If(IFrgT1.eq.6.or.IFrgT1.eq.9) then
            IFrgTp = -1
          else
            IFrgTp = 1
            endIf
          endIf
        Call FrgDim(IOut,IPrint,IFrgTp,NAtoms,NBasis,NAE,NBE,IAn,IAtTyp,
     $    IAtFrg,AtmChg,AtZNuc,NFrag,V(IIChOn),V(IMulOn),NDNBsC,NumNBC,
     $    V(INBCom),V(IChFrg))
        Call RdWrFr(1,0,NAtoms,NFrag,IHavMO,NACore,NAVal2,NAVal1,NBCore,
     $    NBVal2,NBVal1,IFrStp,V(IIChOn),V(IMulOn),NDNBsC,V(INBCom))
        If(IPrint.ge.MinPFr) then
          Call GetFrg(IOut,IPrint,0,IAbs(IFrgTp),.True.,NAtoms,NBasis,
     $      IAtFrg,1,NFrag,IHavMO,NACore,NAVal2,NAVal1,NBCore,NBVal2,
     $      NBVal1,IFrStp,NDNBsC,NumNBC,IIChOn,IMulOn,INBCom,IMpF2M,
     $      IMpM2F,IIndFr,IIFrSt,IIFrAt,IChFrg,V,MDV)
          Call IMtOuS(IOut,'MapF2M',0,V(IMpF2M),NBasis,1,NBasis,1)
          Call IMtOuS(IOut,'MapM2F',0,V(IMpM2F),NBasis,2,NBasis,2)
          Call IMtOuS(IOut,'IndFr',0,V(IIndFr),NFrag+1,1,NFrag+1,1)
          Call PrPLis(IOut,'IFrgAt',NFrag,1,V(IIFrSt),V(IIFrAt))
          endIf
        endIf
C
C     BEBO parameters.
C
      Call RdBEBO(IOut,IPrint,IDump,MaxEl,IBEBO,IRwBEB,V,MDV)
C
C     Time-dependent external fields.
C
      Call RdPuls(In,IOut,IPrint,IEMPul)
C
C     Remove files which should be sorted out later.
      Call CkBadN(IOut,.True.,0,0)
      Call ChainX(0,V,MDV)
      Return
      End
