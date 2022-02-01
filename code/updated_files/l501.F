*Deck SimSCF
      Subroutine UHFOpn(V,MDV)
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
C     This link does a simple SCF and supports NEO SCF. 
C     It obeys the standard IOps for SCF links documented in l502.
C
#include "commonb2.inc"
      Parameter (LRwGen=1000,IRwGen=501,IRwH=515,IRwEig=522,IRwCA=524,
     $  IRwCB=526,IRwPA=528,IRwPB=530,IRwFA=536,IRwFB=538,IRwS=514,
     $  MOmega=DEFMOMEGA,NOmega=DEFNOMEGA,
     $  NScaDF=DEFNSCADF,MxGWIn=DEFMXGWIN,Two=2.0d0,F100=100.0d0)
Cxslis
      Parameter (IRwNEC=830,IRwNED=831,IRwNEE=843,IRwSGs=851)
Cxslie
      Dimension V(*), ScaDFX(NScaDF), DGen(LRwGen), MulTab(8,8), NS(9),
     $  Omega(MOmega,NOmega), CelVec(3,3), NS2(9), NS4(9)
      Integer FMFlag, FMFlg1, DirSc0, IGWInf(MxGWIn), ChrTbl(8,0:8)
      Logical DirSCF,FrcJEn,SvSort,DFT,AllowP(50),FMM,DoPrnt,DoPrn1,
     $  Halt,CBfn,Spinor,HarTst,OKSAO,No2E,ClShel,SimulO,DoDIIS
      Character*80 MethNm
      Integer MaxIOp, IOp
      Parameter (MaxIOp=DEFMAXIOP)
      Common /IOp/ IOp(MaxIOp)
#include "commonmol.inc"
#include "repall.inc"
      Common/IO/In,IOut,IPunch
 1110 Format(' SCF ',A,':  E(',A,') = ',G19.12,12H A.U. after ,I4,
     $  ' cycles',/,'            NFock=',I3,'  Conv=',D8.2)
 1130 Format(' E=',G22.15) 
C
      Call Drum(V,MDV)
C
C     IRwGen has various scalars, including the nuclear repulsion.
      Call FileIO(2,-IRwGen,LRwGen,DGen,0)
      ENR = DGen(41)
      jStrt = 1
      IExtp = Mod(IOp(11),10)
      IPrint = IOp(33)
      IFMM = IOp(63)
 
C     ISimul controls whether we want to do simultaneous optimization
C     for NEO. If ISimul != 1, do simultaneous  
C     IDIIS controls whether we want to use DIIS. If IDIIS != 1, do DIIS 
      ISimul = IOp(73)
      IDIIS = IOp(81)

C     SetSCF sets parameters used by most SCF links based on IOps
C     and information set up by overlays 3 and 4.
      Call SetSCF(IOut,IPrint,0,0,0,0,.True.,DoPrnt,DoPrn1,IOpCl,IROHF,
     $  NPDir,ISCRF,NBas6D,MaxCyc,FrcJEn,Halt,SvSort,OKSAO,IntTyp,
     $  IAccFD,IRaf,ISleaz,IAprXC,CnvDef,Acurcy,CelVec,NFull,NFullF,
     $  DirSc0,DirSCF,ICntrl,No2E,IIFock,JSym2E,HarTst,NOp1,NOp2)
      If(MaxCyc.eq.0) MaxCyc = 500
      IDenFit = 0
      IDoSPI = 0
C
C     SetUpF sets of parameters required for computation of Fock
C     matrices for any purpose (not just for SCF).  It allocates
C     arrays out of V starting at JStrt and returns jNext, the first
C     free location in V.
      Call SetUpF(IOut,IPrint,'SimSCF',IOpCl,IROHF,DirSc0,0,0,0,ICntrl,
     $  ICtDFT,1,0,No2E,SvSort,IDoSPI,IFMM,IOp(64),IOp(65),IRaf,IOp(60),
     $  IOp(61),IOp(62),IAccFD,0,ISCRF,Acurcy,CnvDef,OKSAO,IntTyp,
     $  0,JSym2E,NFullF,NAtoms,NBasis,NBas6D,IAn,IAtTyp,AtmChg,C,
     $  MaxAtm,IAtPrM,AtParM,NPDir,NOp1,NOp2,jStrt,MulTab,ChrTbl,IGWInf,
     $  MxGWIn,InStor,AccDes,NDBF6D,IPFlag,AllowP,FMM,FMFlag,FMFlg1,
     $  NFxFlg,IHMeth,Omega,IDenFit,NucTyp,DFT,IExCor,IXCFnc,N6Func,
     $  IRadAn,ScaHFX,ScaDFX,IRanWt,IRanGd,ICorTp,IEmpDi,NSAOBk,NOpUse,
     $  JSym2M,NBT,NPrDim,LenDup,CelVec,PCut,NS,NS2,NS4,ISEAll,LSEAll,
     $  IPMAll,LPMAll,IConct,NDAtP,jIAtP,jRAtP,jIPBC,LIPBC,IBuf2E,ICP,
     $  INEqB,INEqSh,INEqS2,IIJ,IIJMap,IOpClM,IR1,IR2,IR3,jNEqB,jISABF,
     $  jRSABF,jMap,IAtBtD,IIRBtD,IRRBtD,IGWBtD,jPMMA,LenXCU,jXCPU,
     $  jNext,V,MDV)
C
C     CkSCFT decodes IOpCl and checks whether it is consistent with the
C     forms of wavefunction implemented in this link.  NSpBlk is the
C     number of spin-blocks for the electronic wavefunction (1 for
C     closed-shell, 2 for UHF).  NSpBlX is the number of spin-blocks
C     during Fock formation (which is different from NSpBlk for GHF).
C     NSpBNX is the number of proton spin blocks (0 if not using NEOs,
C     1 for all spin-up protons, 2 for unrestricted proton
C     wavefunction).  NRI is 1/2 for real/complex wavefunction.
C     NDimBl is 1/2 for AO/spinorbital basis (1 for spin-aligned,
C     2 for GHF).
C
C     The following routines are called in CkSCFT
C     Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlY,Spinor)
C     NSpBlX = Mod(NSpBlY,10)
C     NSpBNX = Mod(NSpBlY,100)/10
      Call CkSCFT(IOut,'FockUt',IOpCl,.True.,.True.,.False.,.True.,
     $  NSpBlk,NRI,NDimBl,CBfn,NSpBlX,NSpBNX,Spinor)
C     Fetch a name for the method for printing.
      IROHF = 0
      Call NameSC(3,IOpCl,IROHF,IHMeth,IExCor,IXCFnc,IEmpDi,N6Func,
     $  ScaHFX,MethNm)  

C     NBD is the leading dimension of the MOs and the dimension of
C     the matrices being diagonalized.  There are NBsUse linearly
C     independent AO basis functions, which is also the number of MOs of
C     each spin for spin-aligned.  NBas6D is the number of cartesian
C     basis functions, which will be larger than NBasis for a basis
C     which uses pure angular momentum functions.  During Fock formation,
C     the densities are converted from pure to cartesian and the Fock
C     matrices are formed as cartesian and converted back, so these
C     arrays must be allocated large enough for the expanded version.
C     For GHF, the one complex lower triangle of dimension NBasis*2
C     is converted to 8 real lower triangles, 4 in DenA and 4 in DenB,
C     so these must also be large enough for this storage.
      NBD = NBasis*NDimBl
      NTT = (NBasis*(NBasis+1))/2
      NTT6D = (NBas6D*(NBas6D+1))/2
      NTTD = (NBD*(NBD+1))/2
C     NTTDR is the number of real values in a lower triangular
C     matrix.
      NTTDR = NRI*NTTD
      NBDSq = NBD*NBD
      NBDSqR = NRI*NBDSq
      NBasUD = NBsUse*NDimBl
      NOrbs = NBasUD
C
C     RdNERw returns the number of quantum protons, NEO basis functions
C     and offsets into the file holding the NEO 1e integrals.
      Call RdNERw(IOut,IOpCl,0,IRwNEX,NQNuc,NBsNEO,NBsNE6,NBUNEO,
     $  IOfNES,IOfNET,IOfNEH,IOfNER,IOpClE,IOClNP,IOClND)
C 
C     If this is a NEO calculation, get the dimensions      
      If(NQNuc.ge.0) then
C       Make sure some dimensions passed to subroutines are non-zero.
        NBsNEO = Max(NBsNEO,1)
        NBsNE6 = Max(NBsNE6,1)
C       Calculate length of the square matrix of NEO basis
        NEOSq = NBsNEO*NBsNEO
        NEOSqR = NRI*NEOSq
C       Calculate the size of lower triangular matrices to use in SCF      
        NEONTT = (NBsNEO*(NBsNEO+1))/2
        NEONT6 = (NBsNE6*(NBsNE6+1))/2
        endif

C     Allocate space in V for MOs, eigenvalues and Fock matrices,
C     and scratch arrays.  For convenience, we allocate the MO arrays as
C     full square matrices even if NBsUse<NBasis and allocate beta copies
C     even for closed-shell or GHF.
      LenCMO = NRI*NBD*NOrbs
      LDenFk = Max(NTTDR,NRI*NDimBl*NTT6D)
      IEigE = jNext
      IEigN = IEigE + NOrbs*2
      ICMOA = IEigN + NBsNEO*2
      ICMOB = ICMOA + NBDSqR
      IDenA = ICMOB + NBDSqR
      IDenB = IDenA + LDenFk
      IFockA = IDenB + LDenFk
      IFockB = IFockA + LDenFk
      ICMONA = IFockB + LDenFk
      ICMONB = ICMONA + NRI*NBsNEO*NBsNEO
      IDenNA = ICMONB + NRI*NBsNEO*NBsNEO
      If(NQNuc.gt.0) then
        LDFNuc = NRI*NEONT6
      else
        LDFNuc = 0
        endIf
      IDenNB = IDenNA + LDFNuc
      IFokNA = IDenNB + LDFNuc
      IFokNB = IFokNA + LDFNuc
      IScr1 = IFokNB + LDFNuc
C     The scratch space should be the larger basis function squared
      IV = IScr1 + Max(NEOSqR,Max(NBDSqR,NSpBlk*NTTDR))
      
C     Check that things fit.  MDV1 is the length free start at V(IV).
      Call TstCor(IV,MDV,'SimSCF')
      MDV1 = MDV - IV + 1
C
C     Read in the initial MOs.
      Call FileIO(2,-IRwCA,LenCMO,V(ICMOA),0)
      If(NSpBlk.eq.2) Call FileIO(2,-IRwCB,LenCMO,V(ICMOB),0)
Cxslis
      If(NQNuc.gt.0) then
        Call FileIO(2,-IRwNED,NBsNEO*(NBsNEO+1)/2,V(IDenNA),0)
        CAll AClear(LDFNuc,V(IDenNB))
        Call LTOutX(IOut,0,IOClNP,1,0,'Guess NEO Density Matrix',
     $    NBsNEO,NBsNEO,V(IDenNA),V(IDenNA),0)
        Call FileIO(2,-IRwNEC,NBsNEO*NBsNEO,V(ICMONA),0)
        Call AClear(NRI*NBsNEO*NBsNEO,V(ICMONB))
        Call LTOutX(IOut,3,IOCLNP,1,0,'Guess NEO coefficients',NBsNEO,
     $    NBsNEO,V(ICMONA),V(ICMONB),0)
        endIf

C     Set parameters for NEO SCF  
      MaxCyc = 20000
      ETol = 1d-10
      Acurcy = 1d-8
      AccDes = 1d-8
      AccMax = F100*Acurcy
      SimulO = .True.
      If(ISimul.eq.1) SimulO = .False.
      DoDIIS = .True.
      If(IDIIS.eq.1) DoDIIS = .False.
      MaxCom = 10
C
C     Iterate.
      If(NQNuc.eq.0) then
        Call SimCyc(IOut,IPrint,IOpCl,IExtp,NRI,NSpBlk,NSpBlX,NDimBl,
     $    NAtoms,NBas6D,NBasis,NBsUse,NBD,NTTD,NOrbs,ICharg,Multip,NAE,
     $    NBE,IAn,IAtTyp,AtmChg,C,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,
     $    IHMeth,IPFlag,AllowP,LSEAll,V(ISEAll),Omega,ICtDFT,
     $    IExCor,N6Func,ScaHFX,IGWInf,V(IAtBtD),V(IIRBtD),V(IRRBtD),
     $    V(IGWBtD),ScaDFX,ICorTp,NEqAll,V(IEigE),
     $    V(ICMOA),V(ICMOB),V(IDenA),V(IDenB),V(IFockA),
     $    V(IFockB),ENR,Energy,JCycle,KilSCF,MaxCyc,ETol,Acurcy,AccDes,
     $    AccMax,V(IScr1),V(IV),MDV1)
      else  
        Call NEOSCF(IOut,IPrint,IOpCl,IExtp,ICharg,Multip,IAn,IAtTyp,
     $    AtmChg,C,NAE,NBE,NAEX,NRI,NSpBlk,NSpBlX,NSpBNX,NDimBl,NAtoms,
     $    ENR,NBasis,NBD,NBas6D,NBsUse,NOrbs,NTTD,NTTDR,NQNuc,NBsNEO,
     $    NBsNE6,NBUNEO,NEONTT,NEONT6,IHMeth,IPFlag,AllowP,LSEAll,
     $    V(ISEAll),Omega,ICtDFT,IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,
     $    IGWInf,V(IAtBtD),V(IIRBtD),V(IRRBtD),V(IGWBtD),NEqAll,IRwFA,
     $    IRwFB,IRwS,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,IRwNED,IRwNEC,
     $    V(ICMOA),V(ICMOB),V(IDenA),V(IDenB),V(IFockA),V(IFockB),
     $    V(ICMONA),V(ICMONB),V(IDenNA),V(IDenNB),V(IFokNA),V(IFokNB),
     $    V(IEigE),V(IEigN),Energy,JCycle,KilSCF,MaxCyc,ETol,Acurcy,
     $    AccDes,AccMax,SimulO,DoDIIS,MaxCom,V(IScr1),V(IV),MDV1)
        endif

      
C     Print out the result of the SCF      
      If(KilSCF.eq.0) then
        Write(*,*) 'NEO-SCF converged after iteration #',JCycle
        Write(IOut,1130) Energy
      else If(KilSCF.eq.1) then
        Write(*,*) 'Stepwise SCF Failed. Electronic sub-system not 
     $    converged.'
      else If(KilSCF.eq.2) then
        Write(*,*) 'Stepwise SCF Failed. Protonic sub-system not 
     $    converged.'
      else If(KilSCF.eq.3) then
        Write(*,*) 'Stepwise SCF Failed. Total Energy of the system not 
     $    converged.'
      else If(KilSCF.eq.4) then
        Write(*,*) 'Simultaneous SCF Failed. Total E of the system not 
     $    converged.'
        Write(*,*) 'Conventional SCF Failed. Total E of the system not 
     $    converged.'
        endIf

C     Save final values.  Note that for silly historical reasons the
C     density stored for closed-shell wavefunctions in IRwPA is the
C     total rather than the alpha spin density.
      DGen(32) = Energy
      DGen(43) = Energy
      Call FileIO(1,-IRwGen,LRwGen,DGen,0)
      Call FileIO(1,-IRwCA,LenCMO,V(ICMOA),0)
      If(ClShel(IOpCl)) then
        Call AMove(NOrbs,V(IEigE),V(IEigE+NOrbs))
        Call AScale(NTTDR,Two,V(IDenA),V(IDenA))
        endIf
      Call ConWrt(IRwEig,2*NOrbs,V(IEigE))
      Call FileIO(1,-IRwPA,NTTDR,V(IDenA),0)
      If(NSpBlk.eq.2) then
        Call FileIO(1,-IRwCB,LenCMO,V(ICMOB),0)
        Call FileIO(1,-IRwPB,NTTDR,V(IDenB),0)
        endIf
      If(Kill.ne.0) Call GauErr('SCF failed to converge.')
      Call ChainX(0,V,MDV)
      Return
      End



*Deck NEOSCF
      Subroutine NEOSCF(IOut,IPrint,IOpCl,IExtp,ICharg,Multip,IAn,
     $  IAtTyp,AtmChg,C,NAE,NBE,NAEX,NRI,NSpBlk,NSpBlX,NSpBNX,NDimBl,
     $  NAtoms,ENR,NBasis,NBD,NBas6D,NBsUse,NOrbs,NTTD,NTTDR,NQNuc,
     $  NBsNEO,NBsNE6,NBUNEO,NEONTT,NEONT6,IHMeth,IPFlag,AllowP,LSEAll,
     $  RSEAll,Omega,ICtDFT,IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,IGWInf,
     $  IAtBtD,IRdBtD,RRdBtD,RGWBtD,NEqAll,IRwFA,IRwFB,
     $  IRwS,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,IRwNED,IRwNEC,CMOA,CMOB,
     $  DenA,DenB,FockA,FockB,CMONA,CMONB,DenNA,DenNB,FockNA,FockNB,
     $  EigE,EigN,Energy,JCycle,KilSCF,MaxCyc,ETol,Acurcy,AccDes,AccMax,
     $  SimulO,DoDIIS,MaxCom,Scr1,V,MDV)
      Implicit Real*8(A-H,O-Z)

      Parameter (Zero=0.0d0,One=1.0d0,SmalEn=1.d-9,MinPrt=2,
     $  F100=100.0d0)
      Logical AllowP(*), DoneQ, Rises, Damp, Done, SimulO, EDone, PDone,
     $  DoDIIS 
      Character*1 JCG
      Dimension IAn(*), IAtTyp(*), AtmChg(*), C(*), RSEAll(*), Omega(*),
     $  ScaDFX(*), IGWInf(*), IAtBtD(*), IRdBtD(*), RRdBtD(*),
     $  RGWBtD(*), NEqAll(*),EigE(NOrbs,*),EigN(NBsNEO,*),
     $  IndFok(MaxCom),EnCyc(MaxCom),ErrCyc(MaxCom),CMOA(NRI,NBD,NBD),
     $  CMOB(NRI,NBD,NBD), DenA(NRI,NTTD), DenB(NRI,NTTD),
     $  FockA(NRI,NTTD), FockB(NRI,NTTD), CMONA(NRI,NBsNEO,*),
     $  CMONB(NRI,NBsNEO,*), DenNA(NRI,NEONT6), DenNB(NRI,NEONT6),
     $  FockNA(NRI,NEONT6), FockNB(NRI,NEONT6), Scr1(NBD,NBD),
     $  V(MDV), JJ(1), XX(1),ESubV(10),
     $  ECDFit(1), EOne(1),
     $  ETwo(1), EX(1), EC(1), E6(1), ERFP(1), EPMM(1), RNE(1)
C     JJ and XX are dummies for arguments which are unused in the
C     calls made here.
      Save JJ, XX
      Data JJ/0/, XX/0.0d0/
 1120 Format(' E=',G22.15,' Delta-E=',F22.12,' Rises=',L1,' Damp=',L1,
     $  ' RMSDP=',1PD8.2,' MaxDP=',1PD8.2)
 1130 Format(' E=',G22.15)
 1170 Format(' E=','   (Non-Variational)  ')
C     NAEX is the number of occupieds for the first or only spin-block.
C     It is NAE for spin-aligned or the total number of electrons for GHF.
      NAEX = NAE + (NDimBl-1)*NBE
      NTTD = (NBD*(NBD+1))/2
      NTTDR = NRI*NTTD
 
C     Get Information about NEO      
      Call RdNERw(IOut,IOpCl,0,IRwNEX,NQNuc,NBsNEO,NBsNE6,NBUNEO,
     $  IOfNES,IOfNET,IOfNEH,IOfNER,IOpClE,IOClNP,IOClND)
      NEONTT = (NBsNEO*(NBsNEO+1))/2
      NEONT6 = (NBsNE6*(NBsNE6+1))/2

C     Allocate matrices for S and V, will be used in SCFDIS
      ISMat = 1
      IVMat = ISMat + NBasis*NBasis
      ISNMat = IVMat + NBasis*NBasis
      IVNMat = ISNMat + NBsNEO*NBsNEO

C     Allocate scratch spaces to call formV subroutine       
      IScr1 = IVNMat + NBsNEO*NBsNEO
      LenB62 = Max(NBasis*NBasis,NBsNEO*NBsNEO)
      LenB6 = Max(NBasis,NBsNEO)
      IScr2 = IScr1 + LenB62
      IScr3 = IScr2 + LenB6
      IScr4 = IScr3 + LenB6
      IBMat = IScr4 + Max(NBasis*NBasis,NBsNEO*NBsNEO) 

C     Allocate scratch matrices to be used in SCFDIS      
      IScrCof = IBMat + (MaxCom+1)**2
      IV = IScrCof + 7*(MaxCom+1)

C     Modify the length of V and MDV for later subroutines
      Call TstCor(IV,MDV,'NEOSCF')
      MDV1 = MDV - IV + 1
 
C     Initialize pointers to files used in SCFDIS
      IOBDIS = 0                                                        
      IOPAE = 0                                                        
      IOPBE = 0                                                        
      IOFAE = 0                                                        
      IOFBE = 0                                                        
      IOErrE = 0 
      IOPNA = 0                                                        
      IOFNA = 0                                                        
      IOErrN = 0 
      NFkSav = 0
     
C     Create temporary RW files to store S and V  
      IRwV = IGetRW(NBasis*NBasis,0)
      IRwNEF = IGetRw(NBsNEO*NBsNEO,0)
      IRwNES = IGetRw(NBsNEO*NBsNEO,0)
      IRwNEV = IGetRw(NBsNEO*NBsNEO,0)
      
C     Form/Check the S and V matrices
C     For electronic subsystem
      Call FileIO(2,-IRwS,NTTDR,V(ISMat),0)
      Call LTOutN(IOut,'Elec S Matrixx:',0,1,NBasis,V(ISMat),0)
      Call FileIO(1,-IRwS,NTTDR,V(ISMAT),0)
      Call FormV(1,NBD,NBsUse,IRwS,IRwV,V(IVMat),V(IScr1),V(IScr2),
     $  V(IScr3))
      Call LTOutN(IOut,'Elec V Matrixx:',0,1,NBasis,V(IVMat),0)
      Call XVFVMu(NRI,.False.,.False.,.False.,0,NBD,NOrbs,NBD,V(ISMat),
     $  V(IVMat),V(IScr1),V(ISMat),Scr1,NBD*NBD)
      Call LTOutN(IOut,'Elec V+SV',0,1,NBasis,V(ISMat),0)
    
C     For protonic subsystem
      Call FileIO(2,-IRwNEX,NEONTT,V(ISNMat),IOfNES)
      Call LTOutN(IOut,'Prot S Matrixx:',0,1,NBsNEO,V(ISNMat),0)
      Call FileIO(1,-IRwNES,NEONTT,V(ISNMAT),0)

C     Read protonic V from IRWNEX
      Call FileIO(2,-IRwNEX,NEONTT,V(IVNMat),IOfNER)
      Call LTOutN(IOut,'Prot V Matrixx:',0,1,NBsNEO,V(IVNMat),0)
      Call FileIO(1,-IRwNEV,NEONTT,V(IVNMAT),0)
      Call XVFVMu(NRI,.False.,.False.,.False.,0,NBsNEO,NBsNEO,NBsNEO,
     $  V(ISNMat),V(IVNMat),V(IScr1),V(ISNMat),Scr1,NBsNEO*NBsNEO)
      Call LTOutN(IOut,'Prot V+SV',0,1,NBsNEO,V(ISNMat),0)
      Call FileIO(2,-IRwNEX,NEONTT,V(ISNMAT),IOfNES)
      Call FileIO(1,-IRwNES,NEONTT,V(ISNMAT),0)

      Call AClear(LenB62,V(IScr1))
      Call AClear(LenB6,V(IScr2))
      Call AClear(LenB6,V(IScr3))
      Call AClear(LenB62,V(IScr4))
C     Form protonic V myself
      
      Call FileIO(2,-IRwNES,NEONTT,V(ISNMAT),0)
      Call LTOutN(IOut,'Prot S Before formV:',0,1,NBsNEO,V(ISNMat),0)
      Call FormV(1,NBsNEO,NBsNEO,IRWNES,IRwNEV,V(IVNMat),V(IScr1),
     $  V(IScr2),V(IScr3))


      Call LTOutN(IOut,'MyProt V LT:',0,1,NBsNEO,V(IVNMat),0)
      Call OutMtS(IOut,'MyProt V SQ',0,0,V(IVNMat),NBsNEO,NBsNEO,NBsNEO
     $  ,NBsNEO)
      Call XVFVMu(NRI,.False.,.False.,.False.,0,NBsNEO,NBUNEO,NBsNEO,
     $  V(ISNMat),V(IVNMat),V(IScr1),V(ISNMat),Scr1,NBsNEO*NBsNEO)
      Call LTOutN(IOut,'MyProt V+SV',0,1,NBsNEO,V(ISNMat),0)
      
      Call FileIO(2,-IRwNEV,NBsNEO*NBsNEO,V(IVNMat),0)
      
      Call LTOutN(IOut,'MyProt V LT:',0,1,NBsNEO,V(IVNMat),0)
      Call OutMtS(IOut,'MyProt V SQ',0,0,V(IVNMat),NBsNEO,NBsNEO,NBsNEO
     $  ,NBsNEO)


C     Print out the kinetic and potential integral for protonic
      Call FileIO(2,-IRwNEX,NEONTT,FockNA,IOfNET)
      Call LTOutN(IOut,'TN from l501:',0,1,NBsNEO,FockNA,0)
      Call FileIO(2,-IRwNEX,NEONTT,V,IOfNEH)
      Call ASub(NEONTT,V,FockNA,FockNA)
      Call LTOutN(IOut,'VN from l501:',0,1,NBsNEO,FockNA,0)
      
      LenCMO = NRI*NBD*NOrbs

C     Define Variables used in the SCF Procedure:
C     Old and current total energy for the system
      ETot = Zero
      EOTot = Zero
      DETot = Zero
C     Flags, counters for MacroSCF
      IFlagO = 0
      IFCyc = 0
C     These next two values are used in ConOpn to keep track of previous
C     density matrices
      IECnt = 0
      IPCnt = 0
C     These next two values are the number of iterations it took for
C     each subsystem to converge
      IEIter = 0
      IPIter = 0
 
C     SCF Settings:
      Damp = .False.
      DampFc = One
      KilSCF = 0
      MaxCom = 10
C     For Simul Opt
      EDone = .False.
      PDone = .False.
 
C     Creat Rw files for convergence check in Simultaneous Opt
      If (SimulO) then
        LScrE = NSpBlk*NTTDR
        LScrP = NEONTT
        IRwCE1 = IGetRW(LScrE,0)
        IRwCE2 = IGetRW(LScrE,0)
        IRwCE3 = IGetRW(LScrE,0)
        IRwCP1 = IGetRW(LScrP,0)
        IRwCP2 = IGetRW(LScrP,0)
        IRwCP3 = IGetRW(LScrP,0)
        endif

C     If we are doing a stepwise optimization, we need the first
C     protonic guess before entering the macro NEO-SCF iteration      
      If(.NOT. SimulO) then
        Call FormP2(.True.,NRI,NBsNEO,NQNuc,0,Zero,CMONA,DenNA)         
        Call FileIO(1,-IRwNED,NEONTT,DenNA,0)                           
        Call LTOutN(IOut,'Protonic Guess Density:',0,1,NBsNEO,DenNA,0)  
        endif

C     Top of Macro NEO-SCF Iteration Loop.
      Do 500 ITIter = 1, MaxCyc

        If(SimulO) then
C
C       SIMULTANEOUS OPTIMIZATION STARTS HERE:
        Call SMTSCF(IOut,IPrint,IOpCl,IExtp,ICharg,Multip,IAn,IAtTyp,
     $    AtmChg,C,NAE,NBE,NAEX,NRI,NSpBlk,NSpBlX,NSpBNX,NDimBl,NAtoms,
     $    NBasis,NBD,NBas6D,NBsUse,NOrbs,NTTD,NTTDR,NQNuc,NBsNEO,NBsNE6,
     $    NBUNEO,NEONTT,NEONT6,IHMeth,IPFlag,AllowP,LSEAll,RSEAll,Omega,
     $    ICtDFT,IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,IGWInf,IAtBtD,
     $    IRdBtD,RRdBtD,RGWBtD,NEqAll,IRwFA,IRwFB,IRwS,IRwH,IRwCA,IRwCB,
     $    IRwPA,IRwPB,IRwV,IRwNEX,IRwNED,IRwNEC,IRwNEF,IRwNES,IRwNEV,
     $    CMOA,CMOB,DenA,DenB,FockA,FockB,CMONA,CMONB,DenNA,DenNB,
     $    FockNA,FockNB,EigE,EigN,ESub,IECnt,IPCnt,EDone,PDone,Acurcy,
     $    AccDes,AccMax,DETot,ETol,ENR,IRwCE1,IRwCE2,IRwCE3,IRwCP1,
     $    IRwCP2,IRwCP3,MaxCom,IOBDIS,IOPAE,IOPBE,IOFAE,IOFBE,IOErrE,
     $    IOPNA,IOFNA,IOErrN,NFkSav,V(IBMat),V(IScrCof),EnCyc,ErrCyc,
     $    IndFok,Scr1,ITIter,DoDIIS,V(IV),MDV1)
        Write(*,*) "EDONE",EDone
        Write(*,*) "PDONE",PDone
        if(EDone.AND.PDone) then
          Write(*,*)'Electronic Density and Protonic Density
     $      Simultaneously Converged'
        else if(EDone.AND.(.Not. PDone)) then
          Write(*,*)'Electronic Density Converged, But the Protonic
     $      Density Has NOT Converged'
        else if(PDone.AND.(.Not. EDone)) then
          Write(*,*)'Protonic Density Converged, But the Electronic
     $      Density Has NOT Converged'
        else if( (.Not.PDone) .AND. (.Not. EDone)) then
          Write(*,*)'Both Electronic Density and Protonic Density Has
     $      NOT Converged'
          endif
C       SIMULTANEOUS OPTIMIZATION ENDS HERE        
        else  
C          
C       STEPWISE OPTIMIZATION STARTS HERE:
C         First we freeze the protonic density, and converge the
C         electronic density
C         Call SubSCF with a flag 0, indicating we want to converge
C         electronic density
          Call SubSCF(IOut,IPrint,IOpCl,IExtp,ICharg,Multip,IAn,IAtTyp,
     $      AtmChg,C,NAE,NBE,NAEX,NRI,NSpBlk,NSpBlX,NSpBNX,NDimBl,
     $      NAtoms,NBasis,NBD,NBas6D,NBsUse,NOrbs,NTTD,NTTDR,NQNuc,
     $      NBsNEO,NBsNE6,NBUNEO,NEONTT,NEONT6,IHMeth,IPFlag,AllowP,
     $      LSEAll,RSEAll,Omega,ICtDFT,IExCor,N6Func,ScaHFX,ScaDFX,
     $      ICorTp,IGWInf,IAtBtD,IRdBtD,RRdBtD,RGWBtD,NEqAll,IRwFA,
     $      IRwFB,IRwS,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,IRwV,IRwNEX,IRwNED,
     $      IRwNEC,IRwNEF,IRwNES,IRWNEV,CMOA,CMOB,DenA,DenB,FockA,FockB,
     $      CMONA,CMONB,DenNA,DenNB,FockNA,FockNB,EigE,EigN,IECnt,IPCnt,
     $      ESub,ISubCy,Kill,MaxCyc,Acurcy,AccDes,AccMax,ETol,MaxCom,
     $      ENR,0,Scr1,DoDIIS,V(IV),MDV1)
          
          if (Kill.eq.1) then
            KilSCF = 1
            goto 800
            endif
          IEIter = ISubCy
          Write(*,*)'Electronic Density Converged!!!'
          Write(*,*)'The Converged Electronic Subsystem Energy is', ESub

C         Now freeze the converged electronic density, and converge the
C         protonic density
C         Call SubSCF with a flag 1, indicating we want to converge
C         protonic density      
          Call SubSCF(IOut,IPrint,IOpCl,IExtp,ICharg,Multip,IAn,IAtTyp,
     $      AtmChg,C,NAE,NBE,NAEX,NRI,NSpBlk,NSpBlX,NSpBNX,NDimBl,
     $      NAtoms,NBasis,NBD,NBas6D,NBsUse,NOrbs,NTTD,NTTDR,NQNuc,
     $      NBsNEO,NBsNE6,NBUNEO,NEONTT,NEONT6,IHMeth,IPFlag,AllowP,
     $      LSEAll,RSEAll,Omega,ICtDFT,IExCor,N6Func,ScaHFX,ScaDFX,
     $      ICorTp,IGWInf,IAtBtD,IRdBtD,RRdBtD,RGWBtD,NEqAll,IRwFA,
     $      IRwFB,IRwS,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,IRwV,IRwNEX,IRwNED,
     $      IRwNEC,IRwNEF,IRwNES,IRWNEV,CMOA,CMOB,DenA,DenB,FockA,FockB,
     $      CMONA,CMONB,DenNA,DenNB,FockNA,FockNB,EigE,EigN,IECnt,IPCnt,
     $      ESub,ISubCy,Kill,MaxCyc,Acurcy,AccDes,AccMax,ETol,MaxCom,
     $      ENR,1,Scr1,DoDIIS,V(IV),MDV1)

          if (Kill.eq.1) then
            KilSCF = 2
            goto 800
            endif
          IPIter = ISubCy
          Write(*,*) 'Protonic Density Converged!!! '
          Write(*,*)'The Converged Protonic Subsystem Energy is', ESub
C       STEPWISE OPTIMIZATION ENDS HERE 
          endif

C       Retrieve Data and Analyze:
        ETot = ESub
C       Printing out the Energy for This Macro Iteration           
        If (SimulO) then
C         For SMT, just print out the energy          
          Write(*,*) 'The Energy for this Simultaneous Optimization
     $      Iteration is: ', ETot
        else  
C         For STW, If we reach here, then both electronic and protonic densities
C         are stationary. We can print out a quick summary for this
C         macro iteration
          Write(*,*) 'Summary for Macro-SCF Iteration #', ITIter
          Write(*,*) 'Electronic Converged after Iteration #', IEIter
          Write(*,*) 'Protonic Converged after Iteration #', IPIter
          Write(*,*) 'NEO-SCF Energy at This Iteration Step', ETot
          endif

C       Use total energy of the system as the convergence criterion        
        DETot = ETot - EOTot
        Energy = ETot      
        

        If(Abs(DETot).lt.ETol) then
          If (SimulO .AND. ((.NOT.EDone) .Or. (.NOT.PDone))) then
            Write(*,*) 'Total Energy in Simultaneous Opt Has Converged,
     $        But the Subsystem Density Has NOT Converged.'
            goto 500
            endif

          Write(*,*) 'TOTAL ENERGY HAS CONVERGED!!!'
          goto 800
          endif

C       If not converged, we update the old energy and print it out.
C       Then go in the the next iteration of the outer loop

        Write(*,*) 'Total Energy Has NOT Converged!!!'

  500   EOTot = ETot

C     If we reach here, then SCF has failed
      If(SimulO) then
       KillSCF = 4
      else   
C       If we reach here, that means although the energies & densities for
C       both sub-systems have converged, the total energy of the system is
C       not converged, then we exit with flag "KilSCF = 3"
        KilSCF = 3 
        endif
C     End of the routine. Free scratch space and report # of cycles it
C     takes for the energy of the system to converge
  800 JCycle = ITIter 
      Return
      End

*Deck SubSCF
      Subroutine SubSCF(IOut,IPrint,IOpCl,IExtp,ICharg,Multip,IAn,
     $  IAtTyp,AtmChg,C,NAE,NBE,NAEX,NRI,NSpBlk,NSpBlX,NSpBNX,NDimBl,
     $  NAtoms,NBasis,NBD,NBas6D,NBsUse,NOrbs,NTTD,NTTDR,NQNuc,NBsNEO,
     $  NBsNE6,NBUNEO,NEONTT,NEONT6,IHMeth,IPFlag,AllowP,LSEAll,RSEAll,
     $  Omega,ICtDFT,IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,IGWInf,IAtBtD,
     $  IRdBtD,RRdBtD,RGWBtD,NEqAll,IRwFA,IRwFB,IRwS,IRwH,IRwCA,IRwCB,
     $  IRwPA,IRwPB,IRwV,IRwNEX,IRwNED,IRwNEC,IRwNEF,IRwNES,IRwNEV,CMOA,
     $  CMOB,DenA,DenB,FockA,FockB,CMONA,CMONB,DenNA,DenNB,FockNA,
     $  FockNB,EigE,EigN,IECnt,IPCnt,ESub,ISubCy,Kill,MaxCyc,Acurcy,
     $  AccDes,AccMax,ETol,MaxCom,ENR,IS2Con,Scr1,DoDIIS,V,MDV)
     
        Implicit Real*8(A-H,O-Z)
C
C     This subroutine performs the SCF procedure on either the
C     electronic or the protonic sub-system, while freezing the density
C     of the other sub-system
C

C     Input arguments:
C     IOut,IPrint,IOpCl,IExtp: have their usual meanings.

C     Molecule Information:
C     ICharg    ... The total electric charge on the molecule,                
C                   0 for neutral, 1 for cations, etc.
C     Multip    ... The molecule's spin multiplicity (1 for                   
C                   singlets, 2 for doublets, etc.)
C     IAN(I)    ... The atomic number of atom I.
C     IAtTyp(I) ... The type of atom I. 0 for normal atom, 2000 for NEO
C     AtmChg(I) ... The nuclear charge of atom I, accounting for any 
C                   effective core potential.
C     C(3,I)    ... The coordinates of the atoms in atomic units
C     NAE       ... The number of alpha-spin electrons
C     NBE       ... The number of beta-spin electrons
C     NAEX      ... The number of occupieds for the first or only spin-block.
C     It is NAE for spin-aligned or the total number of electrons for GHF.
C     NRI       ... 1/2 for real/complex wavefunction.
C     NSpBlk    ... The number of spin-blocks for the electronic
C                   wavefunction (1 for closed-shell, 2 for UHF).
C     NSpBlX    ... The number of spin-blocks during Fock formation (which
C                   is different from NSpBlk for GHF).
C     NSpBNX    ... The number of proton spin blocks (0 if not using NEOs,
C                   1 for all spin-up protons, 2 for unrestricted proton 
C                   wavefunction).
C     NDimBl    ... 1/2 for AO/spinorbital basis (1 for spin-aligned, 2 for
C                   GHF) 
C
C     Electronic Sub-system Information
C     NAtoms    ... The total number of atoms
C     NBasis    ... The total number of basis functions
C     NBD       ... The leading dimension of the MOs and the leading
C                   dimensionof the matrices being diagonalized (NBasis 
C                   * NDimBl).
C     NBas6D    ... The number of cartesian basis functions, which will be
C                   larger than NBasis for a basis which uses pure angular 
C                   momentum fucntions.
C     NBsUse    ... The total number of linearly independent AO basis
C                   function.
C     NOrbs     ... The number of orbitals (NBsUse * NDimBl).
C     NTTD      ... The size of the lower triangular matrix with
C                NBD
C     NTTDR     ... The number of real values stored in NTTD
C
C     Protonic Sub-system Information:
C     NQNuc  ... The number of quantum nuclei.
C     NBsNEO ... The total number of NEO basis functions.
C     NBsNE6 ... The number of cartesian NEO basis functions.
C     NBUNEO ... The number of linearly independent NEO basis functions
C     NEONTT ... The size of the lower triangular matrix with dimension
C                NBsNEO
C     NEONT6 ... The size of the lower triangular matrix with dimension
C                NBsNE6
C
C
C     Options:    
C     IHMeth ... Method:
C                0 ... Default (1).
C                1 ... Ab initio
C                2 ... CNDO/2
C                3 ... INDO/2
C                4 ... ZINDO/1 (NYI)
C                5 ... ZINDO/S
C                6 ... MINDO/3 (NYI)
C                7 ... MNDO
C                8 ... AM1
C                9 ... PM3
C               10 ... PM3MM
C            11,12 ... Harris functional.  All the terms for this
C                      are done in CalDSu, so this routine just
C                      initializes empty matrices.
C               13 ... Approximate DFT with Harris XC but regular
C                      density fitting J.
C               14 ... Approximate DFT with Harris XC and atomic
C                      densities for Coulomb, but done through
C                      regular DFT code.
C               15 ... Regular SCF with separate K, for testing,
C                      Treated as ab initio here.
C               16 ... Regular SCF with separate NDDO K,
C                      Treated as ab initio here.
C               17 ... K part of model 16, used internally.
C               18 ... DFT-SCTB using interpolated values.
C               19 ... DFT-SCTB using analytic values.
C               20 ... EHT-SC.
C               21 ... 2e terms forced to be zero.
C               22 ... Approximate DFT with Harris XC and DFTB-style
C                      charge iteration.
C               23 ... Approximate DFT with Harris XC and improved
C                      DFTB-style charge iteration.
C               24 ... PM6 with overlap, solving in AO basis.
C               25 ... PM6HXC.
C               26 ... PM6AXC.
C               27 ... NDBEBO.
C            28-39 ... Reserved for non-zdo semi-empirical models.
C               40 ... PM6.
C               41 ... PDDG/PM3.
C               42 ... PM6E/PM6PFD.
C               43 ... PM7.
C               44 ... PM6 with T transformed to OAO.
C               45 ... PM7TS.
C               46 ... PM7MOPAC.
C               47 ... Old (G09) PM6.
C            48-99 ... Reserved for zdo semi-empirical models.
C              101 ... Molecular Mechanics (not used here, but
C                      relevent for PCM in FoFCou).
C     IPFlag ... Flags for Prism, PrfRaf, CalDSu, etc.:
C                1 ... Use AllowP(50) to turn off paths.
C                2 ... If set, force front and back end memory allocations
C                      to match.  Otherwise, allocate independently.
C                4 ... Run a script to execute Prism, PrmRaf, or CalDFT.
C                8 ... Use unformatted I/O on the data for the script process.
C               16 ... Do extra work to use cutoffs better in DFTMat.
C               32 ... Reverse normal choice of diagonal/canonical sampling
C                      in Prism and PrmRaf.  The default is diagonal only
C                      on vector machines.
C               64 ... Trace input and output using Linda/subprocess.
C              128 ... Force single matrix code in CPKS.
C     AllowP ... Logical(50) array for PRISM if IPFlag is set.
C     ISEAll ... The usual array of semi-empirical parameters, which
C                has length LSEAll (From SetupF)
C     RSEALL ... Some semi-empirical related quantity
C     Omega  ... Flags and coefficients for short and long range Coulomb.
C     ICtDFT ... Control flag for DFT.IExCor,N6Func,ScaHFX are the usual XC quadrature settings.
C     IExCor,N6Func,ScaHFX,ScaDFX,ICorTp are the usual XC quadrature settings.
C     IGWInf ... Array for precomputed XC quadrature information.
C                MxGWIn should be its allocated length.
C     IAtBtD,IIRBtD,IRRBtD,IGWBtD are the usual arrays for precomputed
C       XC quadrature information.
C     IAtBtD(1,I) -- Offset of first radial shell for active atom I.
C     IAtBtD(2,I) -- Real atom number for active atom I.
C     IAtBtD(3,I) -- Offset of first elements in V(jGWBtD) for first
C                    grid point on I.
C     IAtBtD(4,I) -- Number of radial shells stored for atom I.  May be less
C                    than the number implied by IRadAn since radial shells which
C                    are too far from atom I to have any signficant shells are
C                    not stored.
C
C     IRdBtD(1,I) -- Number of angular points in radial shell I
C     IRdBtD(2,I) -- Number of stored angular points in radial shell I,
C                    after symmetry.  Currently either NAng or 0 if
C                    there are no significant shells on I.
C     IRdBtD(3,I) -- Number of significant basis functions on shell I, MaxDrB=0
C     IRdBtD(4,I) -- Offset of first elements in GWBtD for the first
C                    point in shell I.
C     IRdBtD(5,I) -- Number of points in radial shell I with significantly
C                    non-zero weights.
C     RRdBtD(1,I) -- Radius of shell I.
C     RRdBtD(2,I) -- Cost of shell I, stored by CalDFT.
C     RGWBtD      -- Sth DFT related I assume
C     NEqAll ... Atom permutations in symmetry operations.

C     Files:
C     IRwH   ... Contains electronic 1e integrals 
C     IRwCA  ... Contains MO Coeff. for alpha spin electrons
C     IRwCB  ... Contains MO Coeff. for beta spin electrons
C     IRwPA  ... Contains density matrix for alpha spin e-
C     IRwPB  ... Contains density matrix for beta spin e-
C     IRWNEX ... Contains NEO 1e integrals
C     IRwNED ... Contains NEO density matrix (right now all spin-up)
C     IRwNEC ... Contains NEO MO Coeff.
C
C
C     Output Arguments:
C     CMOA   ... Alpha-spin MO Coeff.
C     CMOB   ... Beta-spin MO Coeff.
C     DenA   ... Alpha-spin density matrix
C     DenB   ... Beta-spin density matrix
C     FockA  ... Alpha-spin fock matrix
C     FockB  ... Beta-spin fock matrix
C     CMONA  ... Alpha-spin NEO MO Coeff
C     CMONB  ... Beta-spin NEO MO Coeff (Not implemented)
C     DenNA  ... Alpha-spin NEO density matrix
C     DenNB  ... Beta-spin NEO density matrix (Not implemented)
C     FockNA ... Alpha-spin NEO fock matrix
C     FockNB ... Beta-spin NEO fock matrix (Not implemented)
C     ESub   ... Energy of the sus-system (as a scaler)
C     ISubCy ... The number of cycles it takes for sus-system to
C                converge
C     Kill   ... Flag to indicate the convergence of SCF
C
C     SCF options
C     MaxCyc ... Max number of iteration.
C     Acurcy ... Convergence threshold.  If this is non-zero and tighter
C                than CnvDef, the integral accuracy will be adjusted
C                accordingly.
C     AccDes ... Desired accuracy.
C     AccMax ... Maximum allowed change in density matrix elements     
C     ETol   ... Energy tolerance
C     
C     ENR    ... Electronic-Classical Nuclear Repulsion
C     IS2Con ... A label to specify which sub-system to converge.
C                If 0, we freeze NEO density and converge e- density
C                If 1, we freeze e- density and converge NEO density

      Parameter (Zero=0.0d0,One=1.0d0,SmalEn=1.d-9,MinPrt=2,
     $  F100=100.0d0)
      Logical AllowP(*), DoneQ, Rises, Damp, Done, DoDIIS, EDone, PDone
      Character*1 JCG
      Dimension IAn(*), IAtTyp(*), AtmChg(*), C(*), RSEAll(*), Omega(*),
     $  ScaDFX(*), IGWInf(*), IAtBtD(*), IRdBtD(*), RRdBtD(*),
     $  RGWBtD(*), NEqAll(*),EigE(NOrbs,*),EigN(NBsNEO,*),
     $  CMOA(NRI,NBD,NBD),
     $  CMOB(NRI,NBD,NBD), DenA(NRI,NTTD), DenB(NRI,NTTD),
     $  FockA(NRI,NTTD), FockB(NRI,NTTD), CMONA(NRI,NBsNEO,*),
     $  CMONB(NRI,NBsNEO,*), DenNA(NRI,NEONT6), DenNB(NRI,NEONT6),
     $  FockNA(NRI,NEONT6), FockNB(NRI,NEONT6),IndFok(MaxCom),
     $  ErrCyc(MaxCom),EnCyc(MaxCom),
     $  Scr1(NRI,Max(NBD,NBsNEO),Max(NBD,NBsNEO)),
     $  V(MDV), JJ(1), XX(1), ESubV(10),
     $  ECDFit(1), EOne(1),
     $  ETwo(1), EX(1), EC(1), E6(1), ERFP(1), EPMM(1), RNE(1)
C     JJ and XX are dummies for arguments which are unused in the
C     calls made here.
      Save JJ, XX
      Data JJ/0/, XX/0.0d0/
 1120 Format(' E=',G22.15,' Delta-E=',F22.12,' Rises=',L1,' Damp=',L1,
     $  ' RMSDP=',1PD8.2,' MaxDP=',1PD8.2)
 1130 Format(' E=',G22.15)
 1170 Format(' E=','   (Non-Variational)  ')
 1180 Format(' In Electronic SubSCF Loop, Total Energy = ', F22.15, 
     $ ' After Iteration # ', I4)     
 1190 Format(' In Protonic SubSCF Loop, Total Energy = ', F22.15, 
     $ ' After Iteration # ', I4)     
C
      Write(*,*)"Enter SubSCF"

      Call RdNERw(IOut,IOpCl,0,IRwNEX,NQNuc,NBsNEO,NBsNE6,NBUNEO,
     $  IOfNES,IOfNET,IOfNEH,IOfNER,IOpClE,IOClNP,IOClND)
   
C     Previous and current total energy for the sub-system
      ECur = Zero
      EPre = Zero
      DE = Zero

C     Flags, counters for SCF
      IFlagO = 0
      IFlag = 0
      ISubCy = 0
      Damp = .False.
      DampFc = One
      Kill = 0
      EDone = .False.
      PDone = .False.

      LenCMO = NRI*NBD*NOrbs
C     Scratch files used in extrapolation, again make sure it's large
C     enough for both electronic/protonic systems.
      LScr = Max(NBsNEO*(NBsNEO+1)/2,NSpBlk*NTTDR)
      IRwC1 = IGetRZ(LScr,0,V,MDV)
      IRwC2 = IGetRZ(LScr,0,V,MDV)
      IRwC3 = IGetRZ(LScr,0,V,MDV)
C
C     Allocate scratch matrices to be used in SCFDIS      
      MaxCom = 10
      IScrA = 1
      LenADim = Max(Max(NBasis,NBsNEO),MaxCom)
      IScrC = IScrA + LenADim*(LenADim+1)
      IScrD = IScrC + LenADim*LenADim
      LenDDim = Max(Max(NBasis,NBsNEO),(MaxCom+1)**2)
      IBMat = IScrD + LenDDim
      IScrCof = IBMat + (MaxCom+1)**2
      IV = IScrCof + 7*(MaxCom+1)
      Call TstCor(IV,MDV,'SubSCF')
      MDV1 = MDV - IV + 1
      IOBDIS = 0                                                        
      IOPADI = 0                                                        
      IOPBDI = 0                                                        
      IOFADI = 0                                                        
      IOFBDI = 0                                                        
      IOErDI = 0 
      NFESAV = 0
      NFPSAV = 0

C     Top of Sub-system SCF loop.

      Do 600 ISIter = 1, MaxCyc
          
        if(IS2Con.eq.0) then
C         Form the electronic density matrix and save it in IRwPA
          Call FormP2(.True.,NRI,NBD,NAEX,0,Zero,CMOA,DenA)
          If(NSpBlk.eq.2) Call FormP2(.True.,NRI,NBD,NBE,0,Zero,CMOB,
     $      DenB) 
          Call FileIO(1,-IRwPA,NTTDR,DenA,0)
          If(NSpBlk.eq.2) Call FileIO(1,-IRwPB,NTTDR,DenB,0)
          Call LTOutN(IOut,'New Electronic  Density',0,1,NBasis,DenA,0)
C         Check if the electronic density has converged
          Call ConOpn(IOut,IOpCl,NBasis,NBsUse,RMSDPE,IFlag,IExtp,DenA,
     $      FockA,Scr1,Acurcy,AccMax,Zero,DelEE,DelMax,IECnt,EDone,
     $      IRwPA,IRwPB,IRwC1,IRwC2,IRwC3,ETol,1,1,Damp,DampFc)

C         If the electronic density has converged, we will go the endpoint 
C         of this routine. But we need to go through one more SCF
C         iteration to get the energy from our final converged density
          If(EDone) then
            Call FileIO(2,-IRwPA,NTTDR,DenA,0)
            Call MkFok1(IOut,4,IOpCl,.False.,.False.,.False.,.False.,0,
     $        1,1,0,1,NAE,NBE,ICharg,Multip,NBasis,NBas6D,NBsUse,DenA,
     $        DenB,FockA,FockB,1,NOpUse,1,1,.True.,0,IPFlag,AllowP,
     $        .False.,0,0,0,0,0,IHMeth,LSEAll,RSEAll,Omega,NBTI,AccDes,
     $        AccXC,ICtDFT,IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,0,IAn,
     $        IAtTyp,AtmChg,C,NAtoms,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,
     $        .True.,0,1,0,0,.False.,NSaved,1,XX,0,0,.False.,1,JJ,JJ,JJ,
     $        JJ,XX,JJ,JJ,XX,XX,XX,JJ,0,XX,JCG,ESubV,EElec,
     $        ECDFit,EOne,ETwo,
     $        Ex,Ec,E6,ErfP,EPMM,DoneQ,XX,One,XX,JJ,JJ,JJ,NAtoms,NEqAll,
     $        XX,XX,RNE,JJ,JJ,IGWInf,IAtBtD,IRdBtD,RRdBtD,RGWBtD,0,XX,1,
     $        JJ,XX,XX,NBsNEO,DenNA,DenNB,FockNA,FockNB,XX,XX,0,XX,V(IV)
     $        ,MDV1)
            Call FileIO(2,-IRwNEX,NEONTT,FockNA,IOfNEH)
            Call FoFNEO(IOut,IPrint,IUHFH,.False.,.True.,.True.,NBasis,
     $        NBas6D,NBsNEO,NBsNE6,Omega,.False.,0,0,0,
     $        IPFlag,AllowP,ISym2E,NAtoms,IAn,AtmChg,C,IAtTyp,JJ,NOpUse,
     $        NOp2,NEqAll,XX,RotAll,XX,XX,
     $        XX,AccDes,DenA,DenNA,FockA,FockB,FockNA,ESubV,V(IV),MDV1)
            ESub = ENR + ESubV(1) + ESubV(2) + ESubV(3) + ESubV(4)
            goto 900
            endif

C         If electronic density is not converged, we need to form a new
C         form matrix using this density.
C         Note Here:        
C         In ConOpn,DenA will be holding P_{N-1}. Therefore we need to 
C         read the DenA from IRwPA again
          Call FileIO(2,-IRwPA,NTTDR,DenA,0)
        else 
C         Form the protonic density matrix
C
          Call FormP2(.True.,NRI,NBsNEO,NQNuc,0,Zero,CMONA,DenNA)
          Call FileIO(1,-IRwNED,NEONTT,DenNA,0)
          Call LTOutN(IOut,'Protonic Density:',0,1,NBsNEO,DenNA,0)

C         Check the convergence of the protonic density matrix
          Call ConOpn(IOut,IOpCl,NBsNEO,NBUNEO,RMSDPP,IFlag,IExtp,DenNA,
     $      FockNA,Scr1,Acurcy,AccMax,Zero,DE,DelMax,IPCnt,PDone,
     $      IRwNED,IRwNED,IRwC1,IRwC2,IRwC3,ETol,1,1,Damp,DampFc)
C         If protonic density has converged, we jump to the end of the
C         routine. But we still need to go through one last iteration of
C         SCF to get the energy of the system using the final converged
C         density
          If(PDone) then
            Write(*,*)'Energy Right Now'
            Write(*,*)ESub
            Write(*,*)'1ONE LAST ITER AFTER PN CONVERGED'
            Call FileIO(2,-IRwPA,NTTDR,DenA,0)
            Call MkFok1(IOut,4,IOpCl,.False.,.False.,.False.,.False.,0,
     $        1,1,0,1,NAE,NBE,ICharg,Multip,NBasis,NBas6D,NBsUse,DenA,
     $        DenB,FockA,FockB,1,NOpUse,1,1,.True.,0,IPFlag,AllowP,
     $        .False.,0,0,0,0,0,IHMeth,LSEAll,RSEAll,Omega,NBTI,AccDes,
     $        AccXC,ICtDFT,IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,0,IAn,
     $        IAtTyp,AtmChg,C,NAtoms,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,
     $        .True.,0,1,0,0,.False.,NSaved,1,XX,0,0,.False.,1,JJ,JJ,JJ,
     $        JJ,XX,JJ,JJ,XX,XX,XX,JJ,0,XX,JCG,ESubV,EElec,
     $        ECDFit,EOne,ETwo,
     $        Ex,Ec,E6,ErfP,EPMM,DoneQ,XX,One,XX,JJ,JJ,JJ,NAtoms,NEqAll,
     $        XX,XX,RNE,JJ,JJ,IGWInf,IAtBtD,IRdBtD,RRdBtD,RGWBtD,0,XX,1,
     $        JJ,XX,XX,NBsNEO,DenNA,DenNB,FockNA,FockNB,XX,XX,0,XX,V(IV)
     $        ,MDV1)
            Call FileIO(2,-IRwNEX,NEONTT,FockNA,IOfNEH)
            Call FoFNEO(IOut,IPrint,IUHFH,.False.,.True.,.True.,NBasis,
     $        NBas6D,NBsNEO,NBsNE6,Omega,.False.,0,0,0,
     $        IPFlag,AllowP,ISym2E,NAtoms,IAn,AtmChg,C,IAtTyp,JJ,NOpUse,
     $        NOp2,NEqAll,XX,RotAll,XX,XX,
     $        XX,AccDes,DenA,DenNA,FockA,FockB,FockNA,ESubV,V(IV),MDV1)
            ESub = ENR + ESubV(1) + ESubV(2) + ESubV(3) + ESubV(4)
            Write(*,*)'2ONE LAST ITER AFTER PE CONVERGED'
            Write(*,*)ESub
            goto 900
            endif
C         If not converged, we form a using fock matrix using the       
C         current protonic density and diagnolize                       
C         Note:                                                         
C         In ConOpn, DenNA will be holding the value of P_{N-1}, so we  
C         need to read DenNA from IRwNED again                          
          Call FileIO(2,-IRwNED,NEONTT,DenNA,0) 
          endif
        Call FileIO(2,-IRwPA,NTTDR,DenA,0) 
        Call LTOutN(IOut,'E- Density before MKFok1:',0,1,NBasis,DenA,0) 
        Call LTOutN(IOut,'E- Fock before MKFok1:',0,1,NBasis,FockA,0)
        NOpUse = 1
        NBTI = 0
        AccXC = AccDes
        NSaved = 0
        DoneQ = .False.

        Call MkFok1(IOut,4,IOpCl,.False.,.False.,.False.,.False.
     $    ,0,
     $    1,1,0,1,NAE,NBE,ICharg,Multip,NBasis,NBas6D,NBsUse,DenA,
     $    DenB,
     $    FockA,FockB,1,NOpUse,1,1,.True.,0,IPFlag,AllowP,.False.
     $    ,0,0,0,
     $    0,0,IHMeth,LSEAll,RSEAll,Omega,NBTI,AccDes,AccXC,ICtDFT,
     $    IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,0,IAn,IAtTyp,AtmChg,C,
     $    NAtoms,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,.True.,0,1,0,0,.False.,
     $    NSaved,1,XX,0,0,.False.,1,JJ,JJ,JJ,JJ,XX,JJ,JJ,XX,XX,XX,JJ,
     $    0,
     $    XX,JCG,ESubV,EElec,ECDFit,EOne,ETwo,Ex,Ec,E6,ErfP,
     $    EPMM,DoneQ,XX,
     $    One,XX,JJ,JJ,JJ,NAtoms,NEqAll,XX,XX,RNE,JJ,JJ,IGWInf,IAtBtD,
     $    IRdBtD,RRdBtD,RGWBtD,0,XX,1,JJ,XX,XX,NBsNEO,DenNA,DenNB,
     $    FockNA,FockNB,XX,XX,0,XX,V(IV),MDV1)
        Write(*,*)'After MkFok1' 
        Call LTOutN(IOut,'Fee from MkFok1:',0,1,NBasis,FockA,0)
        Write(*,*)IOfNEH 
        Call FileIO(2,-IRwNEX,NEONTT,FockNA,IOfNEH)
 
C         TODO: Symmetry is turned off 
          ISym2E = 0
          NOpUse = 0
          NOp2 = 0
C         TODO: RHF only
          IUHFH = 0

C         Add JNe to Fe, (mn|MN) contracted with DN 
C
C         Input: FockNA contains the classical-quantum nuclear
C         potential energy Hamiltonian
C
C         Output: FockNA contains the quantum nuclear Fock matrix
C         
C         Input: FockA contains the pure electronic Fock matrix
C
C         Output: FockA contains the electronic Fock matrix with
C         contributions from electron-quantum-nuclear interaction

        Call FoFNEO(IOut,IPrint,IUHFH,.False.,.True.,.True.,NBasis,
     $    NBas6D,NBsNEO,NBsNE6,Omega,.False.,0,0,0,
     $    IPFlag,AllowP,ISym2E,NAtoms,IAn,AtmChg,C,IAtTyp,JJ,NOpUse,
     $    NOp2,NEqAll,XX,RotAll,XX,XX,
     $    XX,AccDes,DenA,DenNA,FockA,FockB,FockNA,ESubV,V(IV),MDV1)

        Call LTOutN(IOut,'FN from l501:',0,1,NBsNEO,FockNA,0)
        Call LTOutN(IOut,'Fe from l501 after NEO contribution:',0,1,
     $      NBasis,FockA,0)
        Call FileIO(1,-IRwFA,NTTDR,FockA,0)
        Call FileIO(1,-IRwNEF,NEONTT,FockNA,0)

        Write(*,*) 'Classical nuclear potential = ',ENR
        Write(*,*) 'Quantum-classical nuclear potential = ',ESubV(3)
        Write(*,*) 'Quantum nuclear potential = ',ESubV(4)
        Write(*,*) 'Pure electronic energy =', ESubV(1)
        Write(*,*) 'Electron-quantum-nuclear potential = ',ESubV(2)
        Write(*,*) 'Quantum-nuclear-electron potential = ',ESubV(5)
        ECur = ENR + ESubV(1) + ESubV(2) + ESubV(3) + ESubV(4)

C       Put the total energy at the 6th index so that NEOSCF can access
C       this from the outside
        ESub = ECur
        DE = ECur - EPre

        If (IS2Con.eq.0) then
          Write(IOut, 1180) ECur, ISIter
        else
          Write(IOut, 1190) ECur, ISIter
          endif
        
        If (IS2Con.eq.0) then
      
C         Implementing DIIS for electronic subsystem:
          If (DoDIIS) then
            Call SCFDIS(IOut,IPrint,1,0,MaxCom,ISIter,1,IOpCl,NBasis,
     $        NBsUse,0,0,0,IRwPA,0,IRwFA,0,IRwPB,0,IRwFB,0,IRwS,IRwV,0,
     $        0,0,0,0,0,0,NFESav,IOBDIS,IOPADI,IOPBDI,IOFADI,IOFBDI,
     $        IOErrDI,0,0,0,0,V(IScrA),V(IScrC),V(IScrD),JJ,XX,JJ,
     $        V(IBMat),V(IScrCof),ErrMax,IndFok,Energy,EnCyc,ErrCyc,1,0,
     $        1,0,V(IV),MDV1)
            endif

          Call FileIO(2,-IRwFA,NTTDR,FockA,0)
C         Transform the electronic Fock matrix to the MO basis &&
C         Diagonalize the electronic Fock matrix
C          Call FileIO(2,-IRwV,NBasis*NBasis,CMOA,0)

          Call SimDag(IOut,IPrint,NRI,NBD,NOrbs,NAEX,FockA,EigE(1,1)
     $      ,CMOA,Scr1,V(IV),MDV1)
          If(NSpBlk.eq.2) Call SimDag(IOut,IPrint,NRI,NBD,NOrbs,NBE,
     $      FockB,EigE(1,2),CMOB,Scr1,V(IV),MDV1)

          Call FileIO(1,-IRwCA,LenCMO,CMOA,0)
          If(NSpBlk.eq.2) Call FileIO(1,-IRwCB,LenCMO,CMOB,0)
          
        else

C         Implementing DIIS for protonic subsystem:
          If (DoDIIS) then
            Call SCFDIS(IOut,IPrint,1,0,MaxCom,ISIter,1,IOpCl,NBsNEO,
     $        NBUNEO,0,0,0,IRwNED,0,IRwNEF,0,0,0,0,0,IRwNES,IRwNEV,0,0,
     $        0,0,0,0,0,NFPSav,IOBDIS,IOPADI,IOPBDI,IOFADI,IOFBDI,
     $        IOErrDI,0,0,0,0,V(IScrA),V(IScrC),V(IScrD),JJ,XX,JJ,
     $        V(IBMat),V(IScrCof),ErrMax,IndFok,Energy,EnCyc,ErrCyc,
     $        1,0,1,0,V(IV),MDV1)
            endif
          Call FileIO(2,-IRwNEF,NEONTT,FockNA,0)
          
C         Transform the protonic Fock matrices to the MO basis &&
C         Diagonalize the protonic Fock matrix
          Call SimDag(IOut,IPrint,NRI,NBsNEO,NBUNEO,NQNuc,FockNA,
     $      EigN(1,1),CMONA,Scr1,V(IV),MDV1)
          Call FileIO(1,-IRwNEC,NRI*NBsNEO*NBUNEO,CMONA,0)
          endif

  600   EPre = ECur
C    
C     If we reach here that means the electronic sub-system is not 
C     converged, we exit the outer loop with flag "Kill = 1"
      Kill = 1
       
  900 ISubCy = ISIter 
      Call FreeRw(IRwC1)
      Call FreeRw(IRwC2)
      Call FreeRw(IRwC3)
      Return
      End

*Deck SimDag
      Subroutine SimDag(IOut,IPrint,NRI,NBD,NOrbs,NOcc,Fock,Eig,CMO,
     $  Scr1,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Transform a Fock matrix the MO basis, diagonalize, and update CMO.
C
      Parameter (Zero=0.0d0,One=1.0d0,MinPrt=2)
      Dimension Fock(*), Eig(*), CMO(*), Scr1(*), V(*), NBF(1),
     $  NOccSy(1)
C
      Call XVFVMu(NRI,.False.,.False.,.False.,0,NBD,NOrbs,NBD,Fock,
     $  CMO,Scr1,Fock,V,MDV)
C     IOpClX tells LTOutX whether matrix is real or complex.
      IOpClX = 2*(NRI-1)
      If(IPrint.ge.MinPrt) Call LTOutX(IOut,0,IOpClX,1,0,'MO Fock',
     $  NOrbs,NOrbs,Fock,Fock,0)
      NBF(1) = NOrbs
      NOccSy(1) = NOcc
      Call DiagDS(IOut,IPrint,0,.False.,NRI,1,NBF,NOrbs,NOccSy,Zero,
     $  OVMax,Fock,Scr1,Eig,MDV,V,XX)
      NBDSqR = NRI*NBD*NBD
      Call TstCor(NBDSqR,MDV,'SimDag')
C      Call XGEMM(NRI,'N','N',NBD,NOrbs,NBD,One,Scr1,NBD,CMO,NBD,Zero,V,
C     $  NBD)
      Call XGEMM(NRI,'N','N',NBD,NOrbs,NBD,One,CMO,NBD,Scr1,NBD,Zero,V,
     $  NBD)
      Call AMove(NBDSqR,V,CMO)
      Return
      End


*Deck SMTSCF
      Subroutine SMTSCF(IOut,IPrint,IOpCl,IExtp,ICharg,Multip,IAn,
     $  IAtTyp,AtmChg,C,NAE,NBE,NAEX,NRI,NSpBlk,NSpBlX,NSpBNX,NDimBl,
     $  NAtoms,NBasis,NBD,NBas6D,NBsUse,NOrbs,NTTD,NTTDR,NQNuc,NBsNEO,
     $  NBsNE6,NBUNEO,NEONTT,NEONT6,IHMeth,IPFlag,AllowP,LSEAll,RSEAll,
     $  Omega,ICtDFT,IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,IGWInf,IAtBtD,
     $  IRdBtD,RRdBtD,RGWBtD,NEqAll,IRwFA,IRwFB,IRwS,IRwH,IRwCA,IRwCB,
     $  IRwPA,IRwPB,IRwV,IRwNEX,IRwNED,IRwNEC,IRwNEF,IRwNES,IRwNEV,
     $  CMOA,CMOB,DenA,DenB,FockA,FockB,CMONA,CMONB,DenNA,DenNB,FockNA,
     $  FockNB,EigE,EigN,ESub,IECnt,IPCnt,EDone,PDone,Acurcy,AccDes,
     $  AccMax,DE,ETol,ENR,IRwCE1,IRwCE2,IRwCE3,IRwCP1,IRwCP2,IRwCP3,
     $  MaxCom,IOBDIS,IOPAE,IOPBE,IOFAE,IOFBE,IOErrE,IOPNA,IOFNA,IOErrN,
     $  NFkSav,BMat,ScrCof,EnCyc,ErrCyc,IndFok,Scr1,NIter,DoDIIS,V,MDV)
     
        Implicit Real*8(A-H,O-Z)
C
C     This subroutine simultaneously optimize the
C     electronic or the protonic sub-system
C

      Parameter (Zero=0.0d0,One=1.0d0,SmalEn=1.d-9,MinPrt=2,
     $  F100=100.0d0)
      Logical AllowP(*), DoneQ, Rises, Damp, Done, EDone, PDone, DoDIIS
      Character*1 JCG
      Dimension IAn(*), IAtTyp(*), AtmChg(*), C(*), RSEAll(*), Omega(*),
     $  ScaDFX(*), IGWInf(*), IAtBtD(*), IRdBtD(*), RRdBtD(*),
     $  RGWBtD(*), NEqAll(*),EigE(NOrbs,*),EigN(NBsNEO,*),
     $  IndFok(MaxCom),
     $  EnCyc(MaxCom),ErrCyc(MaxCom),CMOA(NRI,NBD,NBD),
     $  CMOB(NRI,NBD,NBD), DenA(NRI,NTTD), DenB(NRI,NTTD),
     $  FockA(NRI,NTTD), FockB(NRI,NTTD), CMONA(NRI,NBsNEO,*),
     $  CMONB(NRI,NBsNEO,*), DenNA(NRI,NEONT6), DenNB(NRI,NEONT6),
     $  FockNA(NRI,NEONT6), FockNB(NRI,NEONT6),BMat(MaxCom+1,*),
     $  Scr(MaxCom+1,7),
     $  Scr1(NRI,Max(NBD,NBsNEO),Max(NBD,NBsNEO)),
     $  V(MDV), JJ(1), XX(1), ESubV(10),
     $  ECDFit(1), EOne(1),
     $  ETwo(1), EX(1), EC(1), E6(1), ERFP(1), EPMM(1), RNE(1)
C     JJ and XX are dummies for arguments which are unused in the
C     calls made here.
      Save JJ, XX
      Data JJ/0/, XX/0.0d0/
 1120 Format(' E=',G22.15,' Delta-E=',F22.12,' Rises=',L1,' Damp=',L1,
     $  ' RMSDP=',1PD8.2,' MaxDP=',1PD8.2)
 1130 Format(' E=',G22.15)
 1170 Format(' E=','   (Non-Variational)  ')
 1180 Format(' In Electronic SubSCF Loop, Total Energy = ', F22.15, 
     $ ' After Iteration # ', I4)     
 1190 Format(' In Protonic SubSCF Loop, Total Energy = ', F22.15, 
     $ ' After Iteration # ', I4)     
C
      Write(*,*)'Enter Simultaneous Opt'
C     Calculate offsets in V so that we have scratch space for SCFDIS
C     subroutine
      IScrAE = 1
      LADimE = Max(NBasis,MaxCom+1)
      LADimP = Max(NBsNEO,MaxCom+1)
      IScrAP = IScrAE + LADimE*(LADimE+1)
      IScrCE = IScrAP + LADimP*(LADimP+1)
      IScrCP = IScrCE + LADimE*LADimE
      IScrDE = IScrCP + LADimP*LADimP
      LDDimE = Max(NBasis,(MaxCom+1)**2)
      LDDimP = Max(NBsNEO,(MaxCom+1)**2)
      IScrDP = IScrDE + LDDimE
      IV = IScrDP + LDDimP
      Call TstCor(IV,MDV,'SMTSCF')
      MDV1 = MDV - IV + 1

      Call RdNERw(IOut,IOpCl,0,IRwNEX,NQNuc,NBsNEO,NBsNE6,NBUNEO,
     $  IOfNES,IOfNET,IOfNEH,IOfNER,IOpClE,IOClNP,IOClND)
      LenCMO = NRI*NBD*NOrbs

C     Whether elec/prot sub-system has converged
      EDone = .False.
      PDone = .False.

C     Flags for Convergence Check
      IFlagO = 0
      IFlag = 0
      Damp = .False.
      DampFc = One
 
C     Form the electronic density matrix and save it in IRwPA
      Call FormP2(.True.,NRI,NBD,NAEX,0,Zero,CMOA,DenA)
      If(NSpBlk.eq.2) Call FormP2(.True.,NRI,NBD,NBE,0,Zero,CMOB,DenB) 
      Call FileIO(1,-IRwPA,NTTDR,DenA,0)
      If(NSpBlk.eq.2) Call FileIO(1,-IRwPB,NTTDR,DenB,0)


C     Form the electronic density matrix and save it in IRwNED 
      Call FormP2(.True.,NRI,NBsNEO,NQNuc,0,Zero,CMONA,DenNA)
      Call FileIO(1,-IRwNED,NEONTT,DenNA,0)
 
      Call LTOutN(IOut,'New Electronic Density:',0,1,NBasis,DenA,0)
      Call LTOutN(IOut,'New Protonic Density:',0,1,NBsNEO,DenNA,0)
          

C     Check if the electronic density has converged
      Call ConOpn(IOut,IOpCl,NBasis,NBsUse,RMSDPE,IFlag,IExtp,DenA,
     $  FockA,Scr1,Acurcy,AccMax,Zero,DE,DelMax,IECnt,EDone,
     $  IRwPA,IRwPB,IRwCE1,IRwCE2,IRwCE3,ETol,1,1,Damp,DampFc)


C     Check the convergence of the protonic density matrix
      Call ConOpn(IOut,IOpCl,NBsNEO,NBUNEO,RMSDPP,IFlag,IExtp,DenNA,
     $  FockNA,Scr1,Acurcy,AccMax,Zero,DE,DelMax,IPCnt,PDone,
     $  IRwNED,IRwNED,IRwCP1,IRwCP2,IRwCP3,ETol,1,1,Damp,DampFc)

      Write(*,*) "RMSDPE",RMSDPE
      Write(*,*) "RMSDPP",RMSDPP

C     In ConOpn, DenA/DenNA will be holding the P_{n-1}, so we need to
C     read them back in from where they are stored.
      Call FileIO(2,-IRwPA,NTTDR,DenA,0)
      Call FileIO(2,-IRwNED,NEONTT,DenNA,0) 

      NOpUse = 1
      NBTI = 0
      AccXC = AccDes
      NSaved = 0
      DoneQ = .False.

      Call MkFok1(IOut,4,IOpCl,.False.,.False.,.False.,.False.
     $  ,0,
     $  1,1,0,1,NAE,NBE,ICharg,Multip,NBasis,NBas6D,NBsUse,DenA,
     $  DenB,
     $  FockA,FockB,1,NOpUse,1,1,.True.,0,IPFlag,AllowP,.False.
     $  ,0,0,0,
     $  0,0,IHMeth,LSEAll,RSEAll,Omega,NBTI,AccDes,AccXC,ICtDFT,
     $  IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,0,IAn,IAtTyp,AtmChg,C,
     $  NAtoms,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,.True.,0,1,0,0,.False.,
     $  NSaved,1,XX,0,0,.False.,1,JJ,JJ,JJ,JJ,XX,JJ,JJ,XX,XX,XX,JJ,
     $  0,
     $  XX,JCG,ESubV,EElec,ECDFit,EOne,ETwo,Ex,Ec,E6,ErfP,EPMM,DoneQ,XX,
     $  One,XX,JJ,JJ,JJ,NAtoms,NEqAll,XX,XX,RNE,JJ,JJ,IGWInf,IAtBtD,
     $  IRdBtD,RRdBtD,RGWBtD,0,XX,1,JJ,XX,XX,NBsNEO,DenNA,DenNB,
     $  FockNA,FockNB,XX,XX,0,XX,V(IV),MDV1)
      Call LTOutN(IOut,'Fee from MkFok1:',0,1,NBasis,FockA,0)
      Call FileIO(2,-IRwNEX,NEONTT,FockNA,IOfNEH)
 
C     TODO: Symmetry is turned off 
      ISym2E = 0
      NOpUse = 0
      NOp2 = 0
C     TODO: RHF only
      IUHFH = 0

C     Add JNe to Fe, (mn|MN) contracted with DN 
C
C     Input: FockNA contains the classical-quantum nuclear
C     potential energy Hamiltonian
C
C     Output: FockNA contains the quantum nuclear Fock matrix
C     
C     Input: FockA contains the pure electronic Fock matrix
C
C     Output: FockA contains the electronic Fock matrix with
C     contributions from electron-quantum-nuclear interaction

      Call FoFNEO(IOut,IPrint,IUHFH,.False.,.True.,.True.,NBasis,
     $  NBas6D,NBsNEO,NBsNE6,Omega,.False.,0,0,0,
     $  IPFlag,AllowP,ISym2E,NAtoms,IAn,AtmChg,C,IAtTyp,JJ,NOpUse,
     $  NOp2,NEqAll,XX,RotAll,XX,XX,
     $  XX,AccDes,DenA,DenNA,FockA,FockB,FockNA,ESubV,V(IV),MDV1)

      Call LTOutN(IOut,'FN from l501:',0,1,NBsNEO,FockNA,0)
      Call LTOutN(IOut,'Fe from l501 after NEO contribution:',0,1,
     $    NBasis,FockA,0)


      Call FileIO(1,-IRwFA,NTTDR,FockA,0)
      Call FileIO(1,-IRwNEF,NEONTT,FockNA,0)

      Write(*,*) 'Classical nuclear potential = ',ENR
      Write(*,*) 'Quantum-classical nuclear potential = ',ESubV(3)
      Write(*,*) 'Quantum nuclear potential = ',ESubV(4)
      Write(*,*) 'Pure electronic energy =', ESubV(1)
      Write(*,*) 'Electron-quantum-nuclear potential = ',ESubV(2)
      Write(*,*) 'Quantum-nuclear-electron potential = ',ESubV(5)
      ESub = ENR + ESubV(1) + ESubV(2) + ESubV(3) + ESubV(4)


      If (DoDIIS) then
        Call SCFDIS(IOut,IPrint,1,0,MaxCom,NIter,1,IOpCl,NBasis,NBsUse,
     $    NBsNEO,NBuNEO,0,IRwPA,0,IRwFA,0,IRwPB,0,IRwFB,0,IRwS,IRwV,
     $    IRwNED,IRwNEF,IRwNES,IRwNEV,0,0,0,NFkSav,IOBDIS,IOPAE,IOPBE,
     $    IOFAE,IOFBE,IOErrE,IOPNA,IOFNA,IOErrN,0,V(IScrAE),V(IScrCE),
     $    V(IScrDE),V(IScrAP),V(IScrCP),V(IScrDP),BMat,ScrCof,
     $    ErrMax,IndFok,Energy,EnCyc,ErrCyc,1,0,1,0,V(IV),MDV1)
        endif

      Call FileIO(2,-IRwFA,NTTDR,FockA,0)
      Call FileIO(2,-IRwNEF,NEONTT,FockNA,0)

C     We know diagonalize both electronic and protonic fock matrices        
      Call SimDag(IOut,IPrint,NRI,NBD,NOrbs,NAEX,FockA,EigE(1,1)
     $  ,CMOA,Scr1,V(IV),MDV1)
      If(NSpBlk.eq.2) Call SimDag(IOut,IPrint,NRI,NBD,NOrbs,NBE,
     $  FockB,EigE(1,2),CMOB,Scr1,V(IV),MDV1)
      Call FileIO(1,-IRwCA,LenCMO,CMOA,0)
      If(NSpBlk.eq.2) Call FileIO(1,-IRwCB,LenCMO,CMOB,0)
        
      Call SimDag(IOut,IPrint,NRI,NBsNEO,NBUNEO,NQNuc,FockNA,
     $  EigN(1,1),    CMONA,Scr1,V(IV),MDV1)
      Call FileIO(1,-IRwNEC,NRI*NBsNEO*NBUNEO,CMONA,0)

       
      Return
      End

*Deck SimCyc
      Subroutine SimCyc(IOut,IPrint,IOpCl,IExtp,NRI,NSpBlk,NSpBlX,
     $  NDimBl,NAtoms,NBas6D,NBasis,NBsUse,NBD,NTTD,NOrbs,
     $  ICharg,Multip,NAE,NBE,IAn,IAtTyp,AtmChg,C,IRwH,
     $  IRwCA,IRwCB,IRwPA,IRwPB,IHMeth,IPFlag,AllowP,LSEAll,
     $  RSEAll,Omega,ICtDFT,IExCor,N6Func,ScaHFX,IGWInf,IAtBtD,
     $  IRdBtD,RRdBtD,RGWBtD,ScaDFX,ICorTp,NEqAll,Eig,
     $  CMOA,CMOB,DenA,DenB,FockA,FockB,ENR,
     $  Energy,JCycle,KilSCF,MaxCyc,ETol,Acurcy,AccDes,
     $  AccMax,Scr1,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Simple Classical Hartree-Fock SCF loop.
C
      Parameter (Zero=0.0d0,One=1.0d0,SmalEn=1.d-9,MinPrt=2,
     $  F100=100.0d0)
      Logical AllowP(*), DoneQ, Rises, Damp, Done
      Character*1 JCG
      Dimension IAn(*), IAtTyp(*), AtmChg(*), C(*), RSEAll(*), Omega(*),
     $  ScaDFX(*), IGWInf(*), IAtBtD(*), IRdBtD(*), RRdBtD(*),
     $  RGWBtD(*), NEqAll(*), Eig(NOrbs,*), CMOA(NRI,NBD,NBD),
     $  CMOB(NRI,NBD,NBD), DenA(NRI,NTTD), DenB(NRI,NTTD),
     $  FockA(NRI,NTTD), FockB(NRI,NTTD), Scr1(NRI,NBD,NBD),
     $  V(MDV), JJ(1), XX(1), ESum(10), ECDFit(1), EOne(1),
     $  ETwo(1), EX(1), EC(1), E6(1), ERFP(1), EPMM(1), RNE(1)
C     JJ and XX are dummies for arguments which are unused in the
C     calls made here.
      Save JJ, XX
      Data JJ/0/, XX/0.0d0/
 1120 Format(' E=',G22.15,' Delta-E=',F22.12,' Rises=',L1,' Damp=',L1,
     $  ' RMSDP=',1PD8.2,' MaxDP=',1PD8.2)
 1130 Format(' E=',G22.15)
 1170 Format(' E=','   (Non-Variational)  ')
C
      Write(*,*)'CONVENTIONAL SCF'
C     NAEX is the number of occupieds for the first or only spin-block.
C     It is NAE for spin-aligned or the total number of electrons for GHF.
      NAEX = NAE + (NDimBl-1)*NBE
      NTTD = (NBD*(NBD+1))/2
      NTTDR = NRI*NTTD
Cxslis
C
C
Cxslie
      IFlag = 0
      ICount = 0
      Damp = .False.
      DampFc = One
      DelEn = Zero
      EnOld = Zero
      KilSCF = 0
      IFlagO = 0
      LenCMO = NRI*NBD*NOrbs
C     Scratch files used in extrapolation.
      LScr = NSpBlk*NTTDR
      IRwC1 = IGetRZ(LScr,0,V,MDV)
      IRwC2 = IGetRZ(LScr,0,V,MDV)
      IRwC3 = IGetRZ(LScr,0,V,MDV)
C
C     Top of iteration loop.
C
      Do 500 Iter = 1, MaxCyc
C
C       Form electronic density matrices.
C
        Call FormP2(.True.,NRI,NBD,NAEX,0,Zero,CMOA,DenA)
        If(NSpBlk.eq.2) Call FormP2(.True.,NRI,NBD,NBE,0,Zero,CMOB,DenB)
        Call FileIO(1,-IRwPA,NTTDR,DenA,0)
        If(NSpBlk.eq.2) Call FileIO(1,-IRwPB,NTTDR,DenB,0)
C

C       Check change in density and possibly extrapolate densities.
C
        Call ConOpn(IOut,IOpCl,NBasis,NBsUse,RMSDPE,IFlag,IExtp,DenA,
     $    FockA,Scr1,Acurcy,AccMax,Zero,DelEn,DelMax,ICount,Done,
     $    IRwPA,IRwPB,IRwC1,IRwC2,IRwC3,ETol,1,1,Damp,DampFc)
        Call FileIO(2,-IRwPA,NTTDR,DenA,0)
        If(NSpBlk.eq.2) Call FileIO(2,-IRwPB,NTTDR,DenB,0)
        If(Done) goto 900
C
C       Form Fock matrices.  Here we always compute the 2e integrals
C       and don't worry about symmetry, incremental Fock formation,
C       FMM, etc.
C
        If(IPrint.ge.0) Call LTOutX(IOut,0,IOpCl,1,0,
     $    'Electron Densities',NBasis,NBasis,DenA,DenB,0)
        NOpUse = 1
        NBTI = 0
        AccXC = AccDes
        NSaved = 0
        DoneQ = .False.
        Call MkFok1(IOut,IPrint,IOpCl,.False.,.False.,.False.,.False.,0,
     $    1,1,0,1,NAE,NBE,ICharg,Multip,NBasis,NBas6D,NBsUse,DenA,DenB,
     $    FockA,FockB,1,NOpUse,1,1,.True.,0,IPFlag,AllowP,.False.,0,0,0,
     $    0,0,IHMeth,LSEAll,RSEAll,Omega,NBTI,AccDes,AccXC,ICtDFT,
     $    IExCor,N6Func,ScaHFX,ScaDFX,ICorTp,0,IAn,IAtTyp,AtmChg,C,
     $    NAtoms,IRwH,IRwCA,IRwCB,IRwPA,IRwPB,.True.,0,1,0,0,.False.,
     $    NSaved,1,XX,0,0,.False.,1,JJ,JJ,JJ,JJ,XX,JJ,JJ,XX,XX,XX,JJ,0,
     $    XX,JCG,ESum,EElec,ECDFit,EOne,ETwo,Ex,Ec,E6,ErfP,EPMM,DoneQ,
     $    XX,
     $    One,XX,JJ,JJ,JJ,NAtoms,NEqAll,XX,XX,RNE,JJ,JJ,IGWInf,IAtBtD,
     $    IRdBtD,RRdBtD,RGWBtD,0,XX,1,JJ,XX,XX,NBsNEO,DenNA,DenNB,
     $    FockNA,FockNB,XX,XX,0,XX,V,MDV)
Cxslie
C   
        Energy = EElec + ENR
        Write(*,*)'ITERATION: ', Iter
        Write(*,*)'EElec: ', EElec
        Write(*,*)'ENR:', ENR
        Write(*,*)'ENERGY: ', Energy

C       Diagonalize electronic matrix        
        Call SimDag(IOut,IPrint,NRI,NBD,NOrbs,NAEX,FockA,Eig,CMOA,Scr1,
     $    V,MDV)
        If(NSpBlk.eq.2) Call SimDag(IOut,IPrint,NRI,NBD,NOrbs,NBE,FockB,
     $    Eig(1,2),CMOB,Scr1,V,MDV)
        Call FileIO(1,-IRwCA,LenCMO,CMOA,0)
        If(NSpBlk.eq.2) Call FileIO(1,-IRwCB,LenCMO,CMOB,0)

CAL     
        EnOld = Energy
  500   IFlagO = IFlag
C     Failed to converge.
      Iter = MaxCyc + 1
      Kill = 1
  900 JCycle = Iter - 1
      Call FreeRw(IRwC1)
      Call FreeRw(IRwC2)
      Call FreeRw(IRwC3)
      Return
      End
