*Deck MkFok1
      Subroutine MkFok1(IOut,IPrint,IOpCl,ROHF,RComp,JustJ,IniPCM,
     $  JSym2X,NMatI,NMatSI,NMatTI,NMtPBC,NAE,NBE,ICharg,Multip,NBasis,
     $  NBas6D,NBsUse,PA,PB,FA,FB,NOpAll,NOpUse,NOp1,NOp2,DirSCF,ICntrl,
     $  IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IRaf,KAlg,IHMeth,LSEAll,
     $  RSEAll,Omega,NBTI,AccDes,AccXC,ICtDFT,IExCor,N6Func,ScaHFX,
     $  ScaDFX,ICorTp,ISCRF,IAn,IAtTyp,AtmChg,C,NAtoms,IRwHI,IRwCA,
     $  IRwCB,IRwPA,IRwPB,DoFull,InStor,IDenFit,IODFit,IOFDen,IncFoc,
     $  NSaved,NDimCP,CP,IOPSav,IOFSav,FrcNuc,NSAOBk,NS,NS2,NS4,ISABF,
     $  RSABF,MulTab,MapSAO,R1,R2,R3,IPBC,LIPBC,CelVec,JGC,Energy,EElec,
     $  ECDFit,EOne,ETwo,Ex,Ec,E6,ErfP,EPMM,DoneQ,gFac,fuFac,EField,
     $  NEqShl,NEqSh2,NEqBas,NDimEq,NEqAtm,RotOp,Buf2E,RNE,IJ,IJMap,
     $  IGWInf,IAtBtD,IRdBtD,RRdBtD,RGWBtD,LPMAll,PMMAll,NDAtP,IAtPar,
     $  RAtPar,Conect,NBsNEO,PNA,PNB,FNA,FNB,PMMA,XDMA,LenXCU,XCParU,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Form Fock matrices, given densities, symmetry info, and symmetry
C     settings.
C
#include "molpars.inc"
      Logical DirSCF, AllowP(50), FMM, DFT, HavBRI, DFTSym, HarTHm,
     $  IfDFT, DoFull, ROHF, CBfn, Spinor, IncFoc, FrcNuc, HaveB,
     $  JustJ, GHF, RComp, FAneFB, CopyFB, PCM, Onsagr, SolEPt,
     $  Init, DoneQ, IniPCM, ClShel, PMOHam, DoDFT, IfPE, IfMMP
      Character*1 JGC
      Integer FMFlag, FMFlg1
      Real*8 MDCutO
      Dimension PA(*), PB(*), FA(*), FB(*), IAn(*), IAtTyp(*),
     $  AtmChg(*), C(3,*), R1(*), R2(*), R3(*), IPBC(*), CelVec(3,3),
     $  Energy(*), ScaDFX(*), Omega(*), RSEAll(*), NS(*), NS2(*),
     $  NS4(*), ISABF(*), RSABF(*), MulTab(*), MapSAO(*), NEqAtm(*),
     $  NEqBas(*), NEqShl(*), Ex(*), Ec(*), NEqSh2(*), RotOp(*),
     $  Buf2E(*), RNE(*), IJ(*), IJMap(*), EOne(*), ETwo(*), ECDFit(*),
     $  CP(NDimCP,*), IGWInf(*), IAtBtD(*), IRdBtD(*), RRdBtD(*),
     $  RGWBtD(*), ErfP(*), gfac(*), EField(0:*), PMMAll(*), XX(1),
     $  IAtPar(NDAtP,*), RAtPar(NDAtP,*), Conect(*), XDMA(*), V(*),
     $  EPMM(*), PNA(*), PNB(*), FNA(*), FNB(*), PMMA(*), XCParU(*),
     $  E6(*), JJ(1)
      Parameter (IRwHDf=515)
      Save Zero, Half, One, Two, JJ, XX
      Data Zero/0.0d0/, Half/0.5d0/, One/1.0d0/, Two/2.0d0/, JJ/0/,
     $  XX/0.0d0/
 1000 Format(' IMat=',I3,' EOne=',D19.12,' ETwo=',D19.12,' Energy=',
     $  D19.12)
 1010 Format(' IMat=',I3,' 1/2 <P',A,'(P)>=',F25.12)
 1030 Format(' IMat=',I3,' Ex=',F25.12,' Ec=',F25.12,' E6=',F25.12,/,
     $  ' 1/2 <PJ(P)> + E(ex) + E(corr) + E(6int)=',F25.12)
 1040 Format(' Calling FoFDFT with Acc=',D9.2,' ICntrl=',I8,' IRadAn=',
     $  I8,'.')
 1050 Format(' Restarting incremental Fock formation.')
 1060 Format(' Stop after computing J.')
C
      Call CkSCFT(IOut,'MkFok1',IOpCl,.True.,.True.,.True.,.False.,
     $  NSpBlk,NRI,NDimBl,CBfn,NSpBlX,NSpBNX,Spinor)
      If(IGWInf(1).eq.0) then
        IRadAn = 0
      else
        IRadAn = IGWInf(3)
        endIf
      Call AClear(NMatI,ErfP)
      Call AClear(NMatI,EPMM)
      GHF = Spinor.or.NDimBl.eq.2
      HaveB = NSpBlk.eq.2.or.ROHF
      HavBRI = HaveB.or.(NDimBl.eq.2.and..not.CBfn)
      FAneFB =.not.RComp
      If(HaveB) then
        NSpBk1 = 2
        IOpClY = 1
      else
        NSpBk1 = NSpBlk
        IOpClY = NSpBlX - 1
        endIf
      CopyFB = NSpBk1.eq.1.and.NDimBl.eq.1.and..not.CBfn.and.FAneFB
      IPrt1 = IPrint
      If(IPSave(0).eq.0) IPrt1 = Max(IPrt1,1)
      DFT = IfDFT(IExCor,N6Func).or.HarTHm(IHMeth)
      NBD = NBasis*NDimBl
      NTTD = (NBD*(NBD+1))/2
      NTTDR = NRI*NTTD
      NTTX = NTTDR*NMtPBC
      NTTM = NTTX*NMatI
      IV = 1
      MDV1 = MDV - IV + 1
      IRwSvP = IGetRW(NTTM*NSpBk1,0)
      Call FileIO(1,-IRwSvP,NTTM,PA,0)
      If(HaveB) Call FileIO(1,IRwSvP,NTTM,PB,0)
      Call PcckI(0,Junk,NEqAtm,NDimEq,NOpAll,NAtoms,NOpAll)
      Call IndSE(IHMeth,1,RSEAll,IHMPar,MaxASE,NSEScI,NSEScR,NSEAtI,
     $  NSEAtR,NCorPA,MaxCPA,NSEElI,NSEElR,MxElBE,jRSESc,jISEP,jRSEP,
     $  jSECP,jIElP,jRElP,jBELP,jBEBO,Junk)
C
C     SclInc is used to scale the nuclear contribution to the Fock matrix
C     when it is done in FoFCou so that incremental Fock matrix formation
C     still works.
C
      SclInc = One
      If(IncFoc) then
        If(NMatI.gt.1) Call GauErr('No incremental Fock with NMatI>1.')
        Call TstCor(NTTX,MDV1,'MkFok1-IncCof')
        If(NSaved.ge.(NDimCP-1)) then
          Write(IOut,1050)
          NSaved = 0
          endIf
        Call IncCof(IOut,IPrt1,IOpCl,NMtPBC,NBasis,NBsUse,NTTDR,NSaved,
     $    IOPSav,V(IV),SclInc,NDimCP,CP,PA,PB)
        Call FileIO(1,-IOPSav,NTTX,PA,NSaved*NSpBlk*NTTX)
        If(NSpBlk.eq.2) Call FileIO(1,IOPSav,NTTX,PB,0)
      else
        NSaved = 0
        endIf
      CoulEn = Zero
      If(PMOHam(IHMeth)) Call XToXB(IOut,IPrint,IOpCl,IHMeth,.True.,
     $  .True.,NBasis,NMatI,NMatSI+NMatTI,NTTDR,PA,PB,V(IV),MDV1)
      Call FoF2E(IOut,IPrint,IOpCl,ROHF,DFT,.True.,FAneFB,NMatI,NMatSI,
     $  NMatTI,NMtPBC,NBasis,NBas6D,NAtoms,IAn,IAtTyp,AtmChg,C,
     $  IAtPar(1,MIAtFr),DirSCF,InStor,SclInc,FrcNuc,IDenFit,DoFull,
     $  ICntrl,IHMeth,IRaf,ScaHFX,AccDes,FMM,FMFlag,FMFlg1,NFxFlg,
     $  IPFlag,AllowP,LSEAll,RSEAll,KAlg,Omega,NBTI,IODFit,IOFDen,
     $  JSym2X,NOp1,NOp2,NOpUse,NEqAtm,NEqBas,NEqShl,NEqSh2,RotOp,IJ,
     $  IJMap,NSAOBk,NS,NS2,NS4,ISABF,RSABF,MulTab,MapSAO,JGC,CoulEn,
     $  ECDFit,RNE,PA,PB,FA,FB,NBsNEO,PNA,PNB,FNA,FNB,R1,R2,R3,IPBC,
     $  LIPBC,CelVec,Buf2E,V(IV),MDV1)
      If(PMOHam(IHMeth)) Call XToXB(IOut,IPrint,IOpCl,IHMeth,.False.,
     $  .False.,NBasis,NMatI,NMatSI+NMatTI,NTTDR,FA,FB,V(IV),MDV1)
      Call UnPckI(0,Junk,NEqAtm,NDimEq,NOpAll,NAtoms,NOpAll)
      If(IncFoc) then
        NS1 = NSaved + 1
        If(Abs(CP(NS1,2)).ge.MDCutO(0)) then
          CP(NS1,3) = (CoulEn+ECDFit(1)) / CP(NS1,2)
        else
          CP(NS1,3) = Zero
          endIf
        Call FileIO(1,-IOFSav,NTTX,FA,NSaved*NSpBlk*NTTX)
        If(NSpBk1.eq.2) Call FileIO(1,IOFSav,NTTX,FB,0)
        If(NSaved.gt.0) then
          Call TstCor(NTTX,MDV1,'MkFok1-IncSum')
          Call FileIO(2,-IOFSav,0,CP,0)
          CoulEn = CP(NS1,3) + SProd(NSaved,CP,CP(1,3))
          Do 10 IMat = 1, NSaved
            Call FileIO(2,IOFSav,NTTX,V(IV),0)
            Call ACasB(NTTX,FA,V(IV),FA,CP(IMat,1))
            If(NSpBk1.eq.2) then
              Call FileIO(2,IOFSav,NTTX,V(IV),0)
              Call ACasB(NTTX,FB,V(IV),FB,CP(IMat,1))
              endIf
   10       Continue
          If(CopyFB) Call AMove(NTTX,FA,FB)
          endIf
        NSaved = NS1
        endIf
      Call FileIO(2,-IRwSvP,NTTM,PA,0)
      If(HaveB) Call FileIO(2,IRwSvP,NTTM,PB,0)
      Call TStamp(1,'After 2e')
C
C     Add in BEBO terms if present.
C
      If(MxElBE.gt.0) then
        Call FoFBE(IOut,IPrint,.True.,IOpCl,IHMeth,MxElBE,RSEAll(jBELP),
     $    RSEAll(jBEBO),NMatI,NMatSI,NAtoms,NBasis,IAn,NTTDR,PA,PB,FA,
     $    FB,V(IV),MDV1)
        Call TStamp(1,'After BEBO')
        endIf
C
C     For the DoFull case, add in the core Hamiltonian and compute the energy.
C
      IRwH = IRwHI
      If(IRwH.eq.0) IRwH = IRwHDf
      If(DoFull) then
        Call TstCor(NTTX*NSpBlk,MDV1,'MkFok1-Core')
        Call FileIO(2,-IRwH,NTTX*NSpBlk,V(IV),0)
        IOpClX = IOpCl - Mod(IOpCl,2)
        Do 20 I = 1, NMatI
          Ind = 1 + (I-1)*NTTX
          ETwo(I) = SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,PA(Ind),
     $      XX,FA(Ind),XX)
          If(I.le.NMatSI) then
            EOne(I) = SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,
     $        PA(Ind),XX,V(IV),XX)
            Call AAdd(NTTX,FA(Ind),V(IV),FA(Ind))
          else
            EOne(I) = Zero
            endIf
          If(HaveB) then
            IndH = IV + NTTX*(NSpBlk-1)
            ETwo(I) = Half*(ETwo(I)+SCFTrX(.False.,.False.,IOpClX,
     $        NMtPBC,NBasis,PB(Ind),XX,FB(Ind),XX))
            If(I.le.NMatSI) then
              EOne(I) = EOne(I) + SCFTrX(.False.,.False.,IOpClX,NMtPBC,
     $          NBasis,PB(Ind),XX,V(IndH),XX)
              Call AAdd(NTTX,FB(Ind),V(IndH),FB(Ind))
              endIf
          else if(GHF) then
            ETwo(I) = Half*ETwo(I)
          else
            EOne(I) = Two*EOne(I)
            If(CopyFB) Call AMove(NTTX,FA,FB)
            endIf
          ETwo(I) = ETwo(I) + CoulEn
   20     Continue
        endIf
C
C     Do exchange-correlation integration for DFT and PCM terms.
C
      If(JustJ) then
        Write(IOut,1010) (I,JGC,ETwo(I),I=1,NMatI)
        Write(IOut,1060)
        Call Chain(998,0,V,MDV)
        endIf
      DoneQ =.not.PCM(ISCRF)
      DoDFT = DFT.and..not.HarTHm(IHMeth)
      ITyPMM = IPMMTp(LPMAll,PMMAll)
      If(DoDFT.or.PCM(ISCRF).or.ITyPMM.ne.0) then
        If(JSym2X.eq.1) then
          NOpDFT = NOpUse
        else
          NOpDFT = 1
          endIf
        DFTSym = NOpDFT.gt.1
        If(DoFull) then
          ICtDF1 = (ICtDFT/100)*100 + 11
          If(DFT.and.IPrint.gt.0) Write(IOut,1040) AccXC, ICtDF1, IRadAn
          IRwSvF = IGetRW(NTTM*NSpBk1,0)
          Call FileIO(1,-IRwSvF,NTTM,FA,0)
          If(HaveB) Call FileIO(1,IRwSvF,NTTM,FB,0)
          Do 30 I = 1, NMatSI
            If(IPrint.gt.0) Write(IOut,1010) I, JGC, ETwo(I)
            Init = IniPCM.and.I.eq.1
            Ind = (I-1)*NTTX
            Call FileIO(2,-IRwSvP,NTTX,PA,Ind)
            Call FileIO(2,-IRwSvF,NTTX,FA,Ind)
            If(HaveB) then
              Call FileIO(2,-IRwSvP,NTTX,PB,Ind+NTTM)
              Call FileIO(2,-IRwSvF,NTTX,FB,Ind+NTTM)
              endIf
            Call RIComp(.True.,HaveB,.False.,.False.,.True.,IOpCl,
     $        .False.,NBasis,NMtPBC,1,1,1,0,NMatY,NMatSY,NMatTY,PA,PB,
     $        V(IV),MDV1)
            Call RIComp(.True.,HaveB,.True.,.False.,.True.,IOpCl,
     $        .False.,NBasis,NMtPBC,1,1,1,0,NMatY,NMatSY,NMatTY,FA,FB,
     $        V(IV),MDV1)
            If(DoDFT) then
              Call FoFDFT(IOut,IPrint,IPFlag,IOpCl,HavBRI,ICtDF1,IExCor,
     $          N6Func,ICorTp,ScaDFX,DFTSym,NOpDFT,NOp1,1,NAtoms,NAE,
     $          NBE,NBasis,NBas6D,ICharg,Multip,IRwCA,0,IRwCB,AccXC,IAn,
     $          IAtTyp,AtmChg,C,Ex(I),Ec(I),E6(I),ETwo(I),.False.,IPBC,
     $          CelVec,JJ,PA,XX,XX,PB,XX,XX,JJ,FA,FB,XX,XX,IGWInf,
     $          IAtBtD,IRdBtD,RRdBtD,RGWBtD,Omega,XDMA,LenXCU,XCParU,
     $          V(IV),MDV1)
              If(IPrint.gt.0) Write(IOut,1030) I, Ex(I), Ec(I), E6(I),
     $          ETwo(I)
              If(ITyPMM.ne.0.or.PCM(ISCRF)) then
                Call FileIO(2,-IRwSvP,NTTX,PA,Ind)
                If(HaveB) Call FileIO(2,-IRwSvP,NTTX,PB,Ind+NTTM)
                Call RIComp(.True.,HaveB,.False.,.False.,.True.,IOpCl,
     $            .False.,NBasis,NMtPBC,1,1,1,0,NMatY,NMatSY,NMatTY,PA,
     $            PB,V(IV),MDV1)
                endIf
              endIf
            If(IfPE(ITyPMM)) then
              Call PEFock(IOut,IPrint,.True.,.True.,.False.,LPMAll,
     $          PMMAll,NAtoms,IAn,IAtTyp,AtmChg,C,NDAtP,IAtPar,RAtPar,
     $          Conect,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,
     $          AccDes,IOpClY,NBasis,NBas6D,1,1,0,.False.,.True.,.True.,
     $          .False.,PA,PB,FA,FB,EPMM(I),PMMA,.True.,EField(1),V(IV),
     $          MDV1)
            else if(IfMMP(ITyPMM)) then
              Call MMPQM(IOut,IPrint,IOpClY,Init,.True.,.True.,.False.,
     $          IPFlag,AllowP,IHMeth,NAtoms,NBasis,NBas6D,ErfP(I),DoneQ,
     $          IAn,IAtTyp,C,PMMAll,LPMAll,PA,PB,FA,FB,V(IV),MDV1)
              endIf
            If(PCM(ISCRF)) Call PCMQM(IOut,IPrint,IOpClY,Init,.True.,
     $        .True.,.False.,Init,NOpUse,IPFlag,AllowP,FMM,FMFlag,
     $        FMFlg1,NFxFlg,IHMeth,NAE,NBE,NAtoms,NAtoms,NBasis,NBas6D,
     $        ICharg,.False.,PB2,PX2,ErfP(I),DoneQ,IAn,IAtTyp,AtmChg,C,
     $        PA,PB,FA,FB,V(IV),MDV1)
            Call RIComp(.False.,HaveB,.True.,.False.,.True.,IOpCl,
     $        .False.,NBasis,NMtPBC,1,1,1,0,NMatY,NMatSY,NMatTY,FA,FB,
     $        V(IV),MDV1)
            Call FileIO(1,-IRwSvF,NTTX,FA,Ind)
            If(HaveB) Call FileIO(1,-IRwSvF,NTTX,FB,Ind+NTTM)
   30       Continue
          Call FileIO(2,-IRwSvF,NTTM,FA,0)
          If(HaveB) then
            Call FileIO(2,IRwSvF,NTTM,FB,0)
          else if(CopyFB) then
            Call AMove(NTTM,FA,FB)
            endIf
          Call FreeRw(IRwSvF)
        else
          If(DoDFT) Call AXDFT(IOut,IPrint,300,IOpCl,IPFlag,IExCor,
     $      N6Func,ICorTp,NOpDFT,JSym2X,ScaDFX,AccXC,IRwPA,0,IRwPB,0,
     $      NBasis,NMatSI,NMatSI,NMatSI,NAtoms,IAn,IAtTyp,AtmChg,C,JJ,
     $      XX,XX,XX,PA,PB,FA,FB,XX,XX,XX,IGWInf,IAtBtD,IRdBtD,RRdBtD,
     $      RGWBtD,Omega,JJ,XX,LenXCU,XCParU,V(IV),MDV1)
          If(ITyPMM.ne.0) Call GauErr('PMM AX NYI in MkFok1.')
          If(PCM(ISCRF)) Call AXPCM(IOut,IPrint,IOpClY,IPFlag,AllowP,
     $      FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,LSEAll,RSEAll,NAtoms,IAn,
     $      IAtTyp,AtmChg,C,NBasis,NBas6D,NMatI,NMatSI,NMatTI,.False.,
     $      .True.,.True.,.False.,0,0,RJunk,XX,PA,PB,FA,FB,V(IV),MDV1)
          If(CopyFB) Call AMove(NTTM,FA,FB)
          endIf
        Call TStamp(1,'After DFT XC')
      else if(DoFull.and.IPrint.gt.0) then
        Write(IOut,1010) (I,JGC,ETwo(I),I=1,NMatI)
        endIf
C
C     Non-PCM SCRF term.
C
      If(Onsagr(ISCRF).or.SolEPt(ISCRF)) then
        If(.not.DoFull)
     $    Call GauErr('Fock derivatives and SCRF in MkFok1 NYI.')
        If(NRI.eq.2) Call GauErr('Complex NYI for Onsager')
        Call FileIO(2,-IRwSvP,NTTM,PA,0)
        If(HaveB) Call FileIO(2,IRwSvP,NTTM,PB,0)
        Call TstCor(NTTX,MDV1,'MkFok1-SolPMt')
        Do 50 I = NMatSI, 1, -1
          IndPF = NTTX*(I-1) + 1
          If(ClShel(IOpCl)) Call AScale(NTTX,Two,PA(IndPF),PA(IndPF))
          If(Onsagr(ISCRF)) then
            ISCRF0 = Mod(ISCRF,10000)/1000
            ISCRFT = Mod(ISCRF,1000)/100
            If(ISCRF0.eq.0) ISCRFT = ISCRFT*100
            Call Onsag2(IOut,IPrint,IOpClY,ISCRFT,1,NAtoms,AtmChg,C,
     $        NBasis,1,NTTX,gFac,fuFac,ErfP(I),PA(IndPF),PB(IndPF),
     $        FA(IndPF),FB(IndPF),V(IV),EField)
          else
            Call SolPMt(IOut,IPrint,IOpClY,1,IRwH,NBasis,1,NTTX,ErfP(I),
     $        PA(IndPF),PB(IndPF),FA(IndPF),FB(IndPF),V(IV))
            endIf
          If(ClShel(IOpCl)) Call AScale(NTTX,Half,PA(IndPF),PA(IndPF))
   50     Continue
        If(CopyFB) Call AMove(NTTM,FA,FB)
        Call AScale(NMatI,-Half,ErfP,ErfP)
        endIf
C
C     Compute and print the final energy
C
      If(DoFull) then
        Call AAdd(NMatI,ETwo,ECDFit,ETwo)
        Call AAdd(NMatI,EOne,ETwo,Energy)
        Call AAdd(NMatI,Energy,EPMM,Energy)
        Call AAdd(NMatI,Energy,ErfP,Energy)
        EElec = Energy(1)
        If(IPrint.gt.0) Write(IOut,1000)
     $    (I,EOne(I),ETwo(I),Energy(I),I=1,NMatI)
        endIf
      Call FreeRw(IRwSvP)
      Return
      End
