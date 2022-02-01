*Deck SCFDIS
      Subroutine SCFDIS(IOut,IPrint,IDITyp,IFSel,MaxCom,NIter,NumFS,
     $  IOpCl,NBasis,NBsUse,NBsNEO,NBuNEO,NVQ,IRwPA,IRwPXA,IRwFA,IRwFXA,
     $  IRwPB,IRwPXB,IRwFB,IRwFXB,IRwS,IRwV,IRwNED,IRwNEF,IRwNES,IRwNEV,
     $  IRwErr,IRwErN,IRwErQ,NSaved,IOB,IOPA,IOPB,IOFA,IOFB,
     $  IOErr,IOPNA,IOFNA,IOErrN,UCofMx,A,C,D,AN,CN,DN,BMat,Scr,ErrMax,
     $  IndFok,Energy,EnCyc,ErrCyc,WFeI,WQeI,EWFeI,EWQeI,Scr3,LenSc3)
      Implicit Real*8(A-H,O-Z)
C
C     DIIS extrapolation for SCF.  See
C     Pulay, J. Comp. Chem. Vol. 3 No. 4 556-560 (1982).
C
C     IOut   ... output unit.
C     IPrint ... Print flag:
C                0 ... print nothing.
C                1 ... print a few messages.
C                2 ... print B-matrix.
C                3 ... print everything
C     IDITyp ...   0 ... Default (104).
C                  1 ... Commutator DIIS only.
C                  2 ... Energy DIIS only.
C                  3 ... Energy DIIS when error is high, switching to commutator
C                        when error is low and diagonal BMat is decreasing.
C                  4 ... Same as 3 but with linear combination for switchover.
C                  5 ... Energy DIIS only when energy has increased in last cycle.
C                 1x ... Use energy DIIS when commutator gives huge coefficients.
C                Nxx ... Switch to commutator when error is 10^(-N)
C     IFSel  ... Selection of the next Fock matrix:
C                0 ... Normal DIIS
C                1 ... Update DIIS information but use the raw Fock matrix
C                      with the lowest energy (saved in EnCyc).
C                2 ... Update DIIS information but use the raw Fock matrix
C                      with the lowest error (saved in ErrCyc).
C                3 ... Update DIIS information but use the raw Fock matrix
C                      from the current iteration.
C                4 ... Use the lowest energy Raw Fock matrix so far and discard
C                      the others.
C                5 ... Keep only the current point and the lowest energy
C                      preceeding point.
C     MaxCom ... Max number of fock matrices to save.
C     NIter  ... Number of current iteration, -1 to restart, otherwise just
C                used in debug printing.
C     NumFS  ... Number of Fock matrices in IRwFA and IRwFB.  Only the first
C                is used to compute the error but all are extrapolated, unless
C                IRwFXA,B are non-zero, in which case the single matrices in
C                these files are used to compute the error while all matrices
C                in IRwFA,B are extrapolated.
C     NVQ    ... Number of variables to extrapolate in addition to NTT. The
C                corresponding error vector is in IRwErQ.
C     IRwPA  ... RWF on which current alpha density is stored
C                (only used if IRwErr.eq.0).
C     IRwFA  ... RWF on which current alpha Fock matrix is stored.
C     IRwFXA ... Fock matrix to extrapolate.  If 0, the first matrix in IRwFA
C                is used.
C     IRwPB  ... RWF on which current beta density is stored, 0 if RHF.
C     IRwFB  ... RWF on which current beta Fock matrix is stored,
C                0 if RHF.
C     IRwFXB ... Fock matrix to extrapolate.  If 0, the first matrix in IRwFB
C                is used.
C     IRwS   ... RWF of overlap matrix.  0 if calculation is in an
C                orthonormal basis.
C     IRwV   ... RWF of S**-1/2.  0 if calculation is in an
C                orthonormal basis.
C     IRwErr ... Optional rwf of error matrix in AO basis.  If this is
C                non-zero, this is used instead of [F,P].  The matrix
C                is stored lower triangular, with beta following alpha
C                if uhf.
C     NSaved ... Preserves number of Fock matrices currently stores;
C                Initialized here and should not be modified by caller
C     IOB,IOFA,IOFB,IOErr ... Preserve files; should not be altered by
C                caller.
C     UCofMx ... User-specified value for CofMax (see below).  If Zero,
C                a default is applied.
C     A      ... N*(N+1) scratch matrix, where N=Max(NBasis,MaxCom+1)
C     C      ... N*N scratch matrix, where N=Max(NBasis,MaxCom+1)
C     D      ... Max(NBasis,(MaxCom+1)**2) scratch matrix
C     BMat   ... (MaxCom+1)**2 scratch matrix.
C     Scr    ... 7*(MaxCom+1) scratch matrix.  The coefficients of the
C                NSaved Fock matrices being extrapolated are returned in
C                Scr(I+1,1),I=1,NSaved.  The new Fock matrix is the last.
C     ErrMax ... Returned maximum error in [F,P].
C     Energy ... Energy corresponding to new Fock matrix being provided.
C     EnCyc  ... Queue of saved energies.
C     ErrCyc ... Queue of saved errors.
C     IndFok ... Pointers to saved Fock matrices on disk.
C     WFeI   ... weight of error on F for commutator DIIS.
C     WQeI   ... weight of error on Q for commutator DIIS.
C     EWFeI  ... weight of error on F for energy DIIS.
C     EWQeI  ... weight of error on Q for energy DIIS.
C
C     The Fock matrix rwfs are overwritten with the extrapolated Fock
C     matrices.  The parameter CofMax determines when extrapolation
C     coeficients are considered to be too large (and hence extrapolated
C     density matrices too far from idempotent).  If any coefficients
C     exceed CofMax, the extrapolation is still performed but all Fock
C     matrices except the last are skipped in later extrapolations.  If
C     any coefficients exceed CofMx1, no extrapolation is done, and all
C     but the last Fock matrix is discarded.
C
      Parameter (InvMth=1,NDet1=15,MinPrM=2,MinPrE=3)
      Real*8 MDCutO
      Logical Inv2, OK, EOnBig, UseEnD, DidCom, DidEn, CBfn, RRHF,
     $  Spinor, Inv3, SaveP, DoExtP
      Dimension A(*), C(NBasis,*), D(NBasis,*),AN(*),CN(NBsNEO,*),
     $  DN(NBsNEO,*),BMat(MaxCom+1,*),
     $  Scr(MaxCom+1,7), IndFok(MaxCom), EnCyc(MaxCom), ErrCyc(MaxCom),
     $  Scr3(*)
      Save Det, Pt5, One, Two, DCofMx, Zero, DCofMM, BMFact, BMLim,
     $  DetLim, Ten, WtEnMx, SmlSum
      Data Pt5/0.5d0/, One/1.0d0/, Two/2.0d0/, DCofMx/10.0d0/,
     $  Zero/0.0d0/, Det/0.0d0/, Big/1.d10/, DCofMM/10.0d0/,
     $  BMFact/1.1d0/, BMLim/1.d-10/, EMaxC1/1.d-4/, DetLim/1.d-14/,
     $  Ten/10.0d0/, WtEnMx/1.0d0/, SmlSum/1.0d-6/
 1010 Format(' Coeff-Com: ',6D10.3)
 1020 Format(' Inversion failed.  Reducing to',I3,' matrices.')
 1030 Format(' DIIS: error=',1PD9.2,' at cycle',I4,' NSaved=',I4,'.')
 1035 Format(' DIIS: error=',1PD9.2,1x,1PD9.2,' at cycle',I4,' NSaved=',
     $  I4,'.')
 1040 Format(' Large coefficients: NSaved=',I3,' BigCof=',F8.2,
     $  ' CofMax=',F8.2,' Det=',1PD9.2)
 1050 Format(' Use energy DIIS, commutator BigCof=',F8.2,' CofMax=',
     $  F8.2,' Det=',1PD9.2)
 1060 Format(' Coeff-En:  ',6D10.3)
 1070 Format(' NSaved=',I2,' IEnMin=',I2,' EnMin=',G22.15,' IErMin=',I2,
     $  ' ErrMin=',1PD9.2,/,' ErrMax=',1PD9.2,1x,1PD9.2,' EMaxC=',
     $  1PD9.2,' BMatC=',1PD9.2,' BMatP=',1PD9.2,/,' IDIUse=',I1,
     $  ' WtCom=',1PD9.2,' WtEn=',1PD9.2)
 1080 Format(' EnCoef returned ',G17.10,100F12.8)
 1090 Format(' Coeff:     ',6D10.3)
 1100 Format(' Coeff-sum=',F12.10,' Err=',1PD9.2,' no extrapolation.')
      Write(*,*) 'ENTER STMDIS'
      Write(*,*) 'NSAVEDSTART', NSaved
      Write(*,*) 'IOB',IOB
C
C     CBfn whether basis fxn / integrals are complex
       
      Call CkSCFT(IOut,'FockUt',IOpCl,.True.,.True.,.False.,.True.,
     $  NSpBlk,NRI,NDimBl,CBfn,NSpBlX,NSpBNX,Spinor)
      If(CBfn) Call GauErr('SCFDIS with CBfn.')
      NIterA = Abs(NIter)
      
C     MDCutO: Machine dependent small value compared to 1      
      Small = MDCutO(0)
      RRHF = NRI.eq.1.and.NDimBl.eq.1.and.NSpBlk.eq.1.and..not.CBfn
      NBD = NBasis*NDimBl
      NBUD = NBsUse*NDimBl
      NTT = (NBD*(NBD+1))/2
      NTTR = NRI*NTT
      NEONTT = NRI*(NBsNEO*(NBsNEO+1))/2 
      
      UseEnD = .False.
      
C     IOpClX ... Can indicate whether matrix is real or complex 
      IOpClX = 2*(NRI-1)

C     MaxCom ... Max number of fock matrices to save.          
      MaxCm1 = MaxCom + 1
      MxCmSq = MaxCm1**2

C     IDITyp ... Type of DIIS. Default (104).      
      IDITyX = IDITyp
      If(IDITyX.eq.0) IDITyX = 104
      IDITy0 = Max(Mod(IDITyX,10),1)
      IDITy1 = Mod(IDITyX,100)/10
      IDITy2 = Mod(IDITyX,1000)/100
      EOnBig = IDITy0.eq.3.and.IDITy1.eq.1

C     GFloat: Integer to working-precision conversion.
      EMaxC = GFloat(10)**(-IDITy2)
 
      Write(*,*) "AL 1"
C     UCofMx ... User-specified value for CofMax (see below).  If Zero,
C                a default is applied.
C     CofMax ... The parameter CofMax determines when extrapolation
C     coeficients are considered to be too large (and hence extrapolated
C     density matrices too far from idempotent).  If any coefficients
C     exceed CofMax, the extrapolation is still performed but all Fock
C     matrices except the last are skipped in later extrapolations.  If
C     any coefficients exceed CofMx1, no extrapolation is done, and all
C     but the last Fock matrix is discarded.
      NDet = 15
      If(UCofMx.ne.Zero) then
        CofMax = UCofMx
      else if(IDITy0.eq.1) then
C       DCofMx/10.0d0/ 
        CofMax = DCofMx
      else
C       DCofMM/10.0d0
        CofMax = DCofMM
        endIf

      If(IDITy0.le.1) then
        NB = 1
      else
        NB = 3
        endIf
C
C     If this is the first iteration, initialize.  NIter=-1 means
C     restart.
C

C     NVQ    ... Number of variables to extrapolate in addition to NTT. The
C                corresponding error vector is in IRwErQ.
      DoExtP = NVQ.gt.0
      SaveP = DoExtP.or.IDITy0.ge.2
      Write(*,*)'NITERRRR',NIter
      If(NIter.le.1) then
C       If it's the first iteration, we should not have more than NTT variables
C       to extrapolate and the # of fock matrices saved should be 1
        If(NVQ.gt.0.and.NumFS.ne.1)
     $    Call GauErr('NVQ>0 and NumFS>1 in SCFDIS.')
C       NSaved ... Preserves number of Fock matrices currently stores;
C                Initialized here and should not be modified by caller
        NSaved = 0
        LenF = MaxCom*(NTTR+NVQ)
        LenP = MaxCom*(NTTR+NVQ)
  
        If(NBsNEO.gt.0) then
          LenFN = MaxCom*NEONTT
          LenPN = MaxCom*NEONTT
          endif



C       EnCyc  ... Queue of saved energies.
C       ErrCyc ... Queue of saved errors.
C       Big/1.d10/       
        EnCyc(1) = Big
        ErrCyc(1) = Big

C       WtEnMx/1.0d0/     
        WtEnMx = One

C       MaxCom ... Max number of fock matrices to save.          
        If(MaxCom.lt.2) Call GauErr('Illegal MaxCom in SCFDIS.')

          Write(*,*) "AL 1.4"
          Write(*,*) IOB
C       IOB,IOFA,IOFB,IOErr ... Preserve files; should not be altered by
C       caller.
        If(IOB.eq.0) then
C         Integer Function IGetRW(LEN, lUnitl):
C         Allocates a scratch read-write file of length Abs(LEN) on unit
C         lUnitl. If LEN is positive, this is a routine-volatile file,
C         while if LEN is negative, it is link-volatile. If lUnitl is 0,
C         the file is allocated on the standard scratch unit lUnit(lO).
          IOB = IGetRW(NB*MxCmSq,0)
          IOPA = 0
          Write(*,*) "AL 1.45"
          Write(*,*) IOB

C         NumFS  ... Number of Fock matrices in IRwFA and IRwFB.  Only the first
C                    is used to compute the error but all are extrapolated, unless
C                    IRwFXA,B are non-zero, in which case the single matrices in
C                    these files are used to compute the error while all matrices
C                    in IRwFA,B are extrapolated.
          If(SaveP) IOPA = IGetRW(NumFS*LenP,0)
          Write(*,*) "AL 1.46"
          Write(*,*) IOFA
          Write(*,*) 'SIZE OF IOFA',NUMFS*LenF
          IOFA = IGetRW(NumFS*LenF,0)

          Write(*,*) "AL 1.5"
          Write(*,*) IOFA
          Write(*,*) 'SIZE OF IOFA',NUMFS*LenF
          IOPB = 0
          IOFB = 0
          If(NSpBlk.eq.2) then
            If(SaveP) IOPB = IGetRW(NumFS*LenP,0)
            IOFB = IGetRW(NumFS*LenF,0)
            endIf
          LenErr = MaxCom*(NSpBlk*NTTR+NVQ)
          IOErr = IGetRW(LenErr,0)

          If(NBsNEO.gt.0) then
            IOFNA = IGetRW(NumFS*LenFN,0)
            IOPNA = 0
            If(SaveP) IOPNA = IGetRW(NumFS*LenPN,0)
            IOErrN = IGetRW(MaxCom*NEONTT,0)
            endIf

          endif

        Write(*,*) "AL 2"
        Write(*,*) "AL 2.1"
        Write(*,*) BMat(1,1)
C       AClear(N, A):
C       Clears N elements in array A.
C       BMat   ... (MaxCom+1)**2 scratch matrix.           
        Call AClear(MaxCm1**2,BMat)
        Write(*,*)'MaxComBeforeLoop',MaxCom
        Do 10 I = 1, MaxCom
C         IndFok ... Pointers to saved Fock matrices on disk.
          IndFok(I) = I - 1
          BMat(1,I+1) = -One
   10     BMat(I+1,1) = -One
      
        Write(*,*)'MaxComAfterLoop',MaxCom
        Write(*,*)'Print INDFOK1'
        Do 12 I = 1, MaxCom
   12     Write(*,*) 'INDFOK',IndFok(I)
        
        Do 20 IB = 0, (NB-1)
          Call FileIO(1,-IOB,MxCmSq,BMat,MxCmSq*IB)
   20     Continue

        Call OutMtS(IOut,'INITIALIZE BMAT',0,0,BMat,MaxCm1,MaxCM1,MaxCM1
     $    ,MaxCm1)
        endIf
C     Ending the If (NITer < 1) condition 

      Write(*,*) "AL 3" 
C     Form error matrices.  NSaved is the number of Fock
C     matrices actually available in the queue.  To free space, we replace
C     the Fock matrix corresponding to the highest energy in the queue,
C     except that we always keep the results of the previous iteration.
C
C     GetWE(WeI,WeD,We):
C     Little function to set/default We to WeI if the latter is non-zero
C     or to the default WeD if WeI is negative.
C     WFeI   ... weight of error on F for commutator DIIS.
C     WQeI   ... weight of error on Q for commutator DIIS.
      Call GetWE(WFeI,One,WFe)
      Call GetWE(WQeI,One,WQe)
      NSaved = NSaved + 1
      Write(*,*) "NSaved+1",NSaved
      MaxCom = MaxCm1 -1 
      Write(*,*) "MaxCom",MaxCom

      If(NSaved.gt.MaxCom) Call GauErr('Inconsistency #1 in SCFDIS.')
      IndF = IndFok(NSaved)
C     IRwErr ... Optional rwf of error matrix in AO basis.  If this is
C                non-zero, this is used instead of [F,P].  The matrix
C                is stored lower triangular, with beta following alpha
C                if uhf.
      Write(*,*) "AL 4" 
      If(IRwErr.ne.0) Call FileIO(2,-IRwErr,0,A,0)
      If(IRwErN.ne.0) Call FileIO(2,-IRwErN,0,AN,0)
      Write(*,*) "AL 5" 
      
      NTTU = (NBUD*(NBUD+1))/2
      NTTUR = NRI*NTTU

      NEONTU = (NBuNEO*(NBUNEO+1))/2
      NENTUR = NRI* NEONTU

C     Save the Fock matrix at offset IndF in IOF.
      LenP = NumFS*NTTR + NVQ
      LenF = NumFS*NTTR + NVQ
      

      Write(*,*)'Print INDFOK'
      Do 11 I = 1, MaxCom
   11   Write(*,*) 'INDFOK',IndFok(I)

      Write(*,*) "IOFA",IOFA
      Write(*,*) "IRwFA",IRwFA
      Write(*,*) 'NumFS',NumFS
      Write(*,*) 'IndF',IndF
      Write(*,*) 'LenF',LenF
      Call FileIO(2,-IRwFA,LenF,A,0)
      Write(*,*) 'POSITIONTOWRITEINIOFA', IndF*LenF
      Call FileIO(1,-IOFA,LenF,A,IndF*LenF)
      Write(*,*) 'Before FPSSPF'
      
      If(NBsNEO.gt.0) then
         LenPN = NumFS*NEONTT
         LenFN = NumFS*NEONTT
         Call FileIO(2,-IRwNEF,LenFN,AN,0)
         Call FileIO(1,-IOFNA,LenFN,AN,IndF*LenFN)
         endif
C     Subroutine FPSSPF(NRI,NBasis,NBsUse,IRwP,IRwF,IRwFX,IRwS,IRwV,
C    $  IRwErr,NumFS,A,C,D,Scr,LenScr)
C     Implicit Real*8(A-H,O-Z)
C
C     Form the commutator of F and P (FPS-SPF).  This is returned
C     as a lower triangle (including diagonal) in A.  A is 2*NTT,
C     C is NBasis**2 and D is NBasis.  IRwV is used to
C     transform the commutator to an orthonormal basis.  If IRwS
C     is 0, the basis is assumed to be orthonormal and IRwS and IRwV
C     are not used. If IRwErr is non-zero, the lower triangular error
C     matrix is taken from it.
C
      Call FPSSPF(NRI,NBD,NBUD,IRwPA,IRwFA,IRwFXA,IRwS,IRwV,IRwErr,
     $  NumFS,A,C,D,Scr3,LenSc3)
      If(NBsNEO.gt.0) then
        Call FPSSPF(NRI,NBsNEO,NBuNEO,IRwNED,IRwNEF,0,IRwNES,IRwNEV,
     $    IRwErN,NumFS,AN,CN,DN,Scr3,LenSc3)
        endif

      If(RRHF) Call AScale(NTTUR,Pt5,A,A)
 
      Write(*,*) 'After FPSSPF'
      Write(*,*) 'RRRRHF',RRHF

C     Save A (which contains [F,P]) into IOErr with offset IOffEr
      IOffEr = IndF*(NSpBlk*NTTUR+NVQ)
      Call FileIO(1,-IOErr,NTTUR,A,IOffEr)
      If(IPrint.ge.MinPrE) Call LTOutX(IOut,0,IOpClX,1,NIterA-1,
     $  'Alpha error matrix at iteration',NBUD,NBUD,A,A,0)

      If(NBsNEO.gt.0) then
        IOfErN = IndF*NEONTT
        Call FileIO(1,-IOErrN,NEONTT,AN,IOfErN)
        If(IPrint.ge.MinPrE) Call LTOutX(IOut,0,IOpClX,1,NIterA-1,
     $    'NEO error matrix at iteration',NBuNEO,NBuNEO,AN,AN,0)
        endif


C     Function ArMaxC(NRI, A, N) 
C     Return the maximum magnitude of an element of A, which may be real or complex.
      ErrMax = ArMaxC(NRI,A,NTTU)
      If(NSpBlk.ne.1) then
        Call FileIO(2,-IRwFB,LenF,A,0)
        Call FileIO(1,-IOFB,LenF,A,IndF*LenF)
        Call FPSSPF(NRI,NBD,NBUD,IRwPB,IRwFB,IRwFXB,IRwS,IRwV,IRwErr,
     $    NumFS,A,C,D,Scr3,LenSc3)
        Call FileIO(1,IOErr,NTTUR,A,0)
        If(IPrint.ge.MinPrE) Call LTOutX(IOut,0,IOpClX,1,NIterA-1,
     $    'Beta error matrix at iteration',NBUD,NBUD,A,A,0)
        X = ArMaxC(NRI,A,NTTU)
        ErrMax = Max(ErrMax,X)
        endIf
      If(.not.RRHF) ErrMax = ErrMax*Pt5
      


C     Add the error for the NVQ variables.
      If(NVQ.gt.0) then
        Call FileIO(2,-IRwErQ,NVQ,A,0)
C       Function ARRMax(A, N, IFABS) 
C       This function returns the maximum element of an array. Ifabs is .true, if the absolute values
C       are to be compared, .false, if the signs are significant.
        ErrMxQ = ArrMax(A,NVQ,.True.)
        Call FileIO(1,IOErr,NVQ,A,0)
        If(IPrint.ge.MinPrE) Call OutMtS(IOut,
     $    'Q error vector at iteration',NIterA-1,0,A,NVQ,1,NVQ,1)
      else
        ErrMxQ = Zero
        endIf
      ErrTot = WFe*ErrMax + WQe*ErrMxQ



C     Save densities
      If(SaveP) then
        Call FileIO(2,-IRwPXA,LenP,A,0)
        Call FileIO(1,-IOPA,LenP,A,IndF*LenP)
        If(NSpBlk.ne.1) then
          Call FileIO(2,-IRwPXB,LenP,A,0)
          Call FileIO(1,-IOPB,LenP,A,IndF*LenP)
          endIf
        endIf
      If(IPrint.gt.0) then
        If(NVQ.eq.0) then
          Write(IOut,1030) ErrMax, NIterA, NSaved
        else
          Write(IOut,1035) ErrMax, ErrMxQ, NIterA, NSaved
          endIf
        endIf

      Write(*,*)'AODONG1'
      Write(*,*)'NSaved',NSaved

C     IFSel  ... Selection of the next Fock matrix:
C                0 ... Normal DIIS
C                1 ... Update DIIS information but use the raw Fock matrix
C                      with the lowest energy (saved in EnCyc).
C                2 ... Update DIIS information but use the raw Fock matrix
C                      with the lowest error (saved in ErrCyc).
C                3 ... Update DIIS information but use the raw Fock matrix
C                      from the current iteration.
C                4 ... Use the lowest energy Raw Fock matrix so far and discard
C                      the others.
C                5 ... Keep only the current point and the lowest energy
C                      preceeding point.
      If(IFSel.eq.4.or.IFSel.eq.5) then
        EnCyc(NSaved) = Energy
        NDel = NSaved - IFSel + 3
        Do 100 I = 1, NDel
C         Subroutine DISRmX(NB,MaxCom,NLook,IndFok,IOB,EnCyc,ErrCyc,Scr,
C         $  ISkip,BMat,C)
C         Pick a saved Fock matrix to remove from the DIIS data, updating then
C         three arrays and the B-matrices and returning ISkip for later use.
C         The updated commutator B-matrix is left in BMat.
          Call DISRmX(NB,MaxCom,NSaved,IndFok,IOB,EnCyc,ErrCyc,A,ISkip,
     $      BMat,C)
  100     NSaved = NSaved - 1
        endIf
      Write(*,*)'AODONG2'
      Write(*,*)'NSaved',NSaved
C
C     Update commutator B-matrix.
C
      Write(*,*) 'MAXCOMMM',MaxCom
      Write(*,*) 'MXCMSQQQ',MxCmSq
      Call FileIO(2,-IOB,MxCmSq,BMat,0)
      Fact = -One
Cal     If(RRHF) Fact = Fact*Two
      FactQ = One
      If(.not.RRHF) FactQ = FactQ*Two*Two

      If(NBsNEO.gt.0) then
        Call UpdtB2(IOpCl,NSaved,MaxCom,NBsUse,NTTUR,NBuNEO,NENTUR,NVQ,
     $    One,-FactQ*WQe,IOErr,IOErrN,BMat,A,C,AN,CN,Energy,EnCyc,
     $    ErrTot,ErrCyc,IndFok)
      else
        Call UpdatB(IOpCl,NSaved,MaxCom,NBsUse,NTTUR,NVQ,Fact*WFe,
     $    -FactQ*WQe,IOErr,BMat,A,C,Energy,EnCyc,ErrTot,ErrCyc,IndFok)
        endif

      Call FileIO(1,-IOB,MxCmSq,BMat,0)
      If(IPrint.ge.MinPrM) Call OutMtS(IOut,'B matrix at iteration',
     $  NIterA,1,BMat,MaxCm1,MaxCm1,NSaved+1,NSaved+1)
C
C     Update the energy-based matrices.
C
      If(IDITy0.ge.2) then
        Call GetWE(EWFeI,One,EWFe)
        Call GetWE(EWQeI,One,EWQe)
        Fact = Pt5
        If(RRHF) Fact = Fact*Two
        IOpClX = 2*(NRI-1)
        Call AClear(MxCmSq,BMat)
C       Subroutine UpdtTr(NSaved,MaxCom,IOpClX,NBasis,NTT,NumFS,NVQ,WFe,
C       $  WQe,BMat,A,B,IndFok,IOF,IOP)
C
C       Compute certain traces to run DIIS with energy.
C
C       BMat(I+1,1) = Tr(Fn*Pi)
C       BMat(I+1,2) = Tr(Fi*Pn)
C       A and B are scratch arrays of length len.
        Call UpdtTr(NSaved,MaxCom,IOpClX,NBD,NTTR,NumFS,NVQ,Fact*EWFe,
     $    EWQe,BMat,A,C,IndFok,IOFA,IOPA)
        If(NSpBlk.eq.2) Call UpdtTr(NSaved,MaxCom,IOpClX,NBD,NTTR,NumFS,
     $    NVQ,Fact*EWFe,EWQe,BMat,A,C,IndFok,IOFB,IOPB)
        Call AMove(2*(MaxCom+1),BMat,C)
        Call FileIO(2,-IOB,MxCmSq,BMat,MxCmSq)
        Call FileIO(2,-IOB,MxCmSq,A,2*MxCmSq)
        Call UpdtBE(NSaved,MaxCom,BMat,A,C)
        Call FileIO(1,-IOB,MxCmSq,BMat,MxCmSq)
        Call FileIO(1,-IOB,MxCmSq,A,2*MxCmSq)
        If(IPrint.ge.MinPrM) then
          Call OutMtS(IOut,'Tr matrix at iteration',NIterA,1,A,MaxCm1,
     $      MaxCm1,NSaved+1,NSaved+1)
          Call OutMtS(IOut,'BE matrix at iteration',NIterA,1,BMat,
     $      MaxCm1,MaxCm1,NSaved+1,NSaved+1)
          endIf
        endIf

      
      Write(*,*)'AODONG3'
      Write(*,*)'NSaved',NSaved

      Call FileIO(2,-IOB,MxCmSq,BMat,0)
      BMatC = Abs(BMat(NSaved+1,NSaved+1))
      BMatP = Abs(BMat(2,2))

C     Find the smallest diagonal value in the B matrix
      Do 110 It = 3, NSaved
  110   BMatP = Min(BMatP,Abs(BMat(It,It)))
      EnMin = ArMin1(EnCyc,NSaved,.False.,IEnMin)
      ErrMin = ArMin1(ErrCyc,NSaved,.False.,IErrMin)
      If(IFSel.ge.1.and.IFSel.le.3) then
        IDIUse = 1
      else if(IDITy0.eq.1.or.IDITy0.eq.2) then
        IDIUse = IDITy0
      else if(IDITy0.eq.3) then
        If(NIter.le.2.or.(ErrMax.lt.EMaxC.and.
     $    (BMatC.lt.BMLim.or.BMatC.le.(BMatP*BMFact)))) then
          IDIUse = 1
        else
          IDIUse = 2
          endIf
      else if(IDITy0.eq.4) then
C       If(NIter.le.2) then
C         IDIUse = 1
C       else if(ErrMax.lt.EMaxC1) then
        If(ErrMax.lt.EMaxC1) then
          IDIUse = 1
        else if(BMatC.gt.(BMatP*BMFact)) then
          If(ErrMax.ge.(EMaxC1*Ten*Ten)) then
            IDIUse = 2
          else

            WtEn = Max(Min(RE2Com/(RE2Com+One),WtEnMx),Zero)
            If(WtEn.eq.Zero) IDIUse = 1
            endIf
        else
          IDIUse = 3
          WtEn = Max(Min(One,ErrMax/EMaxC,WtEnMx),Zero)
          If(WtEn.eq.Zero) IDIUse = 1
          endIf
      else
        IDIUse = 1
        Call GauErr('Illegal IDITyp in SCFDIS.')
        endIf


      Write(*,*)'AODONG4'
      Write(*,*)'NSaved',NSaved
      If(IDIUse.eq.1) then
        WtEn = Zero
      else if(IDIUse.eq.2) then
        WtEn = One
      else if(IDIUse.ne.3) then
        Call GauErr('Consistency failure #2 in SCFDIS.')
        endIf
      If(NIter.ge.5) WtEnMx = Min(WtEnMx,WtEn)
      WtCom = One - WtEn
      If(IPrint.gt.0) Write(IOut,1070) NSaved, IEnMin, EnMin, IErrMin,
     $  ErrMin, ErrMax, ErrMxQ, EMaxC, BMatC, BMatP, IDIUse, WtCom, WtEn
C
C     Invert Pulay's B matrix, using a routine which ensures numerical
C     stability even if the matrix is quite ill-conditioned.  If the
C     matrix is algorithmically singular, go back and eliminate another
C     row and column.
C
      Call AClear(2*MaxCm1,Scr(1,6))
      If(Mod(IDIUse,2).eq.1) then
        Call FileIO(2,-IOB,MxCmSq,BMat,0)
  200   If((IFSel.eq.0.or.IFSel.eq.5).and.NSaved.ge.2.and.
     $    ErrTot.ge.Small) then
          NDim = NSaved + 1
C         Subroutine Pcck(K,KK,X,M,N,MM,NN)
C
C         Pack X dimension (M,N) but used (MM,NN) down to the first
C         words in X.  K=0 to pack a rectangular matrix to compressed form;
C         K=1 to pack a symmetric matrix.  KK is set to the length of the
C         compressed matrix.
          Call Pcck(0,KK,BMat,MaxCm1,MaxCm1,NDim,NDim)
          If(InvMth.eq.2) Call PreCon(0,NDim,NDim,NDim,BMat,Scr(1,7))
          If(InvMth.ge.0.and.InvMth.le.2) then
            IU = 0
            If(IPrint.gt.1) IU = IOut
            IOpt = Min(InvMth,1)

C           Logical Function Inv2(IOut, lOpt, A, N, IS, lADl, IAD2, D, Det, Q, R, QR, LenA)
C           Matrix inversion with iterative refinement. Function value: .true, if inversion was successful.

C           IU   ... Non-zero for print
C           IOpt ... l=in memory refinement: Q, R, and QR are N’*^N.
C           A    ... Matrix to invert. Inverse returned here
C           N    ... Dimension
C           IS   ... 2*N scratch vector
C           IAD1 ... N scratch vector               
C           IAD2 ... N scratch vector
C           D    ... N scratch vector 
C           Det  ... Approximate determinant is returned here if matrix is algorithmically singular.
C           Q    ... N*N scratch matrix if IOpt = 1 Can overlap with IS, lADl, IAD2, D.
C           R    ... N*N scratch matrix if IOpt = 1
C           QR   ... N*N scratch matrix if IOpt = 1 Must not overlap with IS,IAD1,IAD2,D.
C           LenA ... Length of A if IOpt = -1
            OK = Inv2(IU,IOpt,BMat,NDim,Scr,Scr(1,3),Scr(1,4),
     $        Scr(1,5),Det,A,C,D,0)

C           BMat now contains the inverse of B
            OK = OK.and.(Abs(Det).gt.DetLim.or.NSaved.le.NDet1)
            If(OK) then

C             Subroutine AMul(L,A,B,C)
C             Vector C = A*B
              If(InvMth.eq.2) Call AMul(NDim,Scr(1,7),BMat,Scr(1,7))
C             Subroutine ANeg(N,A,B)
C             Vector B = -A.
C             The RHS is vector (-1,0,0,0....), we want to left multiply
C             by B^-1. Thus the coeff. is just the negative of the first
C             column of B^-1
C             We save NDim numbers of the first column of B^-1 into
C             Scr(1,6)
              Call ANeg(NDim,BMat,Scr(1,6))
            else
              Call AClear(NDim,Scr(1,6))
              Scr(NSaved+1,6) = One
              endIf
          else if(InvMth.eq.3) then
            Call AClear(NSaved,Scr(2,6))
            Scr(1,6) = -One
            Call SimEqn(NDim,NDim,BMat,1,1,Scr(1,6),IError)
            OK = IError.eq.0
          else if(InvMth.eq.4) then
            ICntrl = 121
            OK = Inv3(IOut,IPrint,ICntrl,NDim,BMat,Scr3,LenSc3)
            If(OK) then
              Call ANeg(NDim,BMat,Scr(1,6))
            else
              Call AClear(NDim,Scr(1,6))
              Scr(NSaved+1,6) = One
              endIf
          else
            Call GauErr('Unrecognized InvMth.')
            OK = .False.
            endIf

      Write(*,*)'AODONG5'
      Write(*,*)'NSaved',NSaved

          If(OK) then
            BigCof = ArrMax(Scr(2,6),NSaved,.True.)
            OK = BigCof.le.CofMax
          else
            BigCof = Zero
            endIf
          If(.not.OK) then
            If(EOnBig) then
              If(IPrint.gt.0) Write(IOut,1050) BigCof, CofMax, Det
              IDIUse = 2
              WtCom = Zero
              WtEn = One
            else
              If(IPrint.gt.0) Write(IOut,1040) NSaved,BigCof,CofMax,Det
              If(NSaved.gt.2) then
                Call DISRmX(NB,MaxCom,NSaved-2,IndFok,IOB,EnCyc,ErrCyc,
     $            A,ISkip,BMat,C)
                NSaved = NSaved - 1
      Write(*,*)'AODONG6'
      Write(*,*)'NSaved',NSaved
                If(IPrint.gt.0) Write(IOut,1020) NSaved
                Goto 200
                endIf
              endIf
            endIf
        else
          IFX = NSaved
          If(IFSel.eq.1) then
            Fact = ArMin1(EnCyc,NSaved,.False.,IFX)
          else if(IFSel.eq.2) then
            Fact = ArMin1(ErrCyc,NSaved,.False.,IFX)
            endIf
          Scr(IFX+1,6) = One
          endIf
        endIf
 
 
      DidCom = Mod(IDIUse,2).eq.1
      DidEn = Mod(IDIUse,4).ge.2
      If(DidEn) then
        Call FileIO(2,-IOB,MxCmSq,BMat,MxCmSq)
        Call FileIO(2,-IOB,MxCmSq,C,2*MxCmSq)
        Call EnCoef(IOut,IPrint,MaxCom,NSaved,EnCyc,BMat,C,Scr,ISkip,
     $    IWarn)
        If(IWarn.eq.1.and.IDIUse.eq.3) then
          WtEn = Zero
          WtCom = One
        else if(IPrint.ge.MinPrE) then
          Write(IOut,1080) (C(KK,1),KK=1,NSaved+1)
          endIf
        Call FileIO(1,-IOB,MaxCm1,C,2*MxCmSq)
        Call FileIO(2,-IOB,MaxCm1,Scr(1,7),2*MxCmSq)
      else
        ISkip = 0
        endIf
      Write(*,*)'AODONG7'
      Write(*,*)'NSaved',NSaved
C
C     Compute the extrapolated Fock matrix.
C
      If(IPrint.ge.MinPrM.or.(DidCom.and.IPrint.eq.1))
     $  Write(IOut,1010) (Scr(I+1,6),I=1,NSaved)
      If(IPrint.ge.MinPrM.or.(DidEn.and.IPrint.eq.1))
     $  Write(IOut,1060) (Scr(I+1,7),I=1,NSaved)
   
C     Subroutine ACSASB(N,A,SA,B,SB,C)
C     Form the linear combination C(I) = SA*A(I) + SB*B(I).
      
C     Put the combined coefficient at (1,1) position of the Scr matrix
      Call ACSASB(NSaved+1,Scr(1,6),WtCom,Scr(1,7),WtEn,Scr)
      
      If(IPrint.ge.1) Write(IOut,1090) (Scr(I+1,1),I=1,NSaved)
      
C     Function ArrSum(N,A)
C     Sum the elements of array A.
      Sum = ArrSum(NSaved,Scr(2,1))
      SumErr = Abs(Sum-One)
      OK = SumErr.le.SmlSum
      
      Write(*,*)'OKAY?', OK
      If(OK) then
      Write(*,*)'AODONG7.5'
      Write(*,*)'NSaved',NSaved
C
C       Subroutine NewF(IOut,IPrint,N,Len,IRwIn,IRwOut,IndFok,Factor,Coef,
C       $  A,B)
C
C       Compute the vector of length Len which is linear combination
C       of the N vectors in IRwIn (after skipping the first NSkip).
C       The output vector is returned in A, and also written to
C       IRwOut if IRwOut is non-zero.  A and B are both vectors of
C       length Len.
C
C       A = Sum (I=1,N)   IRwIn(I) * Factor * Coef(I)
C
        Call NewF(IOut,IPrint,NSaved,LenF,IOFA,IRwFA,IndFok,One,
     $    Scr(2,1),A,C)

        If(NBsNEO.gt.0) then
          Call NewF(IOut,IPrint,NSaved,LenFN,IOFNA,IRwNEF,IndFok,One,
     $      Scr(2,1),AN,CN)
          endif

        If(IPrint.ge.MinPrE) then
          Call LTOutX(IOut,0,IOpClX,1,NIterA-1,
     $    'Extrapolated FA at iteration',NBD,NBD,A,A,0)
          If(NVQ.gt.0) Call OutMtS(IOut,
     $      'Extrapolated dE/dQ at iteration',NIterA-1,0,
     $      A(NumFS*NTTR+1),NVQ,1,NVQ,1)
          endIf
        If(DoExtP) then
          Call NewF(IOut,IPrint,NSaved,LenP,IOPA,IRwPXA,IndFok,One,
     $      Scr(2,1),A,C)
          If(NBsNEO.gt.0) then
            Call NewF(IOut,IPrint,NSaved,LenFN,IOFNA,IRwNEF,IndFok,One,
     $        Scr(2,1),AN,CN)
            endif
          If(IPrint.ge.MinPrE) then
            Call LTOutX(IOut,0,IOpClX,1,NIterA-1,
     $        'Extrapolated PA at iteration',NBD,NBD,A,A,0)
            If(NVQ.gt.0) Call OutMtS(IOut,
     $        'Extrapolated Q at iteration',NIterA-1,0,
     $        A(NumFS*NTTR+1),NVQ,1,NVQ,1)
            endIf
          endIf
        If(NSpBlk.eq.2) then
          Call NewF(IOut,IPrint,NSaved,LenF,IOFB,IRwFB,IndFok,One,
     $      Scr(2,1),A,C)
          If(IPrint.ge.MinPrE) Call LTOutX(IOut,0,IOpClX,1,NIterA-1,
     $      'Extrapolated FB at iteration',NBD,NBD,A,A,0)
          If(DoExtP) then
            Call NewF(IOut,IPrint,NSaved,LenP,IOPB,IRwPXB,IndFok,One,
     $        Scr(2,1),A,C)
            If(IPrint.ge.MinPrE) Call LTOutX(IOut,0,IOpClX,1,NIterA-1,
     $        'Extrapolated PB at iteration',NBD,NBD,A,A,0)
            endIf
          endIf
      else if(IPrint.ge.MinPrE) then
        Write(IOut,1100) Sum, SumErr
        endIf



      Write(*,*)'AODONG7.6'
      Write(*,*)'NSaved',NSaved
      
      If(NSaved.eq.MaxCom) then
        If(DidEn.and.ISkip.ne.0.and.(UseEnD.or..not.DidCom)) then
          NTest = 0
        else
          NTest = NSaved - 1
      Write(*,*)'AODONG8'
      Write(*,*)'NSaved',NSaved
          endIf

C       Subroutine DISRmX(NB,MaxCom,NLook,IndFok,IOB,EnCyc,ErrCyc,Scr,
C      $  ISkip,BMat,C)
C
C       Pick a saved Fock matrix to remove from the DIIS data, updating then
C       three arrays and the B-matrices and returning ISkip for later use.
C       The updated commutator B-matrix is left in BMat.
        Call DISRmX(NB,MaxCom,NTest,IndFok,IOB,EnCyc,ErrCyc,A,ISkip,
     $    BMat,C)
        NSaved = NSaved - 1
      Write(*,*)'AODONG9'
      Write(*,*)'NSaved',NSaved
        endIf
      Write(*,*) "NSavedEND", NSaved
      Write(*,*) 'LEAVE SCFDIS'
      Return
      End
