*Deck ConOpn
      Subroutine ConOpn(IOut,IOpCl,N,NU,RMSDP,IFlag,IExtp,A1,A2,A3,
     $  Crit,CritM,CritE,DelEn,DPMax,ICount,Done,IRwPA,IRwPB,IRwC1,
     $  IRwC2,IRwC3,ETol,IPass,NMtPBC,Damp,DampFc)
      Implicit Real*8(A-H,O-Z)
C
C     Simplified SCF convergence testing.
C
      Dimension A1(*), A2(*), A3(*), XX(1)
      Character*30 ExtMsg(4)
      Logical Done, Damp, Do34, CBfn, Spinor
      Real*8 MDCutO
      Save Zero, Pt99, Pt995, One, OnePt9, Two, Four, SP12, SP22,
     $  ExtMsg, XX
      Data Zero/0.0d0/, Pt99/0.99d0/, Pt995/0.995d0/, One/1.0d0/,
     $  OnePt9/1.9d0/, Two/2.0d0/, Four/4.0d0/, SP12/0.0d0/,
     $  SP22/0.0d0/, ExtMsg/'3-Point extrapolation',
     $  '4-Point extrapolation','3-Point extrapolation skipped',
     $  '4-Point extrapolation skipped'/, XX/0.0d0/
 2000 Format(' Damping current iteration by',1PD9.2)
 2010 Format(1X,A,'.')
C
      Write(*,*) "ENTERCONOPN"
      Small = MDCutO(1)
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      IOpClX = IOpCl - Mod(IOpCl,2)
      NMat = NSpBlk*NMtPBC
      NB = N*NDimBl
      NTT = (NB*(NB+1))/2
      IType = 0
      IFlag = 0
      Done = .False.
      ROp = GFloat(3-NSpBlk)
      NTTX = NMtPBC*NTT
      NTT2 = NSpBlk*NTTX
      NTTXR = NTTX*NRI
      NTT2R = NTT2*NRI
      ICount = ICount + 1
      Write(*,*) 'ICount',ICount 
      Call FileIO(2,-IRwPA,NTTXR,A3,0)
      Write(*,*) 'Reading the newest density matrix into A3'
      If(NSpBlk.eq.2) Call FileIO(2,-IRwPB,NTTXR,A3(1+NTTXR),0)
      Call LTOutN(IOut,'Newest Density',0,1,N,A3,0)
      
      If(ICount.gt.1) then
        Write(*,*) 'Reading the last density in IRwC1 into A1'
        Call FileIO(2,-IRwC1,NTT2R,A1,0)
        Call LTOutN(IOut,'Last Density',0,1,N,A1,0)
        endif
      Call FileIO(1,-IRwC1,NTT2R,A3,0)
      Write(*,*) 'Writing the newest density matrix into A3'
C     Nothing else to do in first cycle or immediately after extrapolation
      If(ICount.eq.1) Return
C     P(N) now in A3 and P(N-1) in A1; form P(N)-P(N-1) in A2.a
      Write(*,*) 'The difference is saved into A2'
      Call ASub(NTT2R,A3,A1,A2)
      Call LTOutN(IOut,'Difference in Density',0,1,N,A2,0)
      DPMax = ArMaxC(NRI,A2,NTT2)
C     Find length DP1
      SP11 = ROp*SCFTrX(.False.,.False.,IOpClX,NMat,N,A2,XX,A2,XX)
      DP1 = Sqrt(SP11/Two)
C     Test for convergence by finding if root mean square dp is less
C     than Acurcy
      RMSDP = DP1 / (GFloat(NU)*Sqrt(GFloat(NDimBl*NMtPBC)))
C     Convergence criterion met ... exit
Cal

      Write(*,*) 'RMSDP', RMSDP
      Write(*,*) 'Crit',  Crit
      Write(*,*) 'DPMax', DPMax
      Write(*,*) 'CritM', CritM
      Write(*,*) 'Abs(DelEn)', Abs(DelEn)
      Write(*,*) 'ETol', ETol
C     Modified the convergence criterion      
      If(RMSDP.lt.Crit.and.DPMax.lt.CritM.and.
     $  Abs(DelEn).lt.ETol) then
        Done = .True.
        Return
        endIf

      Done = .False.
      Return
      End
