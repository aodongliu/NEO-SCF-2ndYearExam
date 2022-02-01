*Deck UpdtB2
      Subroutine UpdtB2(IOpCl,NSaved,MaxCom,NBsUse,NTTUR,NBuNEO,NENTUR,
     $   NVQ,WFe,WQe,IOErrE,IOErrP,BMat,AElec,BElec,AProt,BProt,Energy,
     $   EnCyc,ErrMax,ErrCyc,IndFok)
      Implicit Real*8(A-H,O-Z)
C
C     Update the B matrix for DIIS:
C
C     BMat(NSaved+1,I+1) += Error(I).Error(NSaved) * Factor
C
C     A and B are scratch arrays of length len.
C
      Logical CBfn, Spinor
      Dimension AElec(*),BElec(*),BMat(MaxCom+1,MaxCom+1),EnCyc(MaxCom),
     $  ErrCyc(MaxCom),IndFok(MaxCom)
      Save One
      Data One/1.0d0/
C
      Call CkSCFT(IOut,'FockUt',IOpCl,.True.,.True.,.False.,.True.,
     $  NSpBlk,NRI,NDimBl,CBfn,NSpBlX,NSpBNX,Spinor)
      ErrCyc(NSaved) = ErrMax
      LenErE = NSpBlk*NTTUR + NVQ
      LenErP = NENTUR
      IA = 1
      IB = NTTUR + 1
      IQ = NSpBlk*NTTUR + 1
      
C     Load the newly formed error vector into AElec and AProt       
      Call FileIO(2,-IOErrE,LenErE,AElec,IndFok(NSaved)*LenErE)
      Call FileIO(2,-IOErrP,LenErP,AProt,IndFok(NSaved)*LenErP)
      
C     Function SCFTrX(Square, Transp, lOpCl, NBasis, PA, PB, FA, FB):
C     Trace the product of one or two sets of matrices, square or triangular
C     real or complex, possibly over 2 or 4 spinors.  For complex matrices,
C     only the real part of the trace of the product is returned.
C
C     Transp is true to do Tr(A.Bt) (=Tr(At.B) for square matrices)
C     instead of Tr(A.B), and is true to do Tr(A*.B)=Tr(A.B*) rather than
C     Tr(A.B) for complex Hermetian matrices.

C     Calculate the inner product of the same error vector <ei|ei> 
C     (if we have beta matrix, calculate the beta inner product and add
C     it to alpha
      Fe = SCFTrX(.False.,.False.,IOpCl,1,NBsUse,AElec(IA),AElec(IB),
     $  AElec(IA),AElec(IB))
      Fp = SCFTrC(AProt,AProt,NBuNEO,1)

C     Function SProd(N, A, B):
C     Returns the scalar product of vectors A and B.
      Qe = SProd(NVQ,AElec(IQ),AElec(IQ))
      
C     Simultaneous Opt: add electronic error vector and protonic error
C     vector together      
      BMat(NSaved+1,NSaved+1) = WFe*Fe + WQe*Qe + Fp
      BMat(1,NSaved+1) = -One
      BMat(NSaved+1,1) = -One
      Do 10 I = 1, (NSaved-1)
        Call FileIO(2,-IOErrE,LenErE,BElec,IndFok(I)*LenErE)
        Fe = SCFTrX(.False.,.False.,IOpCl,1,NBsUse,AElec(IA),AElec(IB),
     $    BElec(IA),BElec(IB))

        Call FileIO(2,-IOErrP,LenErP,BProt,IndFok(I)*LenErP)
        Fp = SCFTrC(AProt,BProt,NBuNEO,1)
        
        Qe = SProd(NVQ,AElec(IQ),BElec(IQ))
        BMat(NSaved+1,I+1) = WFe*Fe + WQe*Qe + Fp
   10   BMat(I+1,NSaved+1) = BMat(NSaved+1,I+1)
      Return
      End
