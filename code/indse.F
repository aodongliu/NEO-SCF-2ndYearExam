*Deck IndSE
      Subroutine IndSE(IHMeth,IStart,ISEScl,IHMPar,MaxASE,NSEScI,NSEScR,
     $  NSEAtI,NSEAtR,NCorPA,MaxCPA,NSEElI,NSEElR,MxElBE,jRSESc,jISEP,
     $  jRSEP,jSECP,jIElP,jRElP,jBELP,jBEBO,Next)
      Implicit Integer(A-Z)
C
C     Generate dimensions and pointers into the SE parameter block.
C
#include "separ.inc"
      Dimension ISEScl(*)
C
      ISETyp = ISEScl(ISEPTy)
Cxslis
      ISETyp = 0
Cxslie
      If(ISETyp.eq.0) then
        IHMPar = IHMeth
        MaxASE = 0
        NSEScI = 0
        NSEScR = 0
        NSEAtI = 1
        NSEAtR = 1
        NCorPA = 1
        NSEElI = 1
        NSEElR = 1
        MaxCPA = 1
        MxElBE = 0
        jRSESc = IStart
        jISEP  = jRSESc
        jRSEP  = jISep
        jSECP  = jRSep
        jIElP  = jSECP
        jRElP  = jIElP
        jBELP  = jRElP
        jBEBO  = jBELP
        Next   = IStart
      else
        If(ISETyp.lt.1000.or.ISETyp.gt.CurSEV)
     $    Call GauErr('Bad data passed to IndSE.')
        IHMPar = ISEScl(ISEMod)
        MaxASE = ISEScl(IMaxAS)
        NSEScI = ISEScl(INSScI)
        NSEScR = ISEScl(INSScR)
        NSEAtI = ISEScl(INSESI)
        NSEAtR = ISEScl(INSESR)
        NCorPA = ISEScl(INCrPA)
        NSEElI = ISEScl(INSElI)
        NSEElR = ISEScl(INSElR)
        MaxCPA = ISEScl(IMaxCP)
        MxElBE = ISEScl(IMxABE)
        jRSESc = IStart + InToWP(ISEScl(INSScI))
        jISEP  = jRSESc + ISEScl(INSScR)
        jRSEP  = jISep + InToWP(NSEAtI*MaxASE)
        jSECP  = jRSep + NSEAtR*MaxASE
        jIElP  = jSECP + NCorPA*MaxCPA*MaxASE
        MASETT = (MaxASE*(MaxASE+1))/2
        jRElP  = jIElP + InToWP(NSEElI*MASETT)
        jBELP  = jRElP + NSEElR*MASETT
        If(MxElBE.gt.0) then
          MxElB1 = MxElBE + 1
        else
          MxElB1 = 0
          endIf
        jBEBO  = jBELP + 2*MxElB1
        Next   = jBEBO + MxElB1*(MxElB1+1)
        endIf
      Return
      End
