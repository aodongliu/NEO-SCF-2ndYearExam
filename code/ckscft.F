*Deck CkSCFT
      Subroutine CkSCFT(IOut,Name,IOpCl,CompOK,CBfnOK,SpinOK,NEOOK,
     $  NSpBlk,NRI,NDimBl,CBfn,NSpBlX,NSpBNX,Spinor)
      Implicit Integer(A-Z)
C
C     Decode IOpCl, check whether the type of calculation is supported,
C     and return separate NSpBlX (electronic) and NSpBNX (nuclear).
C
      Character*(*) Name
      Logical CompOK,CBfnOK,SpinOK,NEOOK,CBfn,Spinor
 1000 Format(1X,A,' not supported in ',A,'.')
C
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlY,Spinor)
      NSpBlX = Mod(NSpBlY,10)
      NSpBNX = Mod(NSpBlY,100)/10
      If(NRI.eq.2.and..not.CompOK) then
        Write(IOut,1000) 'Complex orbitals', Name
        Call Lnk1E(0)
        endIf
      If(CBfn.and..not.CBfnOK) then
        Write(IOut,1000) 'Complex basis functions', Name
        Call Lnk1E(0)
        endIf
      If(Spinor.and..not.SpinOK) then
        Write(IOut,1000) 'Spinor basis functions', Name
        Call Lnk1E(0)
        endIf
      If(NSpBNX.gt.0.and..not.NEOOK) then
Cxslis
C        Write(IOut,1000) 'Nuclear orbitals not supported in', Name
C        Call Lnk1E(0)
Cxslie
        endIf
      Return
      End
