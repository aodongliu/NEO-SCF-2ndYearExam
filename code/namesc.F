*Deck NameSC
      Subroutine NameSC(IType,IOpClI,IROHF,IHMeth,IExCor,IXCFnc,IEmpDi,
     $  N6Func,ScaHFX,SCFNam)
      Implicit Integer(A-Z)
C
C     Return the short or long name for a particular SCF procedure.
C     IType = 1 ... Short abbreviation for Open/Closed etc.
C             2 ... Long name for Open/Closed etc.
C             3 ... Full name including method.
C
      Real*8 ScaHFX
      Parameter (MaxIOC=22,MxXCNm=DEFMAXXCNAME,MaxROT=3)
Cxslis
      Parameter (MskNEO=1024)
      Logical NucOrb
Cxslie
      Character SCFNam*(*), NameS(0:MaxIOC)*2, NameL(0:MaxIOC)*36,
     $  OpClNm*2, ExName*(MxXCNm), CoName*(MxXCNm), NameRS(1:MaxROT)*5,
     $  NameRL(1:MaxROT)*36, N6Name*(MxXCNm)
      Save NameS, NameL, NameRS, NameRL
      Data NameS/'R','U','CR','CU',2*' ','G',3*' ','YR','YU',2*' ','YG',
     $  3*' ','SG',3*' ','D'/,NameRS/'RO','GVB','MCSCF'/,NameRL/
     $  'Restricted open shell','Generalized Valence Bond','MCSCF'/
      Data NameL/'Closed shell','UHF open shell','Complex closed shell',
     $  'Complex UHF open shell',2*' ','General Hartree-Fock',3*' ',
     $  'Complex RHF (Ylm basis)','Complex UHF (Ylm basis)',2*' ',
     $  'General Hartree-Fock (Ylm basis)',3*' ',
     $  'General Hartree-Fock (Spinor basis)',3*' ','Dirac-Fock'/
 1000 Format('SCF ROType',I1)
 1010 Format('RO',I1)
C
Cxslis
      IOpCl = Mod(Max(IOpClI,0),MskNEO)
      NucOrb = IAnd(IAbs(IOpClI),MskNEO).ne.0
Cxslie
      If(IOpCl.lt.0.or.IOpCl.gt.MaxIOC)
     $  Call GauErr('IOpCl out of range in NameSC.')
      If(IType.eq.2) then
        If(IROHF.gt.MaxROT) then
          SCFNam = ' '
          Write(SCFNam,1000) IROHF
        else if(IROHF.gt.0) then
          SCFNam = NameRL(IROHF)
        else
          SCFNam = NameL(IOpCl)
          endIf
      else
        If(IROHF.gt.MaxROT) then
          OpClNm = ' '
          Write(OpClNm,1010) IROHF
        else if(IROHF.gt.0) then
          OpClNm = NameRS(IROHF)
        else
          OpClNm = NameS(IOpCl)
          If(OpClNm.eq.' ') Call GauErr('Unrecognized IOpCl in NameSC.')
          endIf
        If(IType.eq.1.or.IROHF.gt.1) then
          SCFNam = OpClNm
        else
          Call DFTNam(IHMeth,IExCor,IXCFnc,IEmpDi,N6Func,ScaHFX,ExName,
     $      CoName,N6Name,MetExc,MetCor)
          SCFNam = OpClNm(1:LinEnd(OpClNm))//ExName(1:LinEnd(ExName))
          If(CoName(1:5).ne.'-None')
     $      SCFNam(LinEnd(SCFNam)+1:) = CoName(1:LinEnd(CoName))
          If(N6Name(1:5).ne.'-None')
     $      SCFNam(LinEnd(SCFNam)+1:) = N6Name(1:LinEnd(N6Name))
          endIf
        endIf
      Return
      End
