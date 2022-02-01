*Deck FoFNEO
      Subroutine FoFNEO(IOut,IPrint,IUHF1,UsePT,DoPurP,DoPurF,NBasis,
     $  NBas6D,NBsNEO,NBsNE6,Omega,FMM,FMFlag,FMFlg1,NFxFlg,IPFlag,
     $  AllowP,ISym2E,NAtoms,IAn,AtmChg,C,IAtTyp,IAtFrg,NOpSet,NOpAb,
     $  NEqAtm,NEqShl,RotOp,NEqBs6,NEqNES,NENEB6,AccDes,PT,PN,FA,FB,FN,
     $  ESum,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Increment FA and FB with the Coulomb interaction with the NEO
C     atoms, and FN with the Coulomb and exhange terms.  If UsePT is
C     true, then PN is also incremented with the Coulomb interaction
C     with the electronic density PT.
C
#include "commonb2.inc"
      Parameter (IHFJ=6)
      Logical UsePT, DoPurP, DoPurF, FMM, AllowP(*), LJJ(1)
      Integer FMFlag, FMFlg1
      Dimension Omega(*), NEqAtm(*), NEqShl(*), RotOp(*), NEqBs6(*),
     $  NEqNES(*), NENEB6(*), PT(*), PN(*), FA(*), FB(*), FN(*), V(*),
     $  XX(1), JJ(1), IAn(*), AtmChg(*), C(*), IAtTyp(*), IAtFrg(*),
     $  ESum(10)
      Save XX, LJJ, One, JJ, Zero, Pt5, Two
      Data XX(1)/0.0d0/, LJJ/.False./, One/1.0d0/, JJ/0/, Zero/0.0d0/,
     $  Pt5/0.5d0/, Two/2.0d0/
C
      NShelA = NShelB
      Call RdWrB(20,0)
      NShelN = NShell
      NShelT = NShelA + NShelN
      NRI = 1
      NTT = (NBasis*(NBasis+1))/2
      NEONTT = (NBsNEO*(NBsNEO+1))/2
      NBTot6 = NBas6D + NBsNE6
      NBTT6 = (NBTot6*(NBTot6+1))/2
      NBTSq6 = NBTot6**2
      INETSh = 1
      INETB6 = INETSh + NOpSet*NShelT
      IPM = INETB6 + NOpAb*NBTot6
      IFJ = IPM + NBTSq6
      IFK = IFJ + NBTSq6
      IV = IFK + NBTSq6
      Call TstCor(IV,MDV,'FoFNEO')
      MDV1 = MDV - IV + 1
      IPurND = I5DB1
      IPurNF = I7FB1
      Call LTOutN(IOut,'NEO Density:',0,1,NBsNEO,PN,0)
      If(DoPurP.and.NBsNEO.ne.NBsNE6) then
        Call AMove(NEONTT,PN,V(IPM))
        Call P2Car1(IOut,IPrint,IPurND,IPurNF,.False.,.True.,.False.,
     $    .False.,.False.,.False.,.False.,.False.,.True.,NBsNEO,NBsNE6,
     $    0,1,1,0,V(IPM),XX,XX,XX,XX,XX,JJ,JJ,XX,LJJ,LJJ,V(IV),MDV1)
        Call Square(V(IPM),V(IFJ),NBsNE6,NBsNE6,3)
      else
        I5DB1 = 1
        I7FB1 = 1
        Call GenAOS(I5DB1,I75B1,NShell,ShellT,ShellC,Junk,AOS)
        Call Square(PN,V(IFJ),NBsNE6,NBsNE6,3)
        endIf
      Call RdWrB(5,0)
      Call RdWrB(3,0)
      IPurAD = I5DB1
      IPurAF = I7FB1
      If(DoPurP.and.NBasis.ne.NBas6D) then
        Call P2Car1(IOut,IPrint,IPurAD,IPurAF,.False.,.True.,.False.,
     $    .False.,.False.,.False.,.False.,.False.,.True.,NBasis,NBas6D,
     $    0,1,1,0,PT,XX,XX,XX,XX,XX,JJ,JJ,XX,LJJ,LJJ,V(IV),MDV1)
      else
        I5DB1 = 1
        I7FB1 = 1
        Call GenAOS(I5DB1,I75B1,NShell,ShellT,ShellC,Junk,AOS)
        endIf
      Call Appbas(1)
      Call AClear(NBTSq6,V(IPM))
C
C     Put the nuclear density in the lower diagonal block of the total
C     (electron + nuclear) density
      Call MatSub(NBTot6,NBTot6,NBas6D+1,NBTot6,NBas6D+1,NBTot6,V(IFJ),
     $  V(IPM))
      Call Linear(V(IPM),V(IPM),NBTot6,NBTot6)
      Call LTOutN(IOut,'Total Density:',0,1,NBTot6,V(IPM),0)
      If(ISym2E.ne.0) then
        Call MrgNEB(NOpAb,NBas6D,NBsNE6,NBTot6,NEqBs6,NENEB6,V(INETB6))
        Call MrgNEB(NOpSet,NShelA,NShelN,NShelT,NEqShl,NEqNES,V(INETSh))
        endIf
      IRaf = 0
      I1Cent = 0
      NBT = 0
      Energy = Zero

C     Here, we scale the protonic density by one half
      Call AScale(NBTT6,Pt5,V(IPM),V(IPM))
      IDoPuP = 2000020

      Write(*,*) 'BEFORE CALLING FOFCOU'      
      Call FoFCou(IOut,IPrint,1,1,500,IRaf,IHFJ,Omega,FMM,FMFlag,FMFlg1,
     $  NFxFlg,IPFlag,AllowP,0,XX,1,NAtoms,.True.,IDoPuP,0,.False.,
     $  I1Cent,AccDes,One,1,1,1,0,1,0,0,NBTot6,ISym2E,NOpSet,NOpAb,
     $  NEqAtm,V(INETSh),JJ,RotOp,V(INETB6),NBT,JJ,V(IPM),XX,XX,XX,XX,
     $  XX,XX,XX,XX,JJ,V(IFJ),V(IFK),JJ,XX,XX,0,0,XX,XX,XX,XX,XX,XX,XX,
     $  XX,XX,XX,XX,NAtoms,IAn,AtmChg,C,IAtTyp,IAtFrg,XX,XX,JJ,XX,
     $  Energy,One,XX,0,0,XX,XX,JJ,JJ,XX,V(IV),V(IV),MDV1)
C     Form NEO 2-JK

C     For the big J matrix returned from fofcou, we have -FeN in the
C     upper diagonal block, and we have JNN in the lower diagonal
C     block. 


C     For the big K matrix returned from fofcou, we have a mixed K term
C     in the upper diagonal block (not important), and we have 0.5 KNN
C     in the lower diagonal block. 

      Call LTOutN(IOut,'Big J Matrix from FoFNEO:',0,1,NBTot6,V(IFJ),0)
      Call LTOutN(IOut,'Big K Matrix from FoFNEO:',0,1,NBTot6,V(IFK),0)

C     The pure protonic fock matrix FNN should be JNN-KNN (+Core HNN)
C     In order to piece these together, we need to do -2 * (0.5 * KNN)
C     + 1 * JNN, which is done here, and the result is saved in the IFK
      Call ASASB(NBTT6,-Two,V(IFK),One,V(IFJ))
C
C     Extract the upper diagonal block of the the big J matrix (stored
C     in IFJ)
C     When the nuclear density is used in FoFCou contraction, the upper 
C     diagonal block contains the FeN contribution to Fe
      Call Square(V(IFJ),V(IPM),NBTot6,NBTot6,3)
      Call SubMat(.False.,NRI,NBTot6,NBTot6,1,NBas6D,1,NBas6D,V(IPM),
     $  V(IFJ))
      Call Linear(V(IFJ),V(IFJ),NBas6D,NBas6D)
C
C     Extract the lower diagonal block the big FNN matrix that we just
C     formed in line 119 (stored in IFK)
C     When the nuclear density is used in FoFCou contraction, the lower
C     diagonal block contains the FNN without the FNe contribution
      Call Square(V(IFK),V(IPM),NBTot6,NBTot6,3)
      Call SubMat(.False.,NRI,NBTot6,NBTot6,NBas6D+1,NBTot6,NBas6D+1,
     $  NBTot6,V(IPM),V(IFK))
      Call Linear(V(IFK),V(IFk),NBsNE6,NBsNE6)
      Call RdWrB(-5,0)
      If(DoPurF.and.(IPurND.ne.1.or.IPurNF.ne.1)) Call F2Pure(IOut,
     $  IPrint,IPurND,IPurNF,1,1,0,NBsNEO,NBsNE6,.True.,.True.,.False.,
     $  .False.,.False.,.False.,V(IFK),XX,XX,XX,V(IV),MDV1)
C

C     We haven't touched FN yet. It's still holding the the protonic
C     core hamiltonian. This allows us to calculate:
C     Classical-quantum nuclear repulsion
      ESum(3) = SCFTrc(PN,FN,NBsNEO,1)
C
C     IFK now is holding the FNN term (JNN - KNN). We can calculate:
C     Quantum nuclear repulsion energy. 
      ESum(4) = Pt5*SCFTrc(PN,V(IFK),NBsNEO,1)
      Call AAdd(NEONTT,V(IFK),FN,FN)
      Call LTOutN(IOut,'FNN from FoFNEO:',0,1,NBsNEO,V(IFK),0)

      Call RdWrB(3,0)
      If(DoPurF.and.(IPurAD.ne.1.or.IPurAF.ne.1)) Call F2Pure(IOut,
     $  IPrint,IPurAD,IPurAF,1,1,0,NBasis,NBas6D,.True.,.True.,.False.,
     $  .False.,.False.,.False.,V(IFJ),XX,XX,XX,V(IV),MDV1)
      Call LTOutN(IOut,'FeN from FoFNEO:',0,1,NBasis,V(IFJ),0)
C      Call FileIO(2,-515,NTT,V(IV),0)
C      Call LTOutN(IOut,'Core Hamiltonian:',0,1,NBasis,V(IV),0)
C      Call ASub(NTT,V(IV),V(IFJ),V(IV))
C      Call LTOutN(IOut,'Core Hamiltonian + NEO:',0,1,NBasis,V(IV),0)
C
C     
C     IFJ is holding the FeN term. We can calculate: 
C     FeN energy. 
      ESum(2) = -Two*SCFTrc(PT,V(IFJ),NBasis,1)
      Call ASub(NTT,FA,V(IFJ),FA)
      If(IUHF1.eq.1) Call ASub(NTT,FB,V(IFJ),FB)
C      Call LTOutN(IOut,'Fe from FoFNEO:',0,1,NBasis,FA,0)
C
C
C     NEXT: we compute JNe contribution to FN
C
      Call LTOutN(IOut,'Electron Density:',0,1,NBas6D,PT,0)
C
C     Put the electron density in the upper diagonal block of the total
C     (electron + nuclear) density
      Call RdWrB(20,0)
      Call RdWrB(5,0)
      Call RdWrB(3,0)
      Call Appbas(1)
      Call AClear(NBTSq6,V(IPM))
      Call Square(PT,V(IFJ),NBas6D,NBas6D,3)
      Call MatSub(NBTot6,NBTot6,1,NBas6D,1,NBas6D,V(IFJ),V(IPM))
      Call Linear(V(IPM),V(IPM),NBTot6,NBTot6)
      Call LTOutN(IOut,'Total Density:',0,1,NBTot6,V(IPM),0)
      IRaf = 0
      I1Cent = 0
      NBT = 0
      Energy = Zero

C     This next line is commented out because there's not need to scale
C     the electronic density by one half
C      Call AScale(NBTT6,Pt5,V(IPM),V(IPM))
      IDoPuP = 2000020
      Call FoFCou(IOut,IPrint,1,1,500,IRaf,IHFJ,Omega,FMM,FMFlag,FMFlg1,
     $  NFxFlg,IPFlag,AllowP,0,XX,1,NAtoms,.True.,IDoPuP,0,.False.,
     $  I1Cent,AccDes,One,1,1,1,0,1,0,0,NBTot6,ISym2E,NOpSet,NOpAb,
     $  NEqAtm,V(INETSh),JJ,RotOp,V(INETB6),NBT,JJ,V(IPM),XX,XX,XX,XX,
     $  XX,XX,XX,XX,JJ,V(IFJ),V(IFK),JJ,XX,XX,0,0,XX,XX,XX,XX,XX,XX,XX,
     $  XX,XX,XX,XX,NAtoms,IAn,AtmChg,C,IAtTyp,IAtFrg,XX,XX,JJ,XX,
     $  Energy,One,XX,0,0,XX,XX,JJ,JJ,XX,V(IV),V(IV),MDV1)
C
C     
C     After this second fofcou call, we will be getting the following
C     In the big J matrix, the upper diagonal block will be holding 2 *
C     Jee, the lower diagonal block will be holding -FNe.

C     In the big K matrix, the upper diagonal block will be holding Kee,
C     and the lower diagonal block will be holding a mixed K term (not
C     important for purpose)

C     Technically we can piece the upper diagonal blocks together and
C     recover the pure electronic fock matrix (minus core H) Fee, but
C     this terms is alreay computed in MkFok1. So we do need to anything
C     with the upper diagonal blocks

C     DEBUG-----------------------------------------      
C     For debugging purpose, I now assemble the Fee matrix and see if it
C     matches the expected value

C     First I do 2Jee - Kee (which will give us Fee) and store it in IFK
      Call ASASB(NBTT6,-One,V(IFK),One,V(IFJ))

C     Then I extract the upper diagonal block
      Call Square(V(IFK),V(IPM),NBTot6,NBTot6,3)
      Call SubMat(.False.,NRI,NBTot6,NBTot6,1,NBas6D,1,NBas6D,V(IPM),
     $  V(IFK))
      Call Linear(V(IFK),V(IFK),NBas6D,NBas6D)
      If(DoPurF.and.(IPurAD.ne.1.or.IPurAF.ne.1)) Call F2Pure(IOut,
     $  IPrint,IPurAD,IPurAF,1,1,0,NBasis,NBas6D,.True.,.True.,.False.,
     $  .False.,.False.,.False.,V(IFK),XX,XX,XX,V(IV),MDV1)

      Call LTOutN(IOut,'DebugFee:',0,1,NBasis,V(IFK),0)

C     DEBUG-----------------------------------------      

C      Call LTOutN(IOut,'J from FoFNEO:',0,1,NBTot6,V(IFJ),0)
C      Call LTOutN(IOut,'K from FoFNEO:',0,1,NBTot6,V(IFK),0)
C
C     Extract the lower diagonal block of the big J matrix (which is
C     stored in IFJ), and store it in IFJ
C     When the electron density is used in FoFCou contraction, the lower
C     diagonal block contains the FNe contribution to FN

C     Note, there's no need to scale anything here!      
      Call Square(V(IFJ),V(IPM),NBTot6,NBTot6,3)
      Call SubMat(.False.,NRI,NBTot6,NBTot6,NBas6D+1,NBTot6,NBas6D+1,
     $  NBTot6,V(IPM),V(IFJ))
      Call Linear(V(IFJ),V(IFJ),NBsNE6,NBsNE6)
      Call LTOutN(IOut,'FNe from FoFNEO:',0,1,NBsNE6,V(IFJ),0)
      Call RdWrB(-5,0)
      If(DoPurF.and.(IPurND.ne.1.or.IPurNF.ne.1)) Call F2Pure(IOut,
     $  IPrint,IPurND,IPurNF,1,1,0,NBsNEO,NBsNE6,.True.,.True.,.False.,
     $  .False.,.False.,.False.,V(IFJ),XX,XX,XX,V(IV),MDV1)
      Call LTOutN(IOut,'FNe from FoFNEO:',0,1,NBsNEO,V(IFJ),0)

C     IFJ is holding FNe, we will be able to calculate:
C     Quantum-nuclear-electronic potential: 
      ESum(5) = -SCFTrc(PN,V(IFJ),NBsNEO,1)
      Call ASub(NEONTT,FN,V(IFJ),FN)
C      Call LTOutN(IOut,'FN from FoFNEO:',0,1,NBsNEO,FN,0)
      Call RdWrB(3,0)
C
      Return
      End
