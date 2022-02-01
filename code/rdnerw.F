*Deck RdNERw
      Subroutine RdNERw(IOut,IOpCl,IFileI,IFileO,NQNuc,NBsNEO,NBsNE6,
     $  NBUNEO,IOfNES,IOfNET,IOfNEH,IOfNER,IOpClE,IOClNP,IOClND)
      Implicit Real*8(A-H,O-Z)
C
C     If NSpBNX>0 read NEO dimension information and return the
C     number of the file holding the one-electron matrix and
C     the offsets to each matrix.  If IFileI is 0, the default
C     file is used.
C
C     IOpCl  ... If positive, used to determine if NEOs are in use,
C                returning zeros if they are not and terminating if
C                NEOs are in use but the file is not found.  -1 to
C                recover values if the file exists and return zeros if
C                it does not.
C     IFileI ... File to read.  If 0, the default is used
C     IFileO ... Returned file to use.
C     NQNuc  ... Number of quantum nuclei.
C     NBsNEO ... Number of NEO basis functions.
C     NBsNE6 ... Number of Cartesian NEO basis functions.
C     NBUNEO ... Number of linearly independent NEO functions.
C     IOfNES ... Offset in file IFileO for the NEO overlap matrix
C                (lower triangular).
C     IOfNET ... Offset in file IFileO for the NEO kinetic energy
C                matrix (lower triangular).
C     IOfNEH ... Offset in file IFileO for the NEO Core Hamiltonian
C                matrix (KE + repulsion from classical nuclei)
C                (lower triangular).
C     IOfNER ... Offset to transformation matrix giving the linear
C                independent NEO functions, (NBsNEO,NBUNEO).
C     IOpClE ... IOpCl for electronic wavefunction.
C     IOCLNP ... IOpCl-like value for printing the NEO orbitals or
C                density.
C     IOCLND ... IOpCl-like value for forming the NEO density.
C
      Parameter (IRwNEX=834,LNEOSc=16,INEONB=1,INEON6=2,INEONU=3,
     $  INEONP=4)
      Logical CBfn, Spinor, HavNEX
      Dimension NEOScl(LNEOSc), RNEOSc(LNEOSc)
      Equivalence (NEOScl(1),RNEOSc(1))
C
      Write(*,*) 'AL a'
      Call IClear(LNEOSc,NEOScl)
      Write(*,*) 'AL b'
      LScl = 0
      Call IndNEO(1,IFileI,IFileO,NEOScl,IOfNES,IOfNET,IOfNEH,IOfNER,
     $  LenNEX)
      Write(*,*) 'AL c'
      HavNEX = ITqry(IFileO).gt.0
      If(IOpCl.ge.0) then
        Write(*,*) 'AL d'
        Call SpIOCl(IOut,'RdNERw',IOpCl,.True.,.True.,.True.,.True.,
     $    NSpBlk,NRI,NDimBl,CBfn,NSpBlX,NSpBNX,Spinor,IOpClE,IOClNP,
     $    IOClND)
        If(NSpBNX.gt.0.and..not.HavNEX)
     $    Call GauErr('NSpBNX>0 but no NEO 1e matrices in RdNERw.')
        endIf
      If(HavNEX) then
        Write(*,*) 'AL e'
        LScl = InToWP(LNEOSc)
        Call FileIO(2,-IFileO,LScl,RNEOSc,0)
        Call IndNEO(1,IFileI,IFileO,NEOScl,IOfNES,IOfNET,IOfNEH,IOfNER,
     $    LenNEX)
        Write(*,*) 'AL f'
      else
        IFileO = 0
        endIf
      NQNuc = NEOScl(INEONP)
      NBsNEO = NEOScl(INEONB)
      NBsNE6 = NEOScl(INEON6)
      NBUNEO = NEOScl(INEONU)
      Return
      End
