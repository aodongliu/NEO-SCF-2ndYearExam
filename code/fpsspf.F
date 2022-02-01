*Deck FPSSPF
      Subroutine FPSSPF(NRI,NBasis,NBsUse,IRwP,IRwF,IRwFX,IRwS,IRwV,
     $  IRwErr,NumFS,A,C,D,Scr,LenScr)
      Implicit Real*8(A-H,O-Z)
C
C     Form the commutator of F and P (FPS-SPF).  This is returned
C     as a lower triangle (including diagonal) in A.  A is 2*NTT,
C     C is NBasis**2 and D is NBasis.  IRwV is used to
C     transform the commutator to an orthonormal basis.  If IRwS
C     is 0, the basis is assumed to be orthonormal and IRwS and IRwV
C     are not used. If IRwErr is non-zero, the lower triangular error
C     matrix is taken from it.
C
      Logical HavScr
      Dimension A(*), C(NBasis,NBasis), D(NBasis), Scr(LenScr)
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
C
      NTT = (NBasis*(NBasis+1))/2
      NTTR = NRI*NTT
      NTTUR = NRI*(NBsUse*(NBsUse+1))/2
      NBSq = NBasis**2
      HavScr = LenScr.ge.(NRI*NBSq)
      LenF = NRI*NumFS*NTT
C     IRwF  ... RWF on which current alpha Fock matrix is stored.
C     IRwFX ... Fock matrix to extrapolate.  If 0, the first matrix in IRwFA
C                is used.
C     Save the Fock matrix in A
      Write(*,*) "AL 6"
      Write(*,*) IRwV
      Call FileIO(2,-IRwV,NBasis*NBsUse*NRI,C,0)

      Write(*,*) IRwF
      Call FileIO(2,-IRwF,LenF,A,0)
      Write(*,*) "AL 7"
      Write(*,*) "AL 7.1"
C     IRwErr ... Optional rwf of error matrix in AO basis.  If this is
C                non-zero, this is used instead of [F,P].  The matrix
C                is stored lower triangular, with beta following alpha
C                if uhf.
C     If IRwErr = 0, we form the commutator of F and P
      If(IRwErr.eq.0) then

        Write(*,*) "AL 7.5"
        Write(*,*) IRwV
C       If IRwFX != 0, then save the alternate Fock matrix in A
        If(IRwFX.ne.0) Call FileIO(2,-IRwFX,NTTR,A,0)

        Write(*,*) "AL 8"
        Write(*,*) IRwV
C       Save the density matrix in A(NTTR+1)
        
        Call FileIO(2,-IRwP,NTTR,A(NTTR+1),0)
 
        Write(*,*) "AL 9"
        Write(*,*) IRwV
C       If we have enough scratch space, we form the square matrices and
C       do matrix operation
        If(HavScr) then
          Write(*,*) "AL 9.1"
          Write(*,*) IRwV
C         XSquare:
C         Places linear array in square form. 0 means hermitian
          Call XSquare(NRI,A(NTTR+1),Scr,NBasis,NBasis,0)
          Call XSquare(NRI,A,A,NBasis,NBasis,0)

C         Subroutine XGEMM(NRI,Str1,Str2,M,N,K,Alpha,A,LDA,B,LDB,BetaI,C,
C         $  LDC):
C         Wrapper for matrix multiplies which handles both real/complex and
C         OpenMP parallelism.
C         C := alpha*op( A  )*op( B  ) + beta*C,
C         where  op( X  ) is one of
C                op( X  ) = X   or   op( X  ) = X**T,
          Call XGEMM(NRI,'N','N',NBasis,NBasis,NBasis,One,A,NBasis,Scr,
     $      NBasis,Zero,C,NBasis)
C         A contains square matrix of F
C         Scr contains square matrix of P
C         This XGEMM call will calculate FP and stores it in C          
          Write(*,*) "AL 9.2"
          Write(*,*) IRwV
        else
C         IWEHT we don't have enough scratch space to do square matrix
C         operations

          Write(*,*) "AL 9.3"
          Write(*,*) IRwV
          If(NRI.eq.2) Call GauErr('NRI=2 and not HavScr in FPSSPF.')

C         MltLL1:                  
C         Multiply lower triangular matrices A and B and store into C.
C         C is a square matrix          
          Call MltLL1(NBasis,A,A(NTT+1),C)
          endIf
C       C now contains the square matrix FP

        Write(*,*) "AL 9.4"
        Write(*,*) IRwV

C       IRwS   ... RWF of overlap matrix.  0 if calculation is in an
C                  orthonormal basis
C       If IRwS != 0, then we need to multiply S to FP 
        If(IRwS.ne.0) then

          Write(*,*) "AL 9.5"
          Write(*,*) IRwV
C         Save S in A and convert it into a square matrix         
          Write(*,*) "AL 10"
          Call FileIO(2,-IRwS,NTTR,A,0)
          Write(*,*) "AL 11"
          Call XSquare(NRI,A,A,NBasis,NBasis,0)

C         If we have enough scratch space, we directly call XGEMM to get
C         FPS
          If(HavScr) then
            Call XGEMM(NRI,'N','N',NBasis,NBasis,NBasis,One,C,NBasis,A,
     $        NBasis,Zero,Scr,NBasis)
C           A contains square matrix of S
C           C contains square matrix of FP
C           This XGEMM call will calculate FPS and stores it in Scr          

C           ASbATC:
C           Form the lower triangular antiHermitean matrix A = D - D(dagger)

C           Scr contains FPS, therefore this ASbATC call will calculate
C           FPS-SPF
            Call ASbATC(NRI,NBasis,Scr,A)
C           A now stores the lower triangular matrix of [F,P]            
          else
C           IWEHT we don't have enough scratch space, call MatRec
C           routine

            If(NRI.eq.2) Call GauErr('NRI=2 and not HavScr in FPSSPF.')

C           Subroutine MatRec(A,B,D,MaxDim,L,M,N,Mode)
C           Matrix multiplication.  All matrices must be allocated square.
C           Mode = 1 ... A = A*B
C           Mode = 2 ... A = A(T)*B
C           Mode = 3 ... A = A*B(T)
C           Mode = 4 ... B = A*B
C           Mode = 5 ... B = A(T)*B
C           Mode = 6 ... B = A*B(T)

C           A contains square matrix of S
C           C contains square matrix of FP
C           This MatRec call will calculate FPS and stores it in C          
            Call MatRec(C,A,D,NBasis,NBasis,NBasis,NBasis,1)

C           C now contains FPS
C           Subroutine ASubAT(NBasis,D,A)
C           Form the lower triangular antisymmetric matrix A = D - Dt.
            Call ASubAT(NBasis,C,A)
C           A now stores the lower triangular matrix of [F,P]            
            endIf

        else if(NRI.eq.2) then
C         IWEHT F and P are already in the orthrnomal basis. No need to
C         consider S

C         NRI = 2, we want to use ASbATC (square matrix - square
C         matrix dagger) to get [F,P] into a lower triangular matix A
          Call ASbATC(NRI,NBasis,C,A)
        else
C         Subroutine ASubAT(NBasis,D,A)
C         Form the lower triangular antisymmetric matrix A = D - Dt.

C         NRI = 1, we want to use ASubAT (square matrix - square
C         matrix transpose) to get [F,P] into a lower triangular matix A
          Call ASubAT(NBasis,C,A)
          endIf
C       A now stores the lower triangular matrix of [F,P]             


        Write(*,*) "AL 11.5"
        Write(*,*) IRwV

C       Make sure [F,P] is in the orthonormal basis
        If(IRwS.ne.0) then

C         Save V (the transformation matrix) in C          
          Write(*,*) "AL 12"
          Write(*,*) IRwV
          Call FileIO(2,-IRwV,NBasis*NBsUse*NRI,C,0)
          Write(*,*) "AL 13"
          If(NRI.eq.2) then
C           Subroutine ZVAVMl(IType,ASymm,Lin,NIn,NOut,FI,V,A)
C           Routine to form matrix product:
C           IType:  0 ...  FI = V(Dagger)*FI*V
C           IType:  1 ...  FI = V*FI*V(Dagger)
C           This routine is somewhat similar to vfvmul for complex matrices, but requires that
C           FI be allocated N**2, even though the input and output may be linear.
C    
C           NOut may be < than NIn
            Call ZVAVMl(0,.True.,.True.,NBasis,NBsUse,A,C,Scr)
          else
            Ind = NBSq - NTT + 1
C           Subroutine AMove1(N,IndF,IndT,A)
C           Move N elements from A(IndF+1) to A(IndT+1), doing the right
C           thing regardless of whether IndF is larger or smaller than IndF.
            Call AMove1(NRI*NTT,0,NRI*(Ind-1),A)
C           Subroutine VFVMul(ASymm,IAIn,IAOut,NO,NIn,NOut,NDim,FI,V,T,FO,S,LenS)
C           Routine to form the matrix product:
C           FO = V(Dagger)*FI*V
            Call VFVMul(.True.,.False.,.False.,0,NBasis,NBsUse,NBasis,
     $        A(Ind),C,A,A,Scr,LenScr)
            endIf
          endIf

      else
        Call FileIO(2,IRwErr,NTTUR,A,0)
        endIf
      Return
      End
