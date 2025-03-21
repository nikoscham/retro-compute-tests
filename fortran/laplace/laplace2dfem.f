      PROGRAM LAPLACE
C     Solves Laplace's equation using FEM with quadratic basis functions
C     Designed for Windows 9x/Me with Fortran 77 compilers (e.g., PowerStation 4.0)
      IMPLICIT NONE
C     Parameter to set problem size (number of elements in x-direction)
C     NEX=10 gives 441 nodes (21x21), NEX=5 gives 121 nodes (11x11)
      INTEGER NEX
      PARAMETER (NEX=10)
C     Derived grid sizes: NNX = 2*NEX+1, NNY = 2*NEY+1, NP = NNX*NNY
      INTEGER NNX, NNY, NP, NEY
      PARAMETER (NEY=NEX)          ! Equal elements in y-direction
      PARAMETER (NNX=2*NEX+1)      ! Nodes in x-direction
      PARAMETER (NNY=2*NEY+1)      ! Nodes in y-direction
      PARAMETER (NP=NNX*NNY)       ! Total nodes
C     Arrays: sized for NEX=10 (441 nodes max), adjust if NEX > 10
      REAL SK(441,441)             ! Stiffness matrix
      REAL X3(441)                 ! Solution vector
      REAL R3(441)                 ! Right-hand side vector
      REAL AXPT(441), AYPT(441)    ! Node coordinates
      INTEGER INDX(441)            ! Pivot indices for solver
      INTEGER NOP(100,9)           ! Element-to-node connectivity
C     Local variables
      INTEGER I
C     Open files for output
      OPEN(1, FILE='results.txt')    ! Final results (x, y, solution)
      OPEN(10, FILE='data_in.txt')   ! Input data (unused here)
      OPEN(20, FILE='data_out.txt')  ! Solution values only
C     Discretize domain and set up grid
      CALL DISCR(AXPT, AYPT, NNX, NNY, NEX, NEY, NOP)
C     Assemble stiffness matrix and RHS
      CALL AXB(SK, NOP, AXPT, AYPT, R3, NNX, NNY, NEX, NEY)
C     Solve system using Gaussian elimination
      CALL LEGS(SK, NNX*NNY, R3, X3, INDX)
C     Write results to files and screen
      DO I = 1, NNX*NNY
        WRITE(20, *) X3(I)
        WRITE(*, *) AXPT(I), AYPT(I), X3(I)
        WRITE(1, *) AXPT(I), AYPT(I), X3(I)
      ENDDO
C     End program
      STOP
      END

C     --------------------------------------------------------------------
      SUBROUTINE DISCR(AXPT, AYPT, NNX, NNY, NEX, NEY, NOP)
C     Sets up the grid and element connectivity
      IMPLICIT NONE
C     Input/Output variables
      REAL AXPT(441), AYPT(441)    ! Node x, y coordinates
      INTEGER NNX, NNY             ! Nodes in x, y directions
      INTEGER NEX, NEY             ! Elements in x, y directions
      INTEGER NOP(100,9)           ! Element-to-node mapping
C     Local variables
      REAL XFIRST, XLAST, DELTAX   ! x-domain parameters
      REAL YFIRST, YLAST, DELTAY   ! y-domain parameters
      INTEGER I, J, NEL, K, L      ! Loop indices
      INTEGER NNODE                ! Node index
C     Define domain: 0 <= x, y <= 1
      XFIRST = 0.0
      YFIRST = 0.0
      XLAST  = 1.0
      YLAST  = 1.0
      DELTAX = (XLAST - XFIRST) / FLOAT(NEX)
      DELTAY = (YLAST - YFIRST) / FLOAT(NEY)
C     Set node coordinates
      AXPT(1) = XFIRST
      AYPT(1) = YFIRST
      DO I = 1, NNX
        NNODE = (I-1)*NNY + 1
        AXPT(NNODE) = XFIRST + FLOAT(I-1)*DELTAX/2.0
        AYPT(NNODE) = YFIRST
        DO J = 2, NNY
          AXPT(NNODE+J-1) = AXPT(NNODE)
          AYPT(NNODE+J-1) = AYPT(NNODE) + FLOAT(J-1)*DELTAY/2.0
        ENDDO
      ENDDO
C     Set element-to-node connectivity (quadratic elements, 9 nodes each)
      NEL = 0
      DO I = 1, NEX
        DO J = 1, NEY
          NEL = NEL + 1
          DO K = 1, 3
            L = 3*K - 2
            NOP(NEL,L)   = NNY*(2*I+K-3) + 2*J - 1
            NOP(NEL,L+1) = NOP(NEL,L) + 1
            NOP(NEL,L+2) = NOP(NEL,L) + 2
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

C     --------------------------------------------------------------------
      SUBROUTINE AXB(SK, NOP, AXPT, AYPT, R3, NNX, NNY, NEX, NEY)
C     Assembles stiffness matrix (SK) and RHS (R3) with boundary conditions
      IMPLICIT NONE
C     Input/Output variables
      REAL SK(441,441)             ! Stiffness matrix
      REAL AXPT(441), AYPT(441)    ! Node coordinates
      REAL R3(441)                 ! RHS vector
      INTEGER NOP(100,9)           ! Element connectivity
      INTEGER NNX, NNY             ! Nodes in x, y
      INTEGER NEX, NEY             ! Elements in x, y
C     Local variables
      REAL PH(9)                   ! Basis functions
      REAL BC(441)                 ! Dirichlet boundary values
      INTEGER I, J, NELL, NE, NP   ! Loop and size variables
      INTEGER NCOD(441)            ! Boundary condition flags
      INTEGER NTOP(100)            ! Neumann boundary flags
C     Compute total elements and nodes
      NE = NEX * NEY
      NP = NNX * NNY
C     Initialize arrays
      DO I = 1, NP
        R3(I) = 0.0
        DO J = 1, NP
          SK(I,J) = 0.0
        ENDDO
      ENDDO
C     Neumann boundary conditions (top edge, y=1)
      DO I = 1, NE
        NTOP(I) = 0
      ENDDO
      DO I = NEX, NE, NEX
        NTOP(I) = 1
      ENDDO
C     Assemble matrix over all elements
      DO NELL = 1, NE
        CALL ABFIND(NELL, NOP, AXPT, AYPT, SK, R3, PH, NTOP)
      ENDDO
C     Dirichlet boundary conditions (y=0 fixed at 0)
      DO I = 1, NP
        NCOD(I) = 0
        BC(I) = 0.0
      ENDDO
      DO I = 1, NP-NNY+1, NNY
        NCOD(I) = 1
        BC(I) = 0.0
      ENDDO
C     Impose Dirichlet conditions
      DO I = 1, NP
        IF (NCOD(I) .EQ. 1) THEN
          R3(I) = BC(I)
          DO J = 1, NP
            SK(I,J) = 0.0
          ENDDO
          SK(I,I) = 1.0
        ENDIF
      ENDDO
      RETURN
      END

C     --------------------------------------------------------------------
      SUBROUTINE TSFUN(X, Y, PH, PHIC, PHIE)
C     Computes quadratic basis functions and their derivatives
      IMPLICIT NONE
C     Input/Output variables
      REAL X, Y                    ! Local coordinates (0 to 1)
      REAL PH(9)                   ! Basis functions
      REAL PHIC(9), PHIE(9)        ! Derivatives w.r.t. xi, eta
C     Local variables
      REAL L1X, L2X, L3X           ! x-direction basis
      REAL DL1X, DL2X, DL3X        ! x-direction derivatives
      REAL L1Y, L2Y, L3Y           ! y-direction basis
      REAL DL1Y, DL2Y, DL3Y        ! y-direction derivatives
C     Compute basis and derivatives in x
      L1X  = 2.0*X*X - 3.0*X + 1.0
      L2X  = -4.0*X*X + 4.0*X
      L3X  = 2.0*X*X - X
      DL1X = 4.0*X - 3.0
      DL2X = -8.0*X + 4.0
      DL3X = 4.0*X - 1.0
C     Compute basis and derivatives in y
      L1Y  = 2.0*Y*Y - 3.0*Y + 1.0
      L2Y  = -4.0*Y*Y + 4.0*Y
      L3Y  = 2.0*Y*Y - Y
      DL1Y = 4.0*Y - 3.0
      DL2Y = -8.0*Y + 4.0
      DL3Y = 4.0*Y - 1.0
C     Assemble basis functions
      PH(1) = L1X * L1Y
      PH(2) = L1X * L2Y
      PH(3) = L1X * L3Y
      PH(4) = L2X * L1Y
      PH(5) = L2X * L2Y
      PH(6) = L2X * L3Y
      PH(7) = L3X * L1Y
      PH(8) = L3X * L2Y
      PH(9) = L3X * L3Y
C     Derivatives w.r.t. xi
      PHIC(1) = L1Y * DL1X
      PHIC(2) = L2Y * DL1X
      PHIC(3) = L3Y * DL1X
      PHIC(4) = L1Y * DL2X
      PHIC(5) = L2Y * DL2X
      PHIC(6) = L3Y * DL2X
      PHIC(7) = L1Y * DL3X
      PHIC(8) = L2Y * DL3X
      PHIC(9) = L3Y * DL3X
C     Derivatives w.r.t. eta
      PHIE(1) = L1X * DL1Y
      PHIE(2) = L1X * DL2Y
      PHIE(3) = L1X * DL3Y
      PHIE(4) = L2X * DL1Y
      PHIE(5) = L2X * DL2Y
      PHIE(6) = L2X * DL3Y
      PHIE(7) = L3X * DL1Y
      PHIE(8) = L3X * DL2Y
      PHIE(9) = L3X * DL3Y
      RETURN
      END

C     --------------------------------------------------------------------
      SUBROUTINE ABFIND(NELL, NOP, AXPT, AYPT, SK, R3, PH, NTOP)
C     Computes element contributions to stiffness matrix and RHS
      IMPLICIT NONE
C     Input/Output variables
      REAL SK(441,441)             ! Stiffness matrix
      REAL AXPT(441), AYPT(441)    ! Node coordinates
      REAL R3(441)                 ! RHS vector
      REAL PH(9)                   ! Basis functions
      INTEGER NELL                 ! Element number
      INTEGER NOP(100,9)           ! Element connectivity
      INTEGER NTOP(100)            ! Neumann boundary flags
C     Local variables
      REAL PHX(9), PHY(9)          ! Derivatives in x, y
      REAL PHIC(9), PHIE(9)        ! Local derivatives
      REAL GP(3), WGP(3)           ! Gauss points and weights
      REAL X, Y, X1, X2, Y1, Y2    ! Mapped coordinates and Jacobians
      REAL DETT                    ! Jacobian determinant
      INTEGER I, K, M, N, M1, N1   ! Loop indices
      INTEGER NGL(9)               ! Global node numbers
C     Set global node numbers for this element
      DO I = 1, 9
        NGL(I) = NOP(NELL,I)
      ENDDO
C     Gauss points and weights (3-point rule)
      GP(1)  = (1.0 - SQRT(3.0/5.0)) / 2.0
      GP(2)  = 0.5
      GP(3)  = (1.0 + SQRT(3.0/5.0)) / 2.0
      WGP(1) = 5.0 / 18.0
      WGP(2) = 8.0 / 18.0
      WGP(3) = 5.0 / 18.0
C     Loop over Gauss points
      DO I = 1, 3
        DO K = 1, 3
          CALL TSFUN(GP(I), GP(K), PH, PHIC, PHIE)
C         Compute mapping and Jacobian
          X  = 0.0
          Y  = 0.0
          X1 = 0.0
          X2 = 0.0
          Y1 = 0.0
          Y2 = 0.0
          DO N = 1, 9
            X  = X  + AXPT(NGL(N)) * PH(N)
            Y  = Y  + AYPT(NGL(N)) * PH(N)
            X1 = X1 + AXPT(NGL(N)) * PHIC(N)
            X2 = X2 + AXPT(NGL(N)) * PHIE(N)
            Y1 = Y1 + AYPT(NGL(N)) * PHIC(N)
            Y2 = Y2 + AYPT(NGL(N)) * PHIE(N)
          ENDDO
          DETT = X1*Y2 - X2*Y1
C         Compute physical derivatives
          DO N = 1, 9
            PHX(N) = (Y2*PHIC(N) - Y1*PHIE(N)) / DETT
            PHY(N) = (X1*PHIE(N) - X2*PHIC(N)) / DETT
          ENDDO
C         Assemble stiffness and RHS
          DO M = 1, 9
            M1 = NGL(M)
            R3(M1) = R3(M1) + WGP(I)*WGP(K)*DETT*PH(M)*0.0
            DO N = 1, 9
              N1 = NGL(N)
              SK(M1,N1) = SK(M1,N1) - WGP(I)*WGP(K)*DETT*
     &                    (PHX(M)*PHX(N) + PHY(M)*PHY(N))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C     Neumann condition on top boundary (y=1)
      IF (NTOP(NELL) .EQ. 1) THEN
        DO I = 1, 3
          CALL TSFUN(GP(I), 1.0, PH, PHIC, PHIE)
          X  = 0.0
          X1 = 0.0
          DO K = 1, 9
            X  = X  + AXPT(NGL(K)) * PH(K)
            X1 = X1 + AXPT(NGL(K)) * PHIC(K)
          ENDDO
          DO M1 = 3, 9, 3
            R3(NGL(M1)) = R3(NGL(M1)) - WGP(I)*X1*PH(M1)*5.0
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END

C     --------------------------------------------------------------------
      SUBROUTINE LEGS(A, N, B, X3, INDX)
C     Solves A*X = B using Gaussian elimination with partial pivoting
      IMPLICIT NONE
C     Input/Output variables
      INTEGER N                    ! Matrix size
      REAL A(441,441)              ! Coefficient matrix
      REAL B(441)                  ! RHS vector
      REAL X3(441)                 ! Solution vector
      INTEGER INDX(441)            ! Pivot indices
C     Local variables
      INTEGER I, J
C     Perform elimination
      CALL ELGS(A, N, INDX)
      DO I = 1, N-1
        DO J = I+1, N
          B(INDX(J)) = B(INDX(J)) - A(INDX(J),I)*B(INDX(I))
        ENDDO
      ENDDO
C     Back substitution
      X3(N) = B(INDX(N)) / A(INDX(N),N)
      DO I = N-1, 1, -1
        X3(I) = B(INDX(I))
        DO J = I+1, N
          X3(I) = X3(I) - A(INDX(I),J)*X3(J)
        ENDDO
        X3(I) = X3(I) / A(INDX(I),I)
      ENDDO
      RETURN
      END

C     --------------------------------------------------------------------
      SUBROUTINE ELGS(A, N, INDX)
C     Performs partial-pivoting Gaussian elimination
      IMPLICIT NONE
C     Input/Output variables
      INTEGER N                    ! Matrix size
      REAL A(441,441)              ! Coefficient matrix
      INTEGER INDX(441)            ! Pivot indices
C     Local variables
      INTEGER I, J, K, ITMP
      REAL C1, PI, PI1, PJ
      REAL C(441)                  ! Scaling factors
C     Initialize pivot indices
      DO I = 1, N
        INDX(I) = I
      ENDDO
C     Compute scaling factors
      DO I = 1, N
        C1 = 0.0
        DO J = 1, N
          C1 = AMAX1(C1, ABS(A(I,J)))
        ENDDO
        C(I) = C1
      ENDDO
C     Eliminate with pivoting
      DO J = 1, N-1
        PI1 = 0.0
        DO I = J, N
          PI = ABS(A(INDX(I),J)) / C(INDX(I))
          IF (PI .GT. PI1) THEN
            PI1 = PI
            K   = I
          ENDIF
        ENDDO
        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
        DO I = J+1, N
          PJ = A(INDX(I),J) / A(INDX(J),J)
          A(INDX(I),J) = PJ
          DO K = J+1, N
            A(INDX(I),K) = A(INDX(I),K) - PJ*A(INDX(J),K)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
