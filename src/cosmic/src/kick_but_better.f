      SUBROUTINE ChangeBasis(Vector, ThetaE, PhiE, PsiE, Result)
* Redefine a vector from one coordinate basis to another using Euler Angles
* Vector is the input vector in the new basis (X', Y', Z')
* Result is the transformed vector in the original basis (X, Y, Z)
* ThetaE, PhiE, PsiE are the Euler angles
  
          DOUBLE PRECISION Vector(3), Result(3)
          DOUBLE PRECISION ThetaE, PhiE, PsiE
          DOUBLE PRECISION cTheta, cPhi, cPsi, sTheta, sPhi, sPsi
          DOUBLE PRECISION rotationMatrix(3,3)
          INTEGER i, j
    
* define trigonometric values
          cTheta = COS(ThetaE)
          sTheta = SIN(ThetaE)
          cPhi   = COS(PhiE)
          sPhi   = SIN(PhiE)
          cPsi   = COS(PsiE)
          sPsi   = SIN(PsiE)
    
* define the Rotation Matrix
          rotationMatrix(1,1) = cPhi * cPsi - sPhi * cTheta * sPsi
          rotationMatrix(1,2) = -cPhi * sPsi - sPhi * cTheta * cPsi
          rotationMatrix(1,3) = sTheta * sPhi
          rotationMatrix(2,1) = sPhi * cPsi + cPhi * cTheta * sPsi
          rotationMatrix(2,2) = -sPhi * sPsi + cPhi * cTheta * cPsi
          rotationMatrix(2,3) = -sTheta * cPhi
          rotationMatrix(3,1) = sTheta * sPsi
          rotationMatrix(3,2) = sTheta * cPsi
          rotationMatrix(3,3) = cTheta
    
* initialize the result to zero
          DO i = 1, 3
             Result(i) = 0.0D0
          END DO
    
* apply rotation to the vector
          DO i = 1, 3
             DO j = 1, 3
                Result(i) = Result(i) + Vector(j) * rotationMatrix(i, j)
             END DO
          END DO
    
          RETURN
          END
