C************************************************************
C  Program for Natural Convection From Heated Slender Body  *
C			                							  *	
C******************* MAIN PROGRAM ***************************
      PROGRAM MAIN
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON SHI(900,700),VORTI(900,700),ERR(900,700),U(900,700),
     +V(900,700),XU(900),YV(700),X(900),Y(700),T(900,700)
      COMMON XL,YL,H,RF,RE,RI,EMAX,RA,PR
      COMMON W,YL1,KBS1,KBE1,KBHS1,KBHE1,KBHS11,KBHE11,KW,KL,KL1,
     +KBS11,KBE11
      COMMON /INDEX/M,N,K,ITER,LAST,IMAX,JMAX,INTV,KW1
C***********************************************
      CALL GRID
      CALL START
	!CALL RESTART
   5  ITER=ITER+1
      CALL BOUND
      CALL SOLVE

	 IF(ITER.EQ.100000)CALL RESULT
		IF(ITER.EQ.200000)CALL RESULT
	IF(ITER.EQ.300000)CALL RESULT

      WRITE(*,30)ITER,EMAX,SHI(11,35),VORTI(14,35),T(177,47)
      IF(ITER.LT.LAST)GO TO 5
      CALL RESULT
  30  FORMAT(2X,I6,4(1X,E16.10)) 
      STOP
      END
C**************************************************
C*************SUBPROGRAM FOR GRID GENERATION*******
      SUBROUTINE GRID 
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON SHI(900,700),VORTI(900,700),ERR(900,700),U(900,700),
     +V(900,700),XU(900),YV(700),X(900),Y(700),T(900,700)
      COMMON XL,YL,H,RF,RE,RI,EMAX,RA,PR
      COMMON W,YL1,KBS1,KBE1,KBHS1,KBHE1,KBHS11,KBHE11,KW,KL,KL1,
     +KBS11,KBE11
      COMMON /INDEX/M,N,K,ITER,LAST,IMAX,JMAX,INTV,KW1
C**************************************************
      WRITE(*,*)'MAXIMUM No.OF ITERATIONS'
      READ(*,*)LAST
      ITER=0000000
      WRITE(*,*)'RELAXATION FACTOR'
      READ(*,*)RF
      XL=8.0
	W=1.0
	H=W
      YL=4.0
	M=480
	K=240
	RE=1.0
	RA=50.0
	PR=0.71
      POWERX=1.0
      POWERY=1.0 
      XU(1)=0.0
      YV(1)=0.0
      X(1)=0.0
      Y(1)=0.0
      DO I=2,M
      XU(I)=XU(1)+((FLOAT(I-1)/FLOAT(M-1))**POWERX)*XL
      ENDDO
      DO I=2,M
      X(I)=XU(I)-XU(I-1)
      ENDDO
      DO J=2,K
      YV(J)=YV(1)+((FLOAT(J-1)/FLOAT(K-1))**POWERY)*YL
      ENDDO
      DO J=2,K
      Y(J)=YV(J)-YV(J-1)
      ENDDO
	RETURN
	END
C*************************************************
C************SUBPROGRAM FOR STARTING**************
      SUBROUTINE START  
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON SHI(900,700),VORTI(900,700),ERR(900,700),U(900,700),
     +V(900,700),XU(900),YV(700),X(900),Y(700),T(900,700)
      COMMON XL,YL,H,RF,RE,RI,EMAX,RA,PR
      COMMON W,YL1,KBS1,KBE1,KBHS1,KBHE1,KBHS11,KBHE11,KW,KL,KL1,
     +KBS11,KBE11
      COMMON /INDEX/M,N,K,ITER,LAST,IMAX,JMAX,INTV,KW1
C*************************************************
      DO I=1,M
      DO J=1,K
      VORTI(I,J)=0.0
      SHI(I,J)=0.0
      ERR(I,J)=0.0
      U(I,J)=0.0
	V(I,J)=0.0
	T(I,J)=0.0
      ENDDO
      ENDDO 
      RETURN
      END
C*****************************************************
C***********SUBPROGRM FOR BOUNDARY CONDITION**********
      SUBROUTINE BOUND  
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON SHI(900,700),VORTI(900,700),ERR(900,700),U(900,700),
     +V(900,700),XU(900),YV(700),X(900),Y(700),T(900,700)
      COMMON XL,YL,H,RF,RE,RI,EMAX,RA,PR
      COMMON W,YL1,KBS1,KBE1,KBHS1,KBHE1,KBHS11,KBHE11,KW,KL,KL1,
     +KBS11,KBE11
      COMMON /INDEX/M,N,K,ITER,LAST,IMAX,JMAX,INTV,KW1
C********************************************************

      K1=K-1
      M1=M-1

C********BOUNDARY CONDITION FOR BOTTOM WALL *********

      DO I=1,M
      SHI(I,1)=-0.5
      VORTI(I,1)=-(.05*RE)/2.0 
      U(I,1)=0.0
      V(I,1)=0.0
	T(I,1)=T(I,2)
      ENDDO
	 
C**********BOUNDERY CONDITION FOR LEFT BOUNDARY*****
      
	DO J=2,K      
	SHI(1,J)=SHI(2,J)      
	VORTI(1,J)=VORTI(2,J)	  
	U(1,J)=U(2,J)
      V(1,J)=V(2,J)
      T(1,J)=T(2,J)
      ENDDO

C***********BOUNDARY CONDITION FOR TOP BOUNDARY ******

      DO I=2,M
	SHI(I,K)=0.5
	VORTI(I,K)=(.05*RE)/2.0 
	U(I,K)=0.0
	V(I,K)=0.0
	T(I,K)=T(I,K1)
      END DO

C**********BOUNDERY CONDITION FOR OUTLET BOUNDARY *****

      DO J=2,K1
      SHI(M,J)=SHI(M1,J)
      VORTI(M,J)=VORTI(M1,J)
	U(M,J)=U(M1,J)
      V(M,J)=V(M1,J)
      T(M,J)=T(M1,J)
      ENDDO

C*************BOUNDARY CONDITION FOR SLENDER BODY *******

      KBS1=210
	KBE1=KBS1+60
	KBHS1=117
	KBHE1=KBHS1+6
	KBS11=KBS1-1
	KBE11=KBE1+1
	KBHS11=KBHS1-1
	KBHE11=KBHE1+1

	DO I=KBS1,KBE1
	DO J=KBHS1,KBHE1
	SHI(I,J)=-0.5
	VORTI(I,J)=-(.05*RE)/2
	T(I,J)=1.0
	U(I,J)=0.0
	V(I,J)=0.0
	ENDDO
	ENDDO

	RETURN 
	END
C*********************************************************
C************SUBPROGRAM FOR CALCULATION*******************
C*********************************************************
      SUBROUTINE SOLVE  
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON SHI(900,700),VORTI(900,700),ERR(900,700),U(900,700),
     +V(900,700),XU(900),YV(700),X(900),Y(700),T(900,700)
      COMMON XL,YL,H,RF,RE,RI,EMAX,RA,PR
      COMMON W,YL1,KBS1,KBE1,KBHS1,KBHE1,KBHS11,KBHE11,KW,KL,KL1,
     +KBS11,KBE11
      COMMON /INDEX/M,N,K,ITER,LAST,IMAX,JMAX,INTV,KW1
C**********************************************************

      M1=M-1
      K1=K-1
C***********************************************************
      DO I=2,M1
      DO J=2,K1
      STF=SHI(I,J)

      TS41=2.0/X(I)**2
      TS42=2.0/Y(J)**2
	TS4=TS41+TS42

      TS51=SHI(I+1,J)/X(I)**2
	TS52=SHI(I-1,J)/X(I)**2
	TS5=TS51+TS52

      TS61=SHI(I,J+1)/Y(J)**2
	TS62=SHI(I,J-1)/Y(J)**2      
	TS6=TS61+TS62

C***************************************************************
C************CALCULATES STREAM FUNCTIONS FOR EACH NODES**********
C******************************************************************
      VORTIH=(VORTI(I,J)*H**2)/SHI(1,1)

	SHI(I,J)=(VORTIH+TS5+TS6)/TS4

C***************************************************************
C************UNDER RELAXED THE STREAM FUNCTION FOR EACH NODES*******
C*******************************************************************
      SHI(I,J)=RF*SHI(I,J)+(1.0-RF)*STF

C**************************************************************
      ENDDO
      ENDDO
C*********************************************************
      DO I=2,M1
      DO J=2,K1

	RFU=0.7
	RFV=0.5
      
	OLDU=U(I,J)
	OLDV=V(I,J)

      U(I,J)=(SHI(I,J+1)-SHI(I,J-1))/(2.0*Y(J))
      V(I,J)=(SHI(I-1,J)-SHI(I+1,J))/(2.0*X(I))

	U(I,J)=RFU*U(I,J)+(1.0-RFU)*OLDU
	V(I,J)=RFV*V(I,J)+(1.0-RFV)*OLDV

      ENDDO
      ENDDO

C*****************************************************************

      DO I=2,M1
      DO J=2,K1

      OMEG=VORTI(I,J)

      PRX=PR/X(I)**2
	PRY=PR/Y(J)**2

	PRX2=2.0*PRX
	PRY2=2.0*PRY

	UX=U(I,J)/X(I)
	VY=V(I,J)/Y(J)
	TV1=PRX2+PRY2-UX-VY
     
      TV211=VORTI(I+1,J)/PRX
	TV212=VORTI(I-1,J)/PRX
      TV21=TV211+TV212

      TV221=VORTI(I,J+1)/PRY
	TV222=VORTI(I,J-1)/PRY
      TV22=TV221+TV222
	   
	TV23=(RA*PR*(T(I+1,J)-T(I,J)))/X(I)	        
      TV2=TV21+TV22-TV23
	
      TV31=(U(I+1,J)*VORTI(I+1,J))/X(I)
	TV32=(V(I,J+1)*VORTI(I,J+1))/Y(J)
      TV3=TV31+TV32 
		   
C******************************************************************
C*********CALCULATES VORTICITY FOR EACH NODES************************   
C*********************************************************************
      
	VORTI(I,J)=(TV2-TV3)/TV1

C**********************************************************************
C*********UNDER RELAXED THE VORTICITY FOR EACH NODES******************
C*********************************************************************
      
	VORTI(I,J)=RF*VORTI(I,J)+(1.0-RF)*OMEG
      IF(VORTI(I,J).EQ.0.0)GO TO 25
      IF(ITER.GT.6)THEN
      ERR(I,J)=(ABS((VORTI(I,J)-OMEG)))/(ABS(VORTI(I,J)))   
      ENDIF
  25  CONTINUE
      ENDDO
      ENDDO      
C***************************************************************
C********************   TEMPERATURE    *************************
C***************************************************************
      DO I=2,M1
	DO J=2,K1

      TMF=T(I,J)

	TS711=2.0/X(I)**2
	TS712=2.0/Y(J)**2
      TS71=TS711+TS712
      UX=U(I,J)/X(I)
	VY=V(I,J)/Y(J)
	TS72=UX+VY
	TS7=TS71-TS72

      TS811=T(I+1,J)/X(I)**2
	TS812=T(I-1,J)/X(I)**2
      TS821=T(I,J+1)/Y(J)**2
      TS822=T(I,J-1)/Y(J)**2
	TS8=TS811+TS812+TS821+TS822

	TS91=((U(I+1,J))*(T(I+1,J)))/X(I)
	TS92=((V(I,J+1))*(T(I,J+1)))/Y(J)
	TS9=TS91+TS92

C****************TEMPERATURE FOR EACH NODE*********************

      T(I,J)=(TS8-TS9)/TS7

C***********UNDER RELAXED THE TEMPERATURE FOR EACH NODE********

      T(I,J)=RF*T(I,J)+(1.0-RF)*TMF

	ENDDO
	ENDDO
	EMAX=0.0
      DO I=2,M1
      DO J=2,K1
      EMAX=MAX(EMAX,ERR(I,J)) 
      IF(EMAX.EQ.ERR(I,J))THEN
      IMAX=I
      JMAX=J
      ENDIF
      END DO
      END DO

	RETURN
      END
C**********************************************************
C*********SUBPROGRAM FOR OUTPUT*********************************
      SUBROUTINE RESULT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON SHI(900,700),VORTI(900,700),ERR(900,700),U(900,700),
     +V(900,700),XU(900),YV(700),X(900),Y(700),T(900,700)
      COMMON XL,YL,H,RF,RE,RI,EMAX,RA,PR
      COMMON W,YL1,KBS1,KBE1,KBHS1,KBHE1,KBHS11,KBHE11,KW,KL,KL1,
     +KBS11,KBE11
      COMMON /INDEX/M,N,K,ITER,LAST,IMAX,JMAX,INTV,KW1
C*********************************************************

	OPEN(UNIT=13,FILE='RESULT.TXT')

C*******************************************************
	K1=K-1
	M1=M-1

	!PI=XU(KBS2)-XU(KBS1)
	RIBH=YV(KBHE1)-YV(KBHS1)
	RIBW=XU(KBE1)-XU(KBS1)
	AR=RIBW/RIBH


	WRITE(13,201)ITER
	WRITE(13,203)RF
	WRITE(13,205)RA
	WRITE(13,207)PR
	WRITE(13,211)XL
	WRITE(13,213)YL
	WRITE(13,215)XU(KBS1)
	WRITE(13,216)RIBH
	WRITE(13,219)RIBW
	WRITE(13,223)AR


C*********************************************************
201	FORMAT(4X,'No. of Iterations =',I8)
203	FORMAT(4X,'Relaxation Factor =',F8.4)
205	FORMAT(4X,'Rayleigh No. =',F10.6)
207	FORMAT(4X,'Prandtl No. =',F8.4)
211	FORMAT(4X,'Domain Length =',F8.4)
213	FORMAT(4X,'Domain Height =',F8.4)
215	FORMAT(4X,'Object Location =',F8.4)
216	FORMAT(4X,'Object Height =',F8.4)
219	FORMAT(4X,'Object Width =',F8.4)
223	FORMAT(4X,'Aspect Ratio =',F8.4,/) 

C*********************************************************

	OPEN(UNIT=12,FILE='INLETCOND.TXT')

C*********************************************************
	DO J=1,K
	WRITE(12,102)J,YV(J),U(1,J),U(3,J),U(7,J),V(7,J)
	ENDDO
102   FORMAT(3X,I3,5(2X,F9.5))
C*********************************************************

	OPEN(UNIT=15,FILE='PARAM.TXT')

C*********************************************************
	DO I=1,M
	DO J=1,K
	WRITE(15,217)SHI(I,J),VORTI(I,J),U(I,J),V(I,J),T(I,J)
	END DO
	END DO
217   FORMAT(2X,5(1X,F12.8))
	CLOSE(15)

C*********************************************************

      OPEN(UNIT=301,FILE='UVel Profile.TXT')
      OPEN(UNIT=303,FILE='Temp Profile.TXT')

	DO J=1,K
	 WRITE(301,317)U(21,J),U(39,J),U(56,J),U(82,J),U(106,J),YV(J)
	WRITE(303,317)T(21,J),T(39,J),T(56,J),T(82,J),T(106,J),YV(J)
	 END DO
317   FORMAT(2X,6(1X,F12.8))
C*********************************************************

	WRITE(13,229)XU(21)
	WRITE(13,231)XU(39)
	WRITE(13,233)XU(56)
	WRITE(13,235)XU(82)
	WRITE(13,237)XU(106)

229	FORMAT(4X,'XU(21) =',F8.4)
231	FORMAT(4X,'XU(39) =',F8.4)
233	FORMAT(4X,'XU(56) =',F8.4)
235	FORMAT(4X,'XU(82) =',F8.4)
237	FORMAT(4X,'XU(106) =',F8.4,/) 

C*********************************************************

      OPEN(UNIT=103,FILE='SHI.TXT')
	OPEN(UNIT=105,FILE='VORTI.TXT')
	OPEN(UNIT=107,FILE='Temp.TXT')
      OPEN(UNIT=109,FILE='UVelo.TXT')
	OPEN(UNIT=111,FILE='VVelo.TXT')

C*********************************************************
	DO I=1,M
	DO J=1,K
	WRITE(103,221)XU(I),YV(J),SHI(I,J)
      WRITE(105,221)XU(I),YV(J),VORTI(I,J)
	WRITE(107,221)XU(I),YV(J),T(I,J)
	WRITE(109,221)XU(I),YV(J),U(I,J)
	WRITE(111,221)XU(I),YV(J),V(I,J)
221   FORMAT(2X,3(1X,E16.10))
	ENDDO
	ENDDO

C*******************************************************
C***********  Slender Body Nusselt Number   ********
C*******************************************************

	OPEN(UNIT=115,FILE='BLFNUX1.TXT')
	OPEN(UNIT=116,FILE='BRFNUX1.TXT')
	OPEN(UNIT=117,FILE='BTFNUX1.TXT')
	OPEN(UNIT=118,FILE='BBFNUX1.TXT')


	BLFNUA1=0.0
	BRFNUA1=0.0

	DO J=KBHS1,KBHE1
	BLFNUX1=(T(KBS1,J)-T(KBS11,J))/X(KBS11)
	BRFNUX1=(T(KBE1,J)-T(KBE11,J))/X(KBE11)
	BLFNUA1=BLFNUA1+(BLFNUX1*Y(J))
	BRFNUA1=BRFNUA1+(BRFNUX1*Y(J))
	WRITE(115,361)YV(J),BLFNUX1
	WRITE(116,361)YV(J),BRFNUX1
	END DO

	BLFNUA2=BLFNUA1/RIBH
	BRFNUA2=BRFNUA1/RIBH

	BTFNUA1=0.0
	BBFNUA1=0.0


	DO I=KBS1,KBE1
	BXW1=XU(I)-XU(KBS1)
	BBFNUX1=(T(I,KBHS1)-T(I,KBHS11))/Y(KBHS11)
	BBFNUA1=BBFNUA1+(BBFNUX1*X(I))
	BTFNUX1=(T(I,KBHE1)-T(I,KBHE11))/Y(KBHE11)
	BTFNUA1=BTFNUA1+(BTFNUX1*X(I))
	WRITE(117,361)BXW1,BTFNUX1
	WRITE(118,361)BXW1,BBFNUX1
	END DO

	BBFNUA2=BBFNUA1/RIBW
	BTFNUA2=BTFNUA1/RIBW

      !WRITE(13,225)BRIBH1
	!WRITE(13,227)BRIBW1
225	FORMAT(4X,'First Rib Height at Bottom =',F8.4)
227	FORMAT(4X,'First Rib Width at Bottom =',F8.4,/)

	WRITE(13,401)BLFNUA2
	WRITE(13,403)BRFNUA2
	WRITE(13,405)BTFNUA2
	WRITE(13,406)BBFNUA2

 401	FORMAT(4X,'Left Face Avg. Nu =',F12.8)
 403	FORMAT(4X,'Right Face Avg. Nu =',F12.8)
 405	FORMAT(4X,'Top Face Avg. Nu =',F12.8)
 406	FORMAT(4X,'Bottom Face Avg. Nu =',F12.8,/)

      CLOSE(115)
	  CLOSE(116)
	     CLOSE(117)
	       CLOSE(118)

C*******************************************************
C***********  Second Bottom Rib Nusselt Number   *******
C*******************************************************

	!OPEN(UNIT=119,FILE='BLFNUX2.TXT')
	!OPEN(UNIT=121,FILE='BRFNUX2.TXT')
	!OPEN(UNIT=123,FILE='BTFNUX2.TXT')

	!KBS21=KBS2-1
	!KBE21=KBE2+1
	!BRIBH2=YV(KBH2)-YV(1)

	!DO J=2,KBH2
	!BLFNUX2=(T(KBS2,J)-T(KBS21,J))/X(KBS21)
	!BRFNUX2=(T(KBE2,J)-T(KBE21,J))/X(KBE21)
	!BLFNUA2=(BLFNUX2*Y(J))/BRIBH2
	!BRFNUA2=(BRFNUX2*Y(J))/BRIBH2
	!WRITE(119,361)YV(J),BLFNUX2
	!WRITE(121,361)YV(J),BRFNUX2
	!END DO


	!KBH21=KBH2+1
	!BRIBW2=XU(KBE2)-XU(KBS2)
	!DO I=KBS2,KBE2
	!BXW2=XU(I)-XU(KBS2)
	!BTFNUX2=(T(I,KBH2)-T(I,KBH21))/Y(KBH21)
	!BTFNUA2=(BTFNUX2*X(I))/BRIBW2
	!WRITE(123,361)BXW2,BTFNUX2
	!END DO

      !WRITE(13,239)BRIBH2
	!WRITE(13,241)BRIBW2
 239	FORMAT(4X,'Second Rib Height at Bottom =',F8.4)
 241	FORMAT(4X,'Second Rib Width at Bottom =',F8.4,/)

	!WRITE(13,407)BLFNUA2
	!WRITE(13,409)BRFNUA2
	!WRITE(13,411)BTFNUA2
 407	FORMAT(4X,'Second Bottom Left Face Avg. Nu =',F8.4)
 409	FORMAT(4X,'Second Bottom Right Face Avg. Nu =',F8.4)
 411	FORMAT(4X,'Second Bottom Top Face Avg. Nu =',F8.4,/)


C*******************************************************
C***********  First Top Rib Nusselt Number   ********
C*******************************************************

	!OPEN(UNIT=125,FILE='TLFNUX1.TXT')
	!OPEN(UNIT=127,FILE='TRFNUX1.TXT')
	!OPEN(UNIT=129,FILE='TTFNUX1.TXT')

	!KTS11=KTS1-1
	!KTE11=KTE1+1
	!KKTOP1=K-KTH1
	!TRIBH1=YV(K1)-YV(KKTOP1)

	!DO J=2,KTH1
	!JJ=K-J
	!YVH1=YV(K1)-YV(JJ)
	!TLFNUX1=(T(KTS1,JJ)-T(KTS11,JJ))/X(KTS11)
	!TRFNUX1=(T(KTE1,JJ)-T(KTE11,JJ))/X(KTE11)
	!TLFNUA1=(TLFNUX1*Y(JJ))/TRIBH1
	!TRFNUA1=(TRFNUX1*Y(JJ))/TRIBH1
	!WRITE(125,361)YVH1,TLFNUX1
	!WRITE(127,361)YVH1,TRFNUX1
	!END DO


	!KTH11=KTH1+1
	!TRIBW1=XU(KTE1)-XU(KTS1)
	!DO I=KTS1,KTE1
	!TXW1=XU(I)-XU(KTS1)
	!TTFNUX1=(T(I,KTH1)-T(I,KTH11))/Y(KTH11)
	!TTFNUA1=(TTFNUX1*X(I))/TRIBW1
	!WRITE(129,361)TXW1,TTFNUX1
	!END DO

      !WRITE(13,243)TRIBH1
	!WRITE(13,245)TRIBW1
243	FORMAT(4X,'First Rib Height at Top =',F8.4)
245	FORMAT(4X,'First Rib Width at Top =',F8.4,/)

	!WRITE(13,413)TLFNUA1
	!WRITE(13,415)TRFNUA1
	!WRITE(13,417)TTFNUA1
 413	FORMAT(4X,'First Top Left Face Avg. Nu =',F8.4)
 415	FORMAT(4X,'First Top Right Face Avg. Nu =',F8.4)
 417	FORMAT(4X,'First Top Top Face Avg. Nu =',F8.4,/)

C*******************************************************
C***********  Second Top Rib Nusselt Number      *******
C*******************************************************

	!OPEN(UNIT=131,FILE='TLFNUX2.TXT')
	!OPEN(UNIT=133,FILE='TRFNUX2.TXT')
	!OPEN(UNIT=135,FILE='TTFNUX2.TXT')

	!KTS21=KTS2-1
	!KTE21=KTE2+1
	!KKTOP2=K-KTH2
	!TRIBH2=YV(K1)-YV(KKTOP2)

	!DO J=2,KTH2
	!JJ=K-J
	!YVH2=YV(K1)-YV(JJ)
	!TLFNUX2=(T(KTS2,JJ)-T(KTS21,JJ))/X(KTS21)
	!TRFNUX2=(T(KTE2,JJ)-T(KTE21,JJ))/X(KTE21)
	!TLFNUA2=(TLFNUX2*Y(JJ))/TRIBH2
	!TRFNUA2=(TRFNUX2*Y(JJ))/TRIBH2
	!WRITE(131,361)YVH2,TLFNUX2
	!WRITE(133,361)YVH2,TRFNUX2
	!END DO


	!KTH21=KTH2+1
	!TRIBW2=XU(KTE2)-XU(KTS2)
	!DO I=KTS2,KTE2
	!TXW2=XU(I)-XU(KTS2)
	!TTFNUX2=(T(I,KTH2)-T(I,KTH21))/Y(KTH21)
	!TTFNUA2=(TTFNUX2*X(I))/TRIBW2
	!WRITE(135,361)TXW2,TTFNUX2
	!END DO

      !WRITE(13,247)TRIBH2
	!WRITE(13,249)TRIBW2
247	FORMAT(4X,'Second Rib Height at Top =',F8.4)
249	FORMAT(4X,'Second Rib Width at Top =',F8.4,/)

	!WRITE(13,419)TLFNUA2
	!WRITE(13,421)TRFNUA2
	!WRITE(13,423)TTFNUA2
 419	FORMAT(4X,'Second Top Left Face Avg. Nu =',F8.4)
 421	FORMAT(4X,'Second Top Right Face Avg. Nu =',F8.4)
 423	FORMAT(4X,'Second Top Top Face Avg. Nu =',F8.4,/)

	BRIBAVG1=(BLFNUA2+BRFNUA2+BTFNUA2+BBFNUA2)/4.0
	!BRIBAVG2=(BLFNUA2+BRFNUA2+BTFNUA2)/3.0
	!TRIBAVG1=(TLFNUA1+TRFNUA1+TTFNUA1)/3.0
	!TRIBAVG2=(TLFNUA2+TRFNUA2+TTFNUA2)/3.0

	WRITE(13,425)BRIBAVG1
	!WRITE(13,427)BRIBAVG2
	!WRITE(13,429)TRIBAVG1
	!WRITE(13,431)TRIBAVG2

 425	FORMAT(4X,'Object Avg. Nu =',F12.8,/)
 427	FORMAT(4X,'Second Bottom Rib Avg. Nu =',F8.4)
 429	FORMAT(4X,'First Top Rib Avg. Nu =',F8.4)
 431	FORMAT(4X,'Second Top Rib Avg. Nu =',F8.4,/)


361	FORMAT(4X,2(2X,F14.8))

C*******************************************************
C***********       Heating Efficiency   ****************
C*******************************************************
   
	TIREL=0.0
      TINREL=0.0
      BIREL=0.0
      BINREL=0.0

      DO I=KBS1,KBE1
	TINREL=TINREL+X(I)*V(I,KBHE11)*T(I,KBHE11)
	TIREL=TIREL+X(I)*V(I,KBHE11)
	BINREL=BINREL+X(I)*V(I,KBHS11)*T(I,KBHS11)
	BIREL=BIREL+X(I)*V(I,KBHS11)
	END DO

      LIREL=0.0
	RIREL=0.0
	LINREL=0.0
	RINREL=0.0


      DO J=KBHS1,KBHE1
	LINREL=LINREL+Y(J)*U(KBS11,J)*T(KBS11,J)
	LIREL=LIREL+Y(J)*U(KBS11,J)
	RINREL=RINREL+Y(J)*U(KBE11,J)*T(KBE11,J)
	RIREL=RIREL+Y(J)*U(KBE11,J)
	END DO

      TI=TIREL+BIREL+LIREL+RIREL
      TIN=TINREL+BINREL+LINREL+RINREL

      KOUR=M-3
		
	TOUL=0.0
	TOUTL=0.0
	TOUR=0.0
	TOUTR=0.0
	K11=K/2+1

	DO J=2,K11
	TOUTL=TOUTL+Y(J)*U(3,J)*T(3,J)
 	TOUL=TOUL+Y(J)*U(3,J)
      TOUTR=TOUTR+Y(J)*U(KOUR,J)*T(KOUR,J)
 	TOUR=TOUR+Y(J)*U(KOUR,J)
	END DO

	ATIN=TIN/TI
	ATOUT=(TOUTR/TOUR)+(TOUTL/TOUL)
	RATIO=ATIN/ATOUT
	EFF=(1.0-RATIO)*100

	HREL=ATIN

	WRITE(13,461)HREL
	WRITE(13,463)TOUTL
	WRITE(13,465)TOUTR

 461	FORMAT(4X,'Amount of Heat Released from the Object =',F12.8)
 463	FORMAT(4X,'Amount of Heat Going out from Left Side =',F12.8)
 465	FORMAT(4X,'Amount of Heat Going out from Right Side =',F12.8,/)

	WRITE(13,333)EFF
 333	FORMAT(4X,'Heating Efficiency =',F8.4,' %')


C********************************************************* 
	OPEN(UNIT=33,FILE='ERROR.TXT')
      REWIND(33)
C*********************************************************
      WRITE(33,309)ITER
309   FORMAT(3X,'ITERATION=',I6,//)
      WRITE(33,311)RF,RA,PR  
311   FORMAT(3X,'RF= ',F6.2,' RA= ',F8.2,' PR= ',F8.2)
	WRITE(33,322) 
322   FORMAT(//,'I, J ,X(I),Y(J) , ERR(I,J) , SHI(I.J) , VORTI(I,J)',/) 
      DO I=1,M
      DO J=1,K
      WRITE(33,52)I,J,X(I),Y(J),ERR(I,J),SHI(I,J),VORTI(I,J)
      ENDDO
      ENDDO
      WRITE(33,*)
      WRITE(33,*)'STREAM FUNCTIONS ON GRIDS'
      WRITE(33,120)((SHI(I,J),I=1,M),J=K,1,-1)
      WRITE(33,*)
      WRITE(33,*)'VORTICITY ON GRIDS'
      WRITE(33,120)((VORTI(I,J),I=1,M),J=K,1,-1)
      WRITE(33,80)
      WRITE(33,120)((U(I,J),I=1,M),J=K,1,-1)
      WRITE(33,81)
      WRITE(33,120)((V(I,J),I=1,M),J=K,1,-1)
 80   FORMAT(//,'=====U-COMPONENT OF VELOCITY=====',/)
 81   FORMAT(//,'=====V-COMPONENT OF VELOCITY=====',/)
 120  FORMAT(10(E12.4))
 52   FORMAT(2X,I3,I3,5(E12.4))
      WRITE(33,123)EMAX,IMAX,JMAX      
 123  FORMAT(/'MAXIMUM ERROR=',D16.8/'ON GRID:','I=',I3,'    J=',I3,/)
      CLOSE(33)
	WRITE(13,123)EMAX,IMAX,JMAX
	CLOSE(13)
      RETURN  
      END

C*************************************************
C************SUBPROGRAM FOR RESTARTING**************
      SUBROUTINE RESTART  
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON SHI(900,700),VORTI(900,700),ERR(900,700),U(900,700),
     +V(900,700),XU(900),YV(700),X(900),Y(700),T(900,700)
      COMMON XL,YL,H,RF,RE,RI,EMAX,RA,PR
      COMMON W,YL1,KBS1,KBE1,KBHS1,KBHE1,KBHS11,KBHE11,KW,KL,KL1,
     +KBS11,KBE11
      COMMON /INDEX/M,N,K,ITER,LAST,IMAX,JMAX,INTV,KW1
C*************************************************
	OPEN(UNIT=15,FILE='PARAM.TXT')
      REWIND(15)
	DO I=1,M
	DO J=1,K
	READ(15,237)SHI(I,J),VORTI(I,J),U(I,J),V(I,J),T(I,J)
	END DO
	END DO
	CLOSE(15)

237   FORMAT(2X,5(1X,F12.8))
C*************************************************
      RETURN
      END
