#include <stdio.h>
#include <math.h>

#define  NMAX  10


     double C[NMAX][NMAX], XM[NMAX][NMAX];
     double R[NMAX][NMAX];
     double DAO[NMAX], RAD[NMAX];

     int P[5][2*NMAX+5];

     int ID,IDES,IF1,IOPTIMAL,ISOU,IX,IY,LT,NR;
     double CC,CT,QT,TT,XINFCC;


void Init() {
  int I,J;
  printf("\n TRANSPORT MODEL\n\n");
  printf(" NUMBER OF SOURCES ? "); scanf("%d", &ISOU);
  printf("\n NUMBER OF DESTINATIONS ? "); scanf("%d", &IDES);
  printf("\n INPUT THE AVAILABLE SOURCE QUANTITIES:\n");
  for (I=1; I<=ISOU; I++) {
    printf("  SOURCE #%d ? ", I); scanf("%lf", &DAO[I]);
  }
  printf("\n INPUT THE REQUIRED DESTINATION QUANTITIES:\n");
  for (I=1; I<=IDES; I++) {
    printf("  DESTINATION #%d ? ", I); scanf("%lf", &RAD[I]);
  }
  printf("\n INPUT TRANSPORT COSTS MATRIX:\n");
  for (I=1; I<=ISOU; I++)
    for (J=1; J<=IDES; J++) {
      printf("  FROM SOURCE #%d TO DESTINATION #%d ? ", I, J);
      scanf("%lf", &C[I][J]);
    }
  printf("\n\n TRANSPORT MODEL\n");
}

  void Corner();
  void Optimal();
  void TCost();
  void Seek_Path(int I,int J);
  void Increase(int I,int J);

void SubMain() {  //STEPPING STONE ALGORITHM
  Corner();
  Optimal();
  TCost();
}

void Corner() {  //N-W CORNER
//Labels: e20, e50, e60
    int I,J;
    I=1; J=1;
e20:if (DAO[I]<=RAD[J]) goto e50;
    XM[I][J] += RAD[J]; DAO[I] -= RAD[J];
    RAD[J]=0; J++;
    goto e60;  //else
e50:XM[I][J] += DAO[I]; RAD[J] -= DAO[I];
    DAO[I]=0; I++;
e60:if (I<=ISOU && J<=IDES) goto e20;
}

void Optimal() {
//Labels: e10, e70, e140, e150
     int I,J;
e10: XINFCC=0.0;
     for (I=1; I<=ISOU; I++)
       for (J=1; J<=IDES; J++) {
         if (XM[I][J] != 0.0) goto e70;
         Seek_Path(I,J);
         Increase(I,J);
e70:;  }
     if (XINFCC>=0.0) {
       IOPTIMAL=1;
       goto e150;
     }
     for (I=1; I<=LT; I++) {
       IX=P[3][I]; IY=P[4][I];
       if (I % 2 == 0) {
         XM[IX][IY] -= TT;
         goto e140;
       }
       XM[IX][IY] += TT;
e140:;}
e150:if (IOPTIMAL==0) goto e10;
}

void Seek_Path(int I, int J) {
//Labels: e70, e160, e260
     int I1,I2;
     for (I1=1; I1<=ISOU; I1++)
       for (I2=1; I2<=IDES; I2++)
         R[I1][I2]=XM[I1][I2];
     for (I1=1; I1<=ISOU; I1++) R[I1][0]=0.0;
     for (I2=1; I2<=IDES; I2++) R[0][I2]=0.0;
     R[I][J]=1.0;
e70: for (I2=1; I2<=IDES; I2++) {
       if (R[0][I2]==1.0) goto e160;
       NR=0;
       for (I1=1; I1<=ISOU; I1++)
         if (R[I1][I2] != 0.0) NR++;
       if (NR!=1) goto e160;
       for (I1=1; I1<=ISOU; I1++) R[I1][I2]=0.0;
       R[0][I2]=1.0; IF1=1;
e160:;}
     for (I1=1; I1<=ISOU; I1++) {
       if (R[I1][0]==1.0) goto e260;
       NR=0;
       for (I2=1; I2<=IDES; I2++)
         if (R[I1][I2] != 0.0) NR++;
       if (NR!=1) goto e260;
       for (I2=1; I2<=IDES; I2++) R[I1][I2]=0.0;
       R[I1][0]=1.0; IF1=1;
e260:;}
     if (IF1==1) {
       IF1=0; goto e70;
     }
}

void Increase(int I, int J) {
//Labels: e20,e70,e130,e170,e180,e230
     int I1,I2;
     P[1][1]=I; P[2][1]=J; IX=I; IY=J; ID=1; CC=0.0; QT=999999.0;
e20: ID++; IF1=0;
     for (I1=1; I1<=ISOU; I1++) {
       if (R[I1][IY]==0.0 || I1==IX) goto e70;
       P[1][ID]=I1; P[2][ID]=IY; IX=I1; CC -= C[IX][IY];
       IF1=1; I1=ISOU;
       if (XM[IX][IY] < QT && XM[IX][IY] > 0.0)  QT=XM[IX][IY];
e70:;}
     if (IF1==0) goto e170;
     ID++; IF1=0;
     for (I2=1; I2<=IDES; I2++) {
       if (R[IX][I2]==0.0 || I2==IY) goto e130;
       P[1][ID]=IX; P[2][ID]=I2; IY=I2; CC += C[IX][IY];
       IF1=1; I2=IDES;
e130:;}
     if (IF1==0) goto e170;
     if (IX!=I || IY!=J) goto e20;
     goto e180;
e170:printf(" DEGENERATE SOLUTION !\n");
     return;
e180:if (CC>0.0 || CC>XINFCC) goto e230;
     TT=QT; XINFCC=CC; ID--; LT=ID;
     for (I1=1; I1<=ID; I1++) {
       P[3][I1]=P[1][I1]; P[4][I1]=P[2][I1];
     }
e230:;}

void TCost() {
    int I,J;
    CT=0.0;
    printf("\n TRANSPORTS:\n");
    for (I=1; I<=ISOU; I++)
      for (J=1; J<=IDES; J++) {
        CT += XM[I][J]*C[I][J];
        if (XM[I][J]==0.0) goto e10;
        printf("    FROM SOURCE #%d TO DESTINATION #%d: %8.2f\n", I, J, XM[I][J]);
e10:;}
    printf("\n TOTAL TRANSPORT COST: %10.1f\n\n", CT);
}


int main()  {
  Init();
  SubMain();
  return 0;
}

