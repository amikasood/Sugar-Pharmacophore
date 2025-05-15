// different places for axial and equatorial were not considered

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

//structure to store top 2 score and matching positions
//score1,position1;score2,position2
struct similar
{
  int score[200][200];
  char position[200][200][12];
};

void calculates(struct similar S, char b[12][200][20], int n, int v);

int main()
{

 int nstr,i,j,k,l;
 printf("Number of structures in input file ");
 scanf("%d", &nstr);

 char ehf[nstr][20],bb[12][200][20],temp;
 temp = ' ';
// char *temp1,*temp2;

 struct similar S1; //best score
 struct similar S2; //second best score

 for(i=0;i<200;i++){
    for(j=0;j<200;j++){
       S1.score[i][j]=0;
       S2.score[i][j]=0;
    }
 }
//read input file in EHF format. 
 FILE *ifp;
 ifp = fopen("input", "rt");
 if (ifp == NULL) {
    fprintf(stderr, "Cannot open input file\n");
    exit(1);
    }

//split each line at "_"
 i=0;
 for(i=0;i<nstr;i++){
  fgets(ehf[i], 20, ifp);
  
  for(j=0;j<20;j++){
    if(ehf[i][j]=='\n'){ehf[i][j]='\0';}
    } 
//  printf("%s\n", ehf[i]);
  j=0;
  k=0;
  l=0;
  while(ehf[i][k] != '\0'){
     if(ehf[i][k]=='_'){j=j+1;l=0;k=k+1;}
     bb[j][i][l]=ehf[i][k];
     k=k+1;
     l=l+1;
     }
//  for(j=0;j<6;j++){printf("%s ",bb[j][i]);}
//  printf("\n");
//if the exocyclic part is axial move it to third place 
  for(j=1;j<6;j++){
     if(bb[j][i][2]=='\0'){ //if there is one exocylcic group at position 1
           if(bb[j][i][1]>='A' && bb[j][i][1]<='Z'){temp=bb[j][i][1]; bb[j][i][2]=temp; bb[j][i][1]=' '; bb[j][i][3]='\0';/*printf("in loop %c\n",temp);*/}
           else if(bb[j][i][1]==' '){bb[j][i][3]='\0';bb[j][i][1]=' ';bb[j][i][2]=' ';}
           else {bb[j][i][2]=' ';bb[j][i][3]='\0';/*printf("in loop eq one %d %c\n",j,bb[j][i][0]);*/}
     }
     else if(bb[j][i][3]=='\0'){//if there is more than one exocyclic group eg. Neu5Ac
         if(bb[j][i][1]>='A' && bb[j][i][1]<='Z'){temp=bb[j][i][2]; bb[j][i][2]=bb[j][i][1];bb[j][i][1]=temp; }
     }
     }
   
//flip the ring and add it to the array
  //strcpy(bb[0][nstr+i],bb[0][i]);
  if(bb[0][i][0]=='O'){ bb[0][nstr+i][0]='o'; }
  else if(bb[0][i][0]=='o'){ bb[0][nstr+i][0]='O'; }
 
  for(j=1;j<6;j++){
     for(k=0;k<3;k++){bb[6-j][nstr+i][k]=bb[j][i][k];}
     if(bb[j][i][0]=='C'){ bb[6-j][nstr+i][0]='c'; }
     else if(bb[j][i][0]=='c'){ bb[6-j][nstr+i][0]='C'; }
//     printf("%c %c\n",(bb[j][i][0]),bb[6-j][nstr+i][0]);
     }

//copy the contents further
  for(j=6;j<12;j++){
     strcpy(bb[j][i],bb[j-6][i]);
     strcpy(bb[j][nstr+i],bb[j-6][nstr+i]);
     }
//  printf("%d %d %s %s\n",i,nstr+i,bb[10][i],bb[10][nstr+i]);
 }
 
// for(i=0;i<2*nstr;i++){for(j=0;j<6;j++){printf("%s ",bb[j][i]);}printf("\n");}
 
//calculate similarity score from top view
 calculates(S1, bb, nstr, 1);
 
 printf("\nSimilarity score with flipped structure. \n\n");
//calculate similarity score from bottom view
 calculates(S2, bb, nstr, 2);

 fclose(ifp);
 return 0;
}

//calculate similarity score
void calculates(struct similar S, char bb[12][200][20], int nstr, int view){

int i,j,k,l,m;
int ts;
char temp[12];
 //compare all of them to one another

//loop through the structures
 for(i=0;i<nstr;i++){ //first structure
    for(j=(view-1)*nstr+i+1;j<view*nstr;j++){ //second structure
       S.score[i][j] = -100;
       for(k=0;k<6;k++){ //start comparison from position k of first structure
           for(m=0;m<12;m++){temp[m]=' ';}
           ts=0;m=0;
           for(l=0;l<6;l++){

              //check all the ring atoms first
              if((bb[k+l][i][0]==bb[l][j][0])){
                 ts=ts+2.0; //same plane and same atom
 
                 //if ring atoms match, i.e. they are both carbon, check the exocyclic atoms
                 //check the equatorial place
                 if(bb[k+l][i][0]=='c' || bb[k+l][i][0]=='C'){
                    if((bb[k+l][i][1]==bb[l][j][1])){
                      ts=ts+1.0; //
                      if(bb[k+l][i][1]!=' '){temp[m]=k+l;m++;temp[m]=l;m++;}
                    }
                 //check axial
                    if((bb[k+l][i][2]==bb[l][j][2])){
                      ts=ts+1.0; //
                      if(bb[k+l][i][1]!=' '){temp[m]=k+l;m++;temp[m]=l;m++;}
                    }
                    printf("m:%d,%s ",m,temp);
                 }
              }

              //ring atoms don't match one is carbon other is oxygen
              //check if they are in the same plane
              else {
                    if((bb[k+l][i][0]=='o' && bb[l][j][0] =='c') || (bb[k+l][i][0]=='c' && bb[l][j][0] =='o')){
                      ts=ts+1.0;/*printf("%d %d %c %c ring\n",i,j,bb[l+k][i][0],bb[l][j][0]);*/} //same plane, different atom
                    else {
                         if ((bb[k+l][i][0]=='O' && bb[l][j][0]=='C') || (bb[k+l][i][0]=='C' && bb[l][j][0]=='O')){
                            ts=ts+1.0;/*printf("%d %d %c %c ring\n",i,j,bb[l+k][i][0],bb[l][j][0]);*/}}
              }
           }
           //find the highest score
           if(ts>S.score[i][j]){
              S.score[i][j]=ts;
              //for(m=0;m<12;m++){printf("%c ",temp[m]);}
              printf("\n");
           }
       }
    printf("best score for %d and %d is %d\n",i,j,S.score[i][j]);
    }
 }
}
