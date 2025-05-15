//different places for axial and exuatorial 
//[1] for equatorial and [2] for axial. small case for eq, capital for ax

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define USAGE "\n\
\tpredict.o input \n"

//structure to store the monosaccharide shape
struct shape
{
   char *bb[12];
   char *eq[12];
   char *ax[12];
};

struct unique_score
{
   int bb[6];
   int ax[6];
   int eq[6];
};

int main(int argc,char *argv[])
{

   int i, j, k, l, ts, plane, nstr;
   char ch[50], *temp[3];
   if(argc<2){printf("\nError:     Insufficient number of arguments on command line.\n %s \n",USAGE);exit(1);}
   if(argc>2){printf("\nError:     Too many arguments on command line.\n %s \n",USAGE);exit(1);}

   //read file with EHF representation of monosaccharides
   FILE *ifp;
   if(argc==2){
      ifp = fopen(argv[1], "rt");
      if (ifp == NULL) {
         fprintf(stderr, "Cannot open input file: %s\n",argv[1]);
         exit(1);
      }
   }

   nstr=0;
   while(!feof(ifp))
   {
     if(fgets(ch,50,ifp)!= NULL)
     {
       nstr++;
     }
   }
   fclose(ifp);

   char input[nstr][50], *aa[nstr][6];

   ifp = fopen(argv[1], "rt");
   i=0;
   while(fgets(input[i], 50, ifp)!=NULL){
      if(input[i][0]!='\0' ){i++;}
   }
   nstr=i;

   if(nstr>1){printf("\n\tPharmacophore will be predicted using %d structures.\n\tIf possible, the structures should be listed in the descending order of affinity.\n\n",nstr);}
   if(nstr==1){printf("\n\tThere is only one structure in the input file.\n\n");exit(1);}
   if(nstr==0){printf("\n\tThere are no structure in the input file.\n\n");exit(1);}

   struct shape monosacc[2*nstr];
   struct shape aligned[nstr];
   struct shape unique[nstr];
   struct shape pharmacophore;

   for(i=0;i<2*nstr;i++){
      for(j=0;j<12;j++){
         monosacc[i].bb[j]=NULL;
         monosacc[i].eq[j]=NULL;
         monosacc[i].ax[j]=NULL;
      }
   }

   for(i=0;i<nstr;i++){
      for(j=0;j<12;j++){
         aligned[i].bb[j]=NULL;
         aligned[i].eq[j]=NULL;
         aligned[i].ax[j]=NULL;
         unique[i].bb[j]=NULL;
         unique[i].eq[j]=NULL;
         unique[i].ax[j]=NULL;
      }
   }

   for(j=0;j<12;j++){
      pharmacophore.bb[j]=NULL;
      pharmacophore.eq[j]=NULL;
      pharmacophore.ax[j]=NULL;
   }

   //divide the input strings and store them in struct array
   for(i=0;i<nstr;i++){
      j=0;
      aa[i][j] = strtok(input[i],"_");
      while(aa[i][j] != NULL){
         j=j+1;
         aa[i][j] = strtok(NULL,"_");
      }
      if(j>6){printf("\nError:     More than 6 atoms in the ring in structure %d. \n",i+1);exit(1);}
      if(j<6){printf("\nError:     Less than 6 atoms in the ring in structure %d. \n",i+1);exit(1);}
      strtok(aa[i][5], "\n");

      //divide the string to find the exocyclic groups
      for(j=0;j<6;j++){
         k=0;
         temp[0]=temp[1]=temp[2]=NULL;
         temp[k] = strtok(aa[i][j],"-");
         while(temp[k] != NULL){
            k=k+1;
            temp[k] = strtok(NULL,"-");
         }
         monosacc[i].bb[j]=temp[0];

         if(monosacc[i].bb[j][0]=='c' || monosacc[i].bb[j][0]=='C' || monosacc[i].bb[j][0]=='o' || monosacc[i].bb[j][0]=='O'){;}
         else{printf("\nError:     Please check the ring atoms for structure %d, %c \n",i+1,monosacc[i].bb[j][0]);exit(1);}

         if(temp[1]!=NULL && temp[1][0]>='a' && temp[1][0] <='z') {
            monosacc[i].ax[j]=temp[1];
         }
         else if (temp[1]!=NULL){
            monosacc[i].eq[j]=temp[1];
         }
         if(temp[2]!=NULL && temp[2][0]>='a' && temp[2][0] <='z') {
            if(monosacc[i].ax[j]!=NULL){printf("\nError:     More than one axial group in structure %d. \n",i+1);exit(1);}
            monosacc[i].ax[j]=temp[2];
         }
         else if (temp[2]!=NULL){
            if(monosacc[i].eq[j]!=NULL){printf("\nError:     More than one equatorial group in structure %d. \n",i+1);exit(1);}
            monosacc[i].eq[j]=temp[2];
         }
         if(monosacc[i].eq[j]==NULL){monosacc[i].eq[j]=" ";}
         if(monosacc[i].ax[j]==NULL){monosacc[i].ax[j]=" ";}
      }
   
//flip the ring and add it to the array beggining nstr+1 to nstr+nstr
//array bb[j][i] from i=0 to nstr-1, input format; i=nstr to 2*nstr-1 flipped
      monosacc[nstr+i].bb[0]=monosacc[i].bb[0];
      if(strcmp(monosacc[i].bb[0],"O")==0){ monosacc[nstr+i].bb[0]="o"; }
      else if(strcmp(monosacc[i].bb[0],"o")==0){ monosacc[nstr+i].bb[0]="O"; }
      monosacc[i].ax[0]=monosacc[i].eq[0]=" ";
      monosacc[i].ax[6]=monosacc[i].eq[6]=" ";
      monosacc[nstr+i].ax[0]=monosacc[nstr+i].eq[0]=" ";
      monosacc[nstr+i].ax[6]=monosacc[nstr+i].eq[6]=" ";

      for(j=1;j<6;j++){
         monosacc[nstr+i].ax[6-j]=monosacc[i].ax[j];
         monosacc[nstr+i].eq[6-j]=monosacc[i].eq[j];
         monosacc[nstr+i].bb[6-j]=monosacc[i].bb[j];
         if(strcmp(monosacc[i].bb[j],"C")==0){ monosacc[nstr+i].bb[6-j]="c"; }
         else if(strcmp(monosacc[i].bb[j],"c")==0){ monosacc[nstr+i].bb[6-j]="C"; }
      }
//copy the contents further
      for(j=6;j<12;j++){
         monosacc[i].bb[j]=monosacc[i].bb[j-6];
         monosacc[nstr+i].bb[j]=monosacc[nstr+i].bb[j-6];
         monosacc[i].ax[j]=monosacc[i].ax[j-6];
         monosacc[nstr+i].ax[j]=monosacc[nstr+i].ax[j-6];
         monosacc[i].eq[j]=monosacc[i].eq[j-6];
         monosacc[nstr+i].eq[j]=monosacc[nstr+i].eq[j-6];
      }
//          printf("%d ",i);
//          for(j=0;j<12;j++){printf("%s-%s-%s_",monosacc[i].bb[j],monosacc[i].ax[j],monosacc[i].eq[j]);}
//          printf("\n");
//          printf("%d ",i+nstr);
//          for(j=0;j<12;j++){printf("%s-%s-%s_",monosacc[nstr+i].bb[j],monosacc[nstr+i].ax[j],monosacc[nstr+i].eq[j]);}
//          printf("\n");
   }//loop i

//Calculate scores 
   int score;
   FILE *sc = fopen("score.out", "w");

   for(l=0;l<6;l++){
      aligned[0].bb[l]=monosacc[0].bb[l];
      aligned[0].ax[l]=monosacc[0].ax[l];
      aligned[0].eq[l]=monosacc[0].eq[l];
   }

   fprintf(sc,"Score\tAligned structures\n");
   for(i=1;i<nstr;i++){ //first structure
      //for(j=0;j<i+1;j++){fprintf(sc,"\t");fprintf(ps,"\t");}

      score=-100;

      for(k=0;k<6;k++){
         ts=0;
         for(l=0;l<6;l++){
            //check all the ring atoms first
            if(strcmp(monosacc[i].bb[k+l],monosacc[0].bb[l])==0){
               ts=ts+2.0; //same plane and same atom
               if(monosacc[i].bb[k+l][0]=='C' || monosacc[i].bb[k+l][0]=='c'){
                  if(strcmp(monosacc[i].ax[k+l],monosacc[0].ax[l])==0){
                     ts=ts+1;
                  } //if axial groups are same
                  if(strcmp(monosacc[i].eq[k+l],monosacc[0].eq[l])==0){
                     ts=ts+1;
                  } //if equatorial groups are same
               }
               else {
                 if(monosacc[i].bb[k+l][0]=='O' || monosacc[i].bb[k+l][0]=='o'){
                    ts=ts+2.0;
                 }
               }
            }
            // if ring atoms don't match, check if they are in the same place
            else {
               plane=0;
               if((strcmp(monosacc[i].bb[k+l],"o")==0 && strcmp(monosacc[0].bb[l],"c")==0) || (strcmp(monosacc[i].bb[k+l],"c")==0 && strcmp(monosacc[0].bb[l],"o")==0)){
                  ts=ts+1.0;plane=1;} //same plane, different atom
               if ((strcmp(monosacc[i].bb[k+l],"O")==0 && strcmp(monosacc[0].bb[l],"C")==0) || (strcmp(monosacc[i].bb[k+l],"C")==0 && strcmp(monosacc[0].bb[l],"O")==0)){
                  ts=ts+1.0;plane=1;}
               if(plane==1){
                  if(strcmp(monosacc[i].ax[k+l],monosacc[0].ax[l])==0){
                     ts=ts+1;}
                  if(strcmp(monosacc[i].eq[k+l],monosacc[0].eq[l])==0){
                     ts=ts+1;}
               }
            }
         }//loop l

         //find the highest score
         if(ts>score){
            score=ts;
            for(l=0;l<6;l++){
               aligned[i].bb[l]=monosacc[i].bb[k+l];
               aligned[i].ax[l]=monosacc[i].ax[k+l];
               aligned[i].eq[l]=monosacc[i].eq[k+l];
            }
         }
      }//loop k
       //score flipped structure
      for(k=0;k<6;k++){
         ts=0;
         for(l=0;l<6;l++){
         //check all the ring atoms first
            if(strcmp(monosacc[i+nstr].bb[k+l],monosacc[0].bb[l])==0){
               ts=ts+2.0; //same plane and same atom
               if(monosacc[i+nstr].bb[k+l][0]=='C' || monosacc[i+nstr].bb[k+l][0]=='c'){ //if its a c ring atom, check exocyclic atoms
                  if(strcmp(monosacc[i+nstr].ax[k+l],monosacc[0].ax[l])==0){
                     ts=ts+1;
                  } //if axial groups are same
                  if(strcmp(monosacc[i+nstr].eq[k+l],monosacc[0].eq[l])==0){
                     ts=ts+1;
                  } //if equatorial groups are same
               }
               else {
                 if(monosacc[i+nstr].bb[k+l][0]=='O' || monosacc[i+nstr].bb[k+l][0]=='o'){
                    ts=ts+2.0;
                 }
               }
            }
            // if ring atoms don't match, check if they are in the same place
            else {
               plane=0;
               if((strcmp(monosacc[0].bb[l],"o")==0 && strcmp(monosacc[i+nstr].bb[k+l],"c")==0) || (strcmp(monosacc[0].bb[l],"c")==0 && strcmp(monosacc[i+nstr].bb[k+l],"o")==0)){
                  ts=ts+1.0;plane=1;} //same plane, different atom
               if((strcmp(monosacc[0].bb[l],"O")==0 && strcmp(monosacc[i+nstr].bb[k+l],"C")==0) || (strcmp(monosacc[0].bb[l],"C")==0 && strcmp(monosacc[i+nstr].bb[k+l],"O")==0)){
                  ts=ts+1.0;plane=1;}//same plane, different atom
               if(plane==1){
                  if(strcmp(monosacc[i+nstr].ax[k+l],monosacc[0].ax[l])==0){
                     ts=ts+1;}
                  if(strcmp(monosacc[i+nstr].eq[k+l],monosacc[0].eq[l])==0){
                     ts=ts+1;}
               }
            }
         }//loop l
         //find the highest score
         if(ts>score){
            score=ts;
            for(l=0;l<6;l++){
               aligned[i].bb[l]=monosacc[i+nstr].bb[k+l];
               aligned[i].ax[l]=monosacc[i+nstr].ax[k+l];
               aligned[i].eq[l]=monosacc[i+nstr].eq[k+l];
            }
         }
      }//loop k
      fprintf(sc,"%d\t",score);
      for(j=0;j<5;j++){
         fprintf(sc,"%s",aligned[i].bb[j]);
         if(strcmp(aligned[i].eq[j]," ")!=0){fprintf(sc,"-%s",aligned[i].eq[j]);}
         if(strcmp(aligned[i].ax[j]," ")!=0){fprintf(sc,"-%s",aligned[i].ax[j]);}
         fprintf(sc,"_");
      }
      fprintf(sc,"%s",aligned[i].bb[5]);
      if(strcmp(aligned[i].eq[5]," ")!=0){fprintf(sc,"-%s",aligned[i].eq[5]);}
      if(strcmp(aligned[i].ax[5]," ")!=0){fprintf(sc,"-%s",aligned[i].ax[5]);}
      fprintf(sc,"\n");
   }//loop i
   fprintf(sc,"\n\n");

//calculate pharmacophore
   struct unique_score us[nstr], n, pharm_score;
   for(l=0;l<6;l++){
      for(i=0;i<nstr;i++){
      us[i].bb[l]=0;
      us[i].eq[l]=0;
      us[i].ax[l]=0;
      }
   } 
   for(l=0;l<6;l++){
     unique[0].bb[l]=aligned[0].bb[l];
     unique[0].ax[l]=aligned[0].ax[l]; 
     unique[0].eq[l]=aligned[0].eq[l];
     us[0].bb[l]=1;
     us[0].ax[l]=1;
     us[0].eq[l]=1;
     n.bb[l]=0;
     n.ax[l]=0;
     n.eq[l]=0;
     pharm_score.bb[l]=0;
     pharm_score.ax[l]=0;
     pharm_score.eq[l]=0;
   }

   for(l=0;l<6;l++){
      for(i=1;i<nstr;i++){
         for(j=0;j<n.bb[l]+1;j++){
            if(strcmp(aligned[i].bb[l],unique[j].bb[l])==0){
               us[j].bb[l]++;
               //printf("test: %d %d %d %s %s \n",l,i,j,aligned[i].bb[l],unique[j].bb[l]);
               break;
            }
            else{
               if(j==n.bb[l]){
                  n.bb[l]++;
                  unique[n.bb[l]].bb[l]=aligned[i].bb[l];
                  us[n.bb[l]].bb[l]++;
                  break;
               }
            }
         }//loop j backbone
         for(j=0;j<n.ax[l]+1;j++){
            if(strcmp(aligned[i].ax[l],unique[j].ax[l])==0){
               us[j].ax[l]++;
               //printf("test: %d %d %d %s %s \n",l,i,j,aligned[i].ax[l],unique[j].ax[l]);
               break;
            }
            else{
               if(j==n.ax[l]){
                  n.ax[l]++;
                  unique[n.ax[l]].ax[l]=aligned[i].ax[l];
                  us[n.ax[l]].ax[l]++;
                  break;
               }
            }
         }//loop j axial
         for(j=0;j<n.eq[l]+1;j++){
            if(strcmp(aligned[i].eq[l],unique[j].eq[l])==0){
               us[j].eq[l]++;
               //printf("test: %d %d %d %s %s \n",l,i,j,aligned[i].eq[l],unique[j].eq[l]);
               break;
            }
            else{
               if(j==n.eq[l]){
                  n.eq[l]++;
                  unique[n.eq[l]].eq[l]=aligned[i].eq[l];
                  us[n.eq[l]].eq[l]++;
                  break;
               }
            }
         }//loop j equatorial
      }
   }
   for(l=0;l<6;l++){
      for(j=0;j<n.bb[l]+1;j++){
         printf("unique bb[%d]:%s ;score:%d \n",l,unique[j].bb[l],us[j].bb[l]);
         if(j==0){pharmacophore.bb[l]=unique[j].bb[l];pharm_score.bb[l]=us[j].bb[l];}
         if(j>0){
            if(us[j].bb[l]>us[j-1].bb[l]){
               pharmacophore.bb[l]=unique[j].bb[l];
               pharm_score.bb[l]=us[j].bb[l];
            }
         }
      }//loop for j backbone
      for(j=0;j<n.ax[l]+1;j++){
         printf("unique ax[%d]:%s ;score:%d \n",l,unique[j].ax[l],us[j].ax[l]);
         if(j==0){pharmacophore.ax[l]=unique[j].ax[l];pharm_score.ax[l]=us[j].ax[l];}
         if(j>0){
            if(us[j].ax[l]>us[j-1].ax[l]){
               pharmacophore.ax[l]=unique[j].ax[l];
               pharm_score.ax[l]=us[j].ax[l];
            }
         }
      }//loop for j axial
      for(j=0;j<n.eq[l]+1;j++){
         printf("unique eq[%d]:%s ;score:%d \n",l,unique[j].eq[l],us[j].eq[l]);
         if(j==0){pharmacophore.eq[l]=unique[j].eq[l];pharm_score.eq[l]=us[j].eq[l];}
         if(j>0){
            if(us[j].eq[l]>us[j-1].eq[l]){
               pharmacophore.eq[l]=unique[j].eq[l];
               pharm_score.eq[l]=us[j].eq[l];
            }
         }
      }//loop for j equatorial
   }

   FILE *ehf = fopen("EHF_inp.plot", "w");
   float s,s1,s2,s3;
   for(l=0;l<6;l++)
   {
      s=pharm_score.bb[l];
      s=s/nstr;
      fprintf(ehf,"P%d = %f\n",l+1,s);
      fprintf(ehf,"P%dN = \"%s\"\n",l+1,pharmacophore.bb[l]);
      if(strcmp(pharmacophore.ax[l]," ")!=0) {
         s=pharm_score.ax[l];
         s=s/nstr;
         fprintf(ehf,"A%d = %f\n",l+1,s);
         fprintf(ehf,"A%dN = \"%s\"\n",l+1,pharmacophore.ax[l]);
      }
      if(strcmp(pharmacophore.eq[l]," ")!=0) {
         s=pharm_score.eq[l];
         s=s/nstr;
         fprintf(ehf,"E%d = %f\n",l+1,s);
         fprintf(ehf,"E%dN = \"%s\"\n",l+1,pharmacophore.eq[l]);
      }
   }
   fclose(ehf);
   system("gnuplot EHF.plot");

   fprintf(sc,"Pharmacophore\n");
   for(l=0;l<5;l++){
   if(strcmp(pharmacophore.bb[l]," ")!=0){
      if(strcmp(pharmacophore.ax[l]," ")!=0){
         if(strcmp(pharmacophore.eq[l]," ")!=0)
            {s1=pharm_score.bb[l]/nstr;
             s2=pharm_score.eq[l]/nstr;
             s3=pharm_score.ax[l]/nstr;
             fprintf(sc,"%s-%s-%s_",pharmacophore.bb[l],pharmacophore.eq[l],pharmacophore.ax[l]);}
         else{s1=pharm_score.bb[l]/nstr;
             s2=pharm_score.ax[l]/nstr;
             fprintf(sc,"%s-%s_",pharmacophore.bb[l],pharmacophore.ax[l]);}
      }
      else{
         if(strcmp(pharmacophore.eq[l]," ")!=0){fprintf(sc,"%s-%s_",pharmacophore.bb[l],pharmacophore.eq[l]);}
         else{fprintf(sc,"%s_",pharmacophore.bb[l]);}
      }
   }
   else{fprintf(sc,"%s_",pharmacophore.bb[l]);}
   }
      if(strcmp(pharmacophore.bb[5]," ")!=0){
      if(strcmp(pharmacophore.ax[5]," ")!=0){
         if(strcmp(pharmacophore.eq[5]," ")!=0)
            {fprintf(sc,"%s-%s-%s",pharmacophore.bb[5],pharmacophore.eq[5],pharmacophore.ax[5]);}
         else{fprintf(sc,"%s-%s",pharmacophore.bb[5],pharmacophore.ax[5]);}
      }
      else{
         if(strcmp(pharmacophore.eq[5]," ")!=0){fprintf(sc,"%s-%s",pharmacophore.bb[5],pharmacophore.eq[5]);}
         else{fprintf(sc,"%s",pharmacophore.bb[5]);}
      }
      }
      else{fprintf(sc,"%s",pharmacophore.bb[5]);}
   fprintf(sc,"\n");
  
   fclose(sc);
   fclose(ifp);
   return 0;
}
