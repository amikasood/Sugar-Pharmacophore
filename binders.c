//different places for axial and exuatorial 
//[1] for equatorial and [2] for axial. small case for eq, capital for ax

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define USAGE "\n\
\tbinders.o input_file pharmacophore \n"

//structure to store the monosaccharide shape
struct shape
{
   char *bb[12];
   char *eq[12];
   char *ax[12];
};

int main(int argc,char *argv[])
{

   int i, j, k, l, nstr, pharm_size;
   char ch[50], *temp[3], *pharm;
   if(argc<3){printf("\nError:     Insufficient number of arguments on command line.\n %s \n",USAGE);exit(1);}
   if(argc>3){printf("\nError:     Too many arguments on command line.\n %s \n",USAGE);exit(1);}

   //read file with EHF representation of monosaccharides
   FILE *ifp;
   if(argc==3){
      ifp = fopen(argv[1], "rt");
      if (ifp == NULL) {
         fprintf(stderr, "Cannot open input file: %s\n",argv[1]);
         exit(1);
      }
      pharm=argv[2]; 
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

   char input[nstr][50], *aa[nstr][6], backup[nstr][50];

   ifp = fopen(argv[1], "rt");
   i=0;
   while(fgets(input[i], 50, ifp)!=NULL){
      if(input[i][0]!='\0' ){i++;}
   }
   nstr=i;

   if(nstr>1){printf("\n\tPharmacophore will searched using %d structures.\n\n",nstr);}
   if(nstr==1){printf("\n\tThere is only one structure in the input file.\n\n");}
   if(nstr==0){printf("\n\tThere are no structure in the input file.\n\n");exit(1);}

   struct shape monosacc[2*nstr];
   struct shape pharmacophore;

   for(i=0;i<2*nstr;i++){
      for(j=0;j<12;j++){
         monosacc[i].bb[j]=NULL;
         monosacc[i].eq[j]=NULL;
         monosacc[i].ax[j]=NULL;
      }
      strcpy(backup[i],input[i]);
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
   
//divide the pharmacophore
   j=0;
   aa[0][j] = strtok(pharm,"_");
   while(aa[0][j] != NULL){
      j=j+1;
      aa[0][j] = strtok(NULL,"_");
   }
   if(j>6){printf("\nError:     More than 6 atoms in the ring in structure %d. \n",i+1);exit(1);}
   strtok(aa[0][j-1], "\n");
   pharm_size=j;
         //divide the string to find the exocyclic groups
   for(j=0;j<pharm_size;j++){
      k=0;
      temp[0]=temp[1]=temp[2]=NULL;
      temp[k] = strtok(aa[0][j],"-");
      while(temp[k] != NULL){
         k=k+1;
         temp[k] = strtok(NULL,"-");
      }
      pharmacophore.bb[j]=temp[0];
      if(pharmacophore.bb[j][0]=='c' || pharmacophore.bb[j][0]=='C' || pharmacophore.bb[j][0]=='o' || pharmacophore.bb[j][0]=='O'){;}
      else{printf("\nError:     Please check the ring atoms of the pharmacophore %c \n",pharmacophore.bb[j][0]);exit(1);}
      if(temp[1]!=NULL && temp[1][0]>='a' && temp[1][0] <='z') {
         pharmacophore.ax[j]=temp[1];
      }
      else if (temp[1]!=NULL){
         pharmacophore.eq[j]=temp[1];
      }
      if(temp[2]!=NULL && temp[2][0]>='a' && temp[2][0] <='z') {
         if(pharmacophore.ax[j]!=NULL){printf("\nError:     More than one axial group in the pharmacophore. \n");exit(1);}
         pharmacophore.ax[j]=temp[2];
      }
      else if (temp[2]!=NULL){
         if(pharmacophore.eq[j]!=NULL){printf("\nError:     More than one equatorial group in the pharmacophore. \n");exit(1);}
         pharmacophore.eq[j]=temp[2];
      }
      if(pharmacophore.eq[j]==NULL){pharmacophore.eq[j]=" ";}
      if(pharmacophore.ax[j]==NULL){pharmacophore.ax[j]=" ";}
   }
//   for(j=0;j<pharm_size;j++){printf("%s-%s-%s_",pharmacophore.bb[j],pharmacophore.ax[j],pharmacophore.eq[j]);}
//   printf("\n");

//Find pharmacophore in the structures in the input file 
   FILE *sc = fopen("score.out", "w");
   int match;

   for(i=0;i<nstr;i++){ //first structure
      for(k=0;k<pharm_size;k++){
         match=1;
         for(l=0;l<pharm_size;l++){
            //check all the ring atoms first
            if(strcmp(pharmacophore.bb[l]," ")!=0){
               if(strcmp(monosacc[i].bb[k+l],pharmacophore.bb[l])==0){
                  if(strcmp(pharmacophore.ax[l]," ")!=0){
                     if(strcmp(monosacc[i].ax[k+l],pharmacophore.ax[l])==0){match=1;}else{match=0;}
                  } //if axial groups are same
                  if(strcmp(pharmacophore.eq[l]," ")!=0){
                     if(strcmp(monosacc[i].eq[k+l],pharmacophore.ax[l])==0){match=1;}else{match=0;}
                  } //if equatorial groups are same
               }
               else{match=0;}
            }
            if(match==0)break;
         }//loop l
         if(match==1)break;
      }//loop k
      if(match==1){printf("%d\t%s",i+1,backup[i]);}
      
//find pharmacophore in flipped structure
      if(match==0){
      for(k=0;k<pharm_size;k++){
         match=1;
         for(l=0;l<pharm_size;l++){
         //check all the ring atoms first
            if(strcmp(pharmacophore.bb[l]," ")!=0){
               if(strcmp(monosacc[i+nstr].bb[k+l],pharmacophore.bb[l])==0){
                  if(strcmp(pharmacophore.ax[l]," ")!=0){
                     if(strcmp(monosacc[i+nstr].ax[k+l],pharmacophore.ax[l])==0){match=1;}else{match=0;}
                  } //if axial groups are same
                  if(strcmp(pharmacophore.eq[l]," ")!=0){
                     if(strcmp(monosacc[i+nstr].eq[k+l],pharmacophore.eq[l])==0){match=1;}else{match=0;}
                  } //if equatorial groups are same
               }
               else{match=0;}
            }
            if(match==0)break;
         }//loop l
         if(match==1)break;
      }//loop k
//      fprintf(sc,"\n");
//          printf("%d ",i);
//          printf("\n");
      if(match==1){
         printf("%d\t%s",i+1,backup[i]);
      }
      }//if match=0
   }//loop i
   
   fclose(sc);
   fclose(ifp);
   return 0;
}
