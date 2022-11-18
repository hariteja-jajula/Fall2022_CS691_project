#include "mpi.h"
#include <stdio.h>
#include<math.h>
#define MATSIZE 2
float ct[MATSIZE][MATSIZE] = {{1, 2}, {4, 5}};
float pt[MATSIZE][MATSIZE]= {{7,0},{17,8}};
float determinant(float at[][MATSIZE], int k)
{
  float s = 1, det = 0;
  float bt[MATSIZE][MATSIZE];
  int i, j, m, n, c;
  if (k == 1)
    {
     return (at[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                bt[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   bt[m][n] = at[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * (at[0][c] * determinant(bt, k - 1));
          s = -1 * s;
          
          }
    }
  
    return (det);
}


float transpose(float num[MATSIZE][MATSIZE], float fac[MATSIZE][MATSIZE], float r)
{
  int i, j,k;
 float b[MATSIZE][MATSIZE], inverse[MATSIZE][MATSIZE],mul[MATSIZE][MATSIZE],d;
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num, r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        inverse[i][j] = b[i][j] / d;
        }
    }
//   printf("\n\n\nThe inverse of matrix is : \n");

/*   for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         printf("\t%f", inverse[i][j]);
        }
    printf("\n");
     }
*/

for(i=0;i<MATSIZE;i++)    
{    
for(j=0;j<MATSIZE;j++)    
{    
mul[i][j]=0;    
for(k=0;k<MATSIZE;k++)    
{    
mul[i][j]+=inverse[i][k]*ct[k][j];    
}    
}    
}

int flag=0;    
for(i=0;i<MATSIZE;i++)
{
for(j=0;j<MATSIZE;j++)
{
if(mul[i][j]-pt[i][j]<=0.000001 && mul[i][j]-pt[i][j]>=-0.000001)
flag=1;
else{
flag=0;break;}
}
}

if(flag==1) 
{
for(i=0;i<MATSIZE;i++)    
{    
for(j=0;j<MATSIZE;j++)    
{    
printf("%f\t",mul[i][j]);    
}    
printf("\n");    
}

}

return 0; 
}

float cofactor(float num[MATSIZE][MATSIZE], int f)
{
 float b[MATSIZE][MATSIZE], fac[MATSIZE][MATSIZE];
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
float a;
a= transpose(num, fac, f);
return a;
}

void main(int argc, char *argv[])  {
        int numtasks, rank, next, prev, buf[2], tag1=1, tag2=2;
     
      

        MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double t1,t2;
	t1=MPI_Wtime();
	int n = pow(26,MATSIZE*MATSIZE); 
	int key[MATSIZE][MATSIZE];
	float KEY[MATSIZE][MATSIZE];
	int k = rank;
	int p;
	float ap;
	float det;
	int pt[MATSIZE][MATSIZE]; 
	while(k<n){
	p=k;
		for(int i = 0;i<MATSIZE;i++)
		{
			for(int j = 0;j<MATSIZE;j++)
			{
			key[i][j]=p%26;
			KEY[i][j]=(float)(key[i][j]);
			p=p/26;
			}	
		}
	/*calculate det of key */ 
	det = determinant(KEY,MATSIZE);
	if(det!=0)
	{
		ap=cofactor(KEY,MATSIZE);
		}	


	k=k+numtasks;	
	}


	MPI_Status status;
        t2=MPI_Wtime();
	printf("%f",t2-t1);
	MPI_Finalize(); 
}
