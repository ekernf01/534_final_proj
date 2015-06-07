/*
 This program computes the R2 of several regressions in parallel.
 Compile the program using the makefile provided.
 
 Run the program using the command:

 mpirun -np 10 parallel_graph 
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <mpi.h>
#include <iomanip>

// For MPI communication
#define GETR2	1
#define DIETAG	0

// Used to determine MASTER or SLAVE
static int myrank;

// Global variables
int nobservations = 40;
int nvariables = 1000;
double* Y = NULL;
double** X = NULL;

double ssy = 0.0;	// used in R2 calculation

// Function Declarations
void NormStand();
void master();
void slave(int slavename);

//graph functions
int** allocgraph(int nvertices);
void freegraph(int**& graph,int nvertices);
int** readgraph(char* filename,int& nvertices);
void printgraph(char* filename,int** graph,int nvertices);
int main(int argc, char* argv[])
{
   ///////////////////////////
   // START THE MPI SESSION //
   ///////////////////////////
   MPI_Init(&argc, &argv);

   /////////////////////////////////////
   // What is the ID for the process? //   
   /////////////////////////////////////
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   
   //find the connected component of vertex i+1
   findConComp(i,graph,nvertices);
   
   // Branch off to master or slave function
   // Master has ID == 0, the slaves are then in order 1,2,3,...
   
   if(myrank==0)
   {
      master();
   }
   else
   {
      slave(myrank);
   }
   
   //free memory
   freegraph(graph,nvertices);
   
   // Finalize the MPI session
   MPI_Finalize();
   
   return(1);
}

void master()
{
   int ntasks;		// the total number of slaves
   int jobsRunning;	// how many slaves we have working
   int work[1];		// information to send to the slaves
   int workresults[2]; // info received from the slaves
   MPI_Status status;	// MPI information

   // Read in the data
   char graphfile[] = "data.txt";
   int nvertices = -1;
   int** graph = readgraph(graphfile,nvertices);
   
   // Find out how many slaves there are
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

   fprintf(stdout, "Total Number of processors = %d\n",ntasks);

   // Loop through the vertices and compute the connected components in
   // parallel
   jobsRunning = 1;

   for(my_vertex=0; my_vertex<nvertices; my_vertex++)
   {
      // This will tell the slave which variable to work on
      work[0] = my_vertex;

      if(jobsRunning < ntasks) // Do we have an available processor?
      {
         // Send out a work request
         MPI_Send(&work, 	     // the vector with the vertex
                  1,             // the size of the vector
                  MPI_INT,       // the type of the vector
                  jobsRunning,	 // the ID of the slave to use
                  GETR2,	     // tells the slave to work, not die
                  MPI_COMM_WORLD); // send the request out to anyone who is available
         
         printf("Master sends out work request [%d] to slave [%d]\n", work[0],jobsRunning);

         // Increase the # of processors in use
         jobsRunning++;

      }
      else // all the processors are in use! Wait, receive, send.
      {
         MPI_Recv(workresults,	         // where to store the results
                  nvertices,		     // the size of the vector
                  MPI_DOUBLE,	         // the type of the vector
                  MPI_ANY_SOURCE,
                  MPI_ANY_TAG, 
                  MPI_COMM_WORLD,
                  &status);              // lets us know which processor returned these results

         printf("Master has received the result of work request [%d] from slave [%d]\n",
                (int) workresults[0],status.MPI_SOURCE);
 
         // Print out the results
         fprintf(fout, "%d\t%f\n", (int)workresults[0]+1, workresults[1]);

         printf("Master sends out work request [%d] to slave [%d]\n",
                work[0],status.MPI_SOURCE);

         // Send out a new work order to the processors that just returned
         MPI_Send(&work,
                  1,
                  MPI_INT,
                  status.MPI_SOURCE, // the slave that just returned
                  GETR2,
                  MPI_COMM_WORLD); 
      } // end of else
   } // this is the loop over vertices


   ///////////////////////////////////////////////////////////////
   // NOTE: we still have some work requests out that need to be
   // collected. Collect those results now!
   ///////////////////////////////////////////////////////////////

   // loop over all the slaves
   for(rank=1; rank<jobsRunning; rank++)
   {
      MPI_Recv(workresults,
               2,
               MPI_DOUBLE,
               MPI_ANY_SOURCE,	// whoever is ready to report back
               MPI_ANY_TAG,
               MPI_COMM_WORLD,
               &status);

       printf("Master has received the result of work request [%d]\n", (int) workresults[0]);
 
      //save the results received
      fprintf(fout, "%d\t%f\n", (int)workresults[0]+1, workresults[1]);
   }

   printf("Tell the slave to die\n");

   // Shut down the slave processes
   for(rank=1; rank<ntasks; rank++)
   {
      printf("Master is killing slave [%d]\n",rank);
      MPI_Send(0,
               0,
               MPI_INT,
               rank,		// which node to kill
               DIETAG,		// tell it to die
               MPI_COMM_WORLD);
   }

   printf("got to the end of Master code\n");
   
   //print the unique connected components
   eiuevgeoahv
   
   //free memory
   freegraph(graph,nvertices);

   // return to the main function
   return;
}

void slave(int slavename)
{
   int work[1];			               // the inputs from the master
   double workresults[nvertices];	   // the outputs for the master
   MPI_Status status;		           // for MPI communication

   // the slave listens for instructions...
   int notDone = 1;
   while(notDone)
   {
      printf("Slave %d is waiting\n",slavename);
      MPI_Recv(&work,		// the inputs from the master
               1,		       // the size of the inputs
               MPI_INT,		   // the type of the inputs
               0,		   // from the MASTER node (rank=0)
               MPI_ANY_TAG,	// any type of order is fine
               MPI_COMM_WORLD,
               &status);
      printf("Slave %d just received smth\n",slavename);

      // switch on the type of work request
      switch(status.MPI_TAG)
      {
         case GETR2:
            // Get conn component

            printf("Slave %d has received work request [%d]\n", slavename,work[0]);
            workresults = avleuabe; 

            // Send the results
            MPI_Send(&workresults,
                     nvertices,
                     MPI_DOUBLE,
                     0,		// send it to the master
                     0,		// doesn't need a TAG
                     MPI_COMM_WORLD);

            printf("Slave %d finished processing work request [%d]\n",
                   slavename,work[0]);

            break;

         case DIETAG:
            printf("Slave %d was told to die\n",slavename);
            return;

         default:
            notDone = 0;
            printf("The slave code should never get here.\n");
            return;
      }
   }

   // No memory to clean up, so just return to the main function
   return;
}

// Data must have zero mean and unit variance
// This is only for 1 variable regressions without an intercept
double GetR2(int v)
{
   int	i;
   double tmp;

   tmp = 0.0;
   for(i=0; i<nobservations; i++)
   {
      tmp += (X[i][v]*Y[i]);
   }
   tmp = tmp*tmp / ((double)nobservations - 1.0);

   return(tmp / ssy);
}


void NormStand()
{
   int i, j;
   double tmp;

   for(j=0; j<nvariables; j++)
   {
      tmp = 0.0;
      for(i=0; i<nobservations; i++)
      {
         tmp += X[i][j];
      }
      for(i=0; i<nobservations; i++)
      {
         X[i][j] -= tmp/((double)nobservations);
      }
   }

   tmp = 0.0;
   for(i=0; i<nobservations; i++)
   {
      tmp += Y[i];
   }
   for(i=0; i<nobservations; i++) Y[i] -= tmp/((double)nobservations);

   // Now make the data have unit sample variance
   for(j=0; j<nvariables; j++)
   {
      tmp = 0.0;
      for(i=0; i<nobservations; i++)
      {
         tmp += X[i][j]*X[i][j];
      }

      tmp = sqrt(tmp / ((double)nobservations-1.0));

      for(i=0; i<nobservations; i++)
      {
         X[i][j] = X[i][j]/tmp;
      }
   }   

   // Do the same for Y
   tmp = 0.0;
   for(i=0; i<nobservations; i++)
   {
      tmp += Y[i]*Y[i];
   }
   tmp = sqrt( tmp / ((double)nobservations - 1.0));

   for(i=0; i<nobservations; i++)
   {
      Y[i] = Y[i]/tmp;
   }

   return;
}

#include "graph.h"
#include "iostream"
using namespace std;
//Eric's contribution. Prints vertices in the connected component of myvertex.
int* findConComp(int myvertex,int** graph,int nvertices)
{
   //vfwyntctn means "vertices for which you need to check the neighbors."
   //Actually, part of this array will have already had its neighbors checked.
   //And by the end of the loop, it will have a simpler meaning:
   // it will list the cc of myvertex.
   int* vfwyntctn = new int[nvertices];
   int vfwyntctn_size = 1;
   vfwyntctn[0] = myvertex;
   
   //in_cc is a vector whose [v] element will say whether v is in myvertex's connected component.
   bool* in_cc = new bool[nvertices];
   //Initialize to have zeros except at myvertex.
   for(int v=0; v<nvertices;v++)
   {
      in_cc[v]=0;
   }
   in_cc[myvertex] = 1;
   
   int vwnabc = 0;
   cout << endl << "The vertex " << myvertex << " (indexing from 0) has this connected component :" << endl;
   //For all vertices whose neighbors need to be checked
   for(int i=0; i<vfwyntctn_size;i++)
   {
      //"vwnabc" is "vertex whose neighbors are being checked"
      //It is a vertex already known to be in the CC of myvertex.
      vwnabc = vfwyntctn[i];
      cout << vwnabc << "  ";

      //Check the potential neighbors.
      for(int maybe_neighbor=0; maybe_neighbor<nvertices; maybe_neighbor++)
      {
         //If you find one that isn't listed in the cc yet,
         // but that is adjacent to vwnabc, then add it to cc,
         // and make sure to check its neighbors, i.e. put it in vfwyntctn
         if (graph[vwnabc][maybe_neighbor]&&!in_cc[maybe_neighbor])
         {
            vfwyntctn[vfwyntctn_size] = maybe_neighbor;
            in_cc[maybe_neighbor] = 1;
            vfwyntctn_size++;
         }
      }
   }
   cout << endl ;
   return in_cc;
}


//the graph is stored as an incidence matrix
//An incidence matrix is a symmetric square matrix
//with dimension equal to the number of vertices in the graph
//The (i,j) entry in this matrix takes value one
//if there is an edge between vertices i and j,
// and takes value zero otherwise.
int** allocgraph(int nvertices)
{
   int i;
   int** graph = new int*[nvertices];

   for(i=0;i<nvertices;i++)
   {
      graph[i] = new int[nvertices];
      //initially the graph has no edges
      memset(graph[i],0,nvertices*sizeof(int));
   }

   return(graph);
}

//free the memory for an incidence matrix
void freegraph(int**& graph,int nvertices)
{
   int i;

   for(i=0;i<nvertices;i++)
   {
      delete[] graph[i];
      graph[i] = NULL;
   }
   delete[] graph;
   graph = NULL;

   return;
}

//reads a graph from an input file
//The format of the input file is:
//number of vertices, followed by each edge
//It returns the graph as an incidence matrix
int** readgraph(char* filename,int& nvertices)
{
   FILE* in = fopen(filename,"r");
   if(NULL==in)
   {
      fprintf(stderr,"Cannot open input file [%s]\n",filename);
      exit(1);
   }

   int i,j;
   fscanf(in,"%d",&i);
   nvertices = i;

   int** graph = allocgraph(nvertices);

   while(2==fscanf(in,"%d %d",&i,&j))
   {
      graph[i-1][j-1] = graph[j-1][i-1] = 1;
   }

   fclose(in);

   return(graph);
}

//prints out a graph in Graphviz format
void printgraph(char* filename,int** graph,int nvertices)
{
   int i,j;

   FILE* out = fopen(filename,"w");
   if(NULL==out)
   {
      fprintf(stderr,"Cannot open output file [%s]\n",filename);
      exit(1);
   }

   fprintf(out,"graph G {\n");
   for(i=0;i<nvertices-1;i++)
   {
      for(j=i+1;j<nvertices;j++)
      {
	 if(graph[i][j])
	 {
	    fprintf(out,"%d -- %d;\n",i+1,j+1);
         }
      }
   }
 
   fprintf(out,"}\n");

   return;
}
