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
#include "iostream"

// For MPI communication
#define GETR2	1
#define DIETAG	0
using namespace std;

// Used to determine MASTER or SLAVE
static int myrank;

// Function Declarations
void master();
void slave(int slavename);

//graph functions
int** allocgraph(int nvertices);
void freegraph(int**& graph,int nvertices);
int** readgraph(char* filename,int& nvertices);
void printgraph(char* filename,int** graph,int nvertices);

//This version is modified to return an array with a couple extra items--list length and original vertex.
int* findConComp(int myvertex,int** graph,int nvertices);

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
   
   // Branch off to master or slave function
   // Master has ID == 0, the slaves are then in order 1,2,3,...
   
   if(myrank==0)
   {
      master();
      printf("In main, a master process finished.");
   }
   else
   {
      slave(myrank);
      printf("In main, a slave process finished.");
   }
   
   // Finalize the MPI session
   MPI_Finalize();
   
   return(1);
}

void master()
{
   // Read in the data
   char graphfile[] = "data.txt";
   int nvertices = -1;
   int** graph = readgraph(graphfile,nvertices);
   
   int ntasks;		                        // the total number of slaves
   int jobsRunning;	                        // how many slaves we have working
   int work[1];		                        // information to send to the slaves
   MPI_Status status;	                    // MPI information

   //workresults[0] is the length of the cc and workresults[1] is the input vertex
   //The rest is the actual cc
   int* workresults = new int[nvertices+2];
   for(int i=0; i< nvertices+2; i++){workresults[i] = -1;}
   //all_results[v] stores the cc of vertex v
   int** all_results = new int*[nvertices];         
   for(int i=0; i<nvertices; i++){all_results[i] = new int[nvertices];} 
   //all_results_lens[v] stores the size of the cc of vertex v
   int* all_results_lens = new int[nvertices];
   //seen_vertex_already keeps track of the obvious to help avoid duplicate print-outs at the end.
   bool* seen_vertex_already = new bool[nvertices];
   for(int v = 0; v<nvertices; v++){seen_vertex_already[v] = 0;}


   // Find out how many slaves there are
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

   fprintf(stdout, "Total Number of processors = %d\n",ntasks);

   // Loop through the vertices and compute the connected components in parallel
   jobsRunning = 1;
   for(int v=0; v<nvertices; v++)
   {
      // This will tell the slave which variable to work on
      work[0] = v;

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
         MPI_Recv(workresults,           // where to store the results
                  nvertices+2,		     //buffer length
                  MPI_DOUBLE,	         // the type of the vector
                  MPI_ANY_SOURCE,
                  MPI_ANY_TAG, 
                  MPI_COMM_WORLD,
                  &status);              // lets us know which processor returned these results
         
         //copy the results
         int v = workresults[1];
         printf("135: Copying results on vertex [%d] \n", v);
         all_results_lens[v] = workresults[0];
         for(int i=0; i<all_results_lens[v]; i++)
         {
            all_results[v][i] = workresults[i+2]; 
         }

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

   // loop over all the slaves and collect remaining jobs
   for(int rank=1; rank<jobsRunning; rank++)
   {
      MPI_Recv(workresults,
               nvertices+2,
               MPI_DOUBLE,
               MPI_ANY_SOURCE,	// whoever is ready to report back
               MPI_ANY_TAG,
               MPI_COMM_WORLD,
               &status);
      
      //copy the results
      int v = workresults[1];
      printf("174: Copying results on vertex [%d] \n", v);
      all_results_lens[v] = workresults[0];
      for(int i=0; i<all_results_lens[v]; i++)
      {
         all_results[v][i] = workresults[i+2]; 
      }
    }

   printf("Shutting down slaves.\n");
   // Shut down the slave processes
   for(int rank=1; rank<ntasks; rank++)
   {
      printf("Master is killing slave [%d]\n",rank);
      MPI_Send(0,
               0,
               MPI_INT,
               rank,		// which node to kill
               DIETAG,		// tell it to die
               MPI_COMM_WORLD);
   }

   printf("About to print final results. \n");
   //print the unique connected components
   for(int v = 0; v<nvertices; v++)
   {
      printf("Component of vertex [%d] contains ", v);
      for(int i=0; i<all_results_lens[v]; i++)
      {
         if(seen_vertex_already[all_results[v][i]])
         {
            printf("duplicates. \n");
            break;
         }
         else
         {
            printf(" [%d] ", (int) all_results[v][i]);
            seen_vertex_already[all_results[v][i]] = 1;
         }
         printf("\n");
      }
   }
   
   //free memory
   printf("master freeing lists ca line 217\n");
   delete seen_vertex_already;
   printf("master freeing lists ca line 219\n");
   for(int v = 0; v<nvertices; v++)
   {
      delete all_results[v];
   }
   delete all_results;
   printf("master freeing graph ca line 225\n");
   freegraph(graph,nvertices);
   printf("master has finished.");
   return;
}

void slave(int slavename)
{
   // Read in the data
   char graphfile[] = "data.txt";
   int nvertices = -1;
   int** graph = readgraph(graphfile,nvertices);
   int* workresults;
   int work[1]; 		               // the inputs from the master
   MPI_Status status;		           // for MPI communication

   // the slave listens for instructions...
   int notDone = 1;
   while(notDone)
   {
      printf("Slave %d is waiting\n",slavename);
      MPI_Recv(&work,		     // the inputs from the master
               1,		         // the size of the inputs
               MPI_INT,		     // the type of the inputs
               0,		         // from the MASTER node (rank=0)
               MPI_ANY_TAG,	     // any type of order is fine
               MPI_COMM_WORLD,
               &status);
      printf("Slave %d just received smth\n",slavename);

      // switch on the type of work request
      switch(status.MPI_TAG)
      {
         case GETR2:
            return;
            // Get conn component
            printf("Slave %d has received vertex [%d]\n", slavename,work[0]);
            workresults = findConComp(work[0], graph, nvertices);

            // Send the results
            MPI_Send(workresults,
                     nvertices+2,
                     MPI_DOUBLE,
                     0,		// send it to the master
                     0,		// doesn't need a TAG
                     MPI_COMM_WORLD);

            printf("Slave %d finished processing work request [%d]\n",
                   slavename,work[0]);

            break;

         case DIETAG:
            freegraph(graph,nvertices);
            delete workresults
            printf("Slave %d was told to die\n",slavename);
            return;

         default:
            notDone = 0;
            printf("The slave code should never get here.\n");
            return;
      }
   }
   freegraph(graph,nvertices);
   return;
}

//workresults[0] is the length of the cc and workresults[1] is the input vertex
//The rest is the actual cc
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
   //For all vertices whose neighbors need to be checked
   for(int i=0; i<vfwyntctn_size;i++)
   {
      //"vwnabc" is "vertex whose neighbors are being checked"
      //It is a vertex already known to be in the CC of myvertex.
      vwnabc = vfwyntctn[i];

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
   
   int* workresults = new int[vfwyntctn_size+2];
   workresults[0] = vfwyntctn_size;
   workresults[1] = myvertex;
   for(int i=0; i<vfwyntctn_size; i++)
   {
      workresults[i+2] = vfwyntctn[i];
   }
   delete[] vfwyntctn;
   delete[] in_cc;
   return workresults;
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
