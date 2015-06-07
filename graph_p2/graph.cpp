#include "graph.h"
#include "iostream"
using namespace std;
//Eric's contribution. Prints vertices in the connected component of myvertex.
void findConComp(int myvertex,int** graph,int nvertices)
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
   delete[] vfwyntctn;
   delete[] in_cc;
   return;
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
