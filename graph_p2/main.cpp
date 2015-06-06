#include "graph.h"

int main()
{
   char graphfile[] = "data.txt";
   char dotfile[] = "data.dot";

   int nvertices = -1;
   int** graph = readgraph(graphfile,nvertices);

   //print out the graph in Graphviz format
   printgraph(dotfile,graph,nvertices);

   //free memory
   freegraph(graph,nvertices);

   return(1);
}
