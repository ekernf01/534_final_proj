/*
 This program sorts the elements of a given vector.
*/

#include "vectors.h"
#include "tree.h"
#include <stdio.h>
#include <stdlib.h>
using namespace std;

int main(){
    int i;
    int n = 8;
    
    char InputFile[] = "somenumbers.txt";
    char OutputFile[] = "mybinarytree.dot";

    //initialise a vector of length n
    double* vector = allocvector(n);
    
    readvector(n,vector,InputFile);
    
    printf("Original vector\n");
    printvector(vector,n);
    printf("\n\n");

    //create a binary tree
    LPNode mytree = MakeNewNode(vector[0]);

    //now add all the other elements of the vector
    for(i=1;i<n;i++)
    {
        treeInsert(mytree,vector[i]);
    }

    //sort in increasing order by traversing the tree in inorder
    double* vectorIncreasing = allocvector(n);
        
    //remove all the elements in order, deleting the tree in the process
    for(i=1;i<n;i++)
    {
        double temp_lowest = 0;
        mytree = remove_lowest(mytree, &temp_lowest);
        vectorIncreasing[i] = temp_lowest;
    }
    
    printf("Vector sorted in increasing order\n");
    printvector(vectorIncreasing,n);
    fprintvector(vectorIncreasing,n, OutputFile);
    printf("\n\n");
    
    //free memory
    freevector(vector);
    freevector(vectorIncreasing);

    return(1);
}
