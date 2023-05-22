#include "the1.h"
#include <climits>
#include <iostream>
using namespace std;
//You can add your own helper functions
void insertion_sort(int* arr, int size, long& comparison, long& swap){
    int i=1,j,x;
    while(i<size){
        
        x=arr[i];
        j=i-1;
        comparison++;
        while(j>=0 && arr[j]>x){
            swap++;
            comparison++;
            arr[j+1]=arr[j];
            j--;
            
        }
        
        arr[j+1]=x;
        i++;
        
    }
}

void heapify(int* arr, int* current_part, int i, int size, long& comparison, long& swap){
        cout<<i<<" ";
        if(i>=size)
            return;
        if(2*i<size){
            comparison++;
            if(arr[2*i]<arr[i]){
                swap++;
                int temp=arr[2*i];
                arr[2*i]=arr[i];
                arr[i]=temp;
                temp=current_part[2*i];
                current_part[2*i]=current_part[i];
                current_part[i]=temp;
                
                
            }
        }
        if(2*i-1<size){
            comparison++;
            if(arr[2*i-1]<arr[i]){
                swap++;
                int temp=arr[2*i-1];
                arr[2*i-1]=arr[i];
                arr[i]=temp;
                temp=current_part[2*i-1];
                current_part[2*i-1]=current_part[i];
                current_part[i]=temp;
                heapify(arr,current_part,i+1,size,comparison,swap);
            }
        }
}


int kWayMergeSortWithHeap(int* arr, int K, int size, long& comparison, long& swap){
    
  int number_of_calls =1;
	
	if(size<K){
    
        insertion_sort(arr,size,comparison,swap);
        
    }else{
    
        int part_size=size/K;
        int heap[K];
        int current_part[K];
        int cursors[K];
        
        
        for(int i=0;i<K;i++)
            number_of_calls += kWayMergeSortWithHeap((arr+part_size*i), K, part_size, comparison,swap);
        
        insertion_sort(arr,size,comparison,swap);
        
    }
	return number_of_calls;
}

