#include "the2.h"
#include <cmath>
#include <iostream>
using namespace std;
// You may write your own helper functions here
int offset_val(int x, int offset, int groupSize){
    int scale=(int)pow(10,groupSize);
    int shift=(int)pow(10,offset);
    return (x/shift)%scale;
}
long countingSort(int* arr, bool ascending, int arraySize, int groupSize, int offset){
    
    long numberOfIterations = 0;
    int *output=(int*)malloc(sizeof(int)*arraySize);
    int k=(int)pow(10,groupSize);
    int count_arr[k];
    for(int i=0;i<k;i++){
        count_arr[i]=0;
        
    }
    
    for(int i=0;i<arraySize;i++){
        count_arr[offset_val(arr[i],offset,groupSize)]++;
        numberOfIterations++;
        
    }
    
    for(int i=1;i<k;i++){
        count_arr[i]+=count_arr[i-1];
        numberOfIterations++;
    }
    
    if(ascending){
        for(int i=arraySize-1;i>=0;i--){
            output[count_arr[offset_val(arr[i],offset,groupSize)]-1] = arr[i];
            count_arr[offset_val(arr[i],offset,groupSize)]--;
            numberOfIterations++;
        }    
        for(int i=0;i<arraySize;i++){
            arr[i]=output[i];
            numberOfIterations++;
        }
    }
    else{
        for(int i=0;i<arraySize;i++){
            output[count_arr[offset_val(arr[i],offset,groupSize)]-1] = arr[i];
            count_arr[offset_val(arr[i],offset,groupSize)]--;
            numberOfIterations++;
        }
        for(int i=0;i<arraySize;i++){
            arr[i]=output[arraySize-i-1];
            numberOfIterations++;
        }
        for(int i=0;i<arraySize;i++){
            cout<<arr[i]<<" ";
        }
        cout<<"\n";
    }
    free(output);
    return numberOfIterations;
}

long multiDigitRadixSort(int* arr, bool ascending, int arraySize, int groupSize, int maxDigitLength){
    long numberOfIterations = 0;

    for(int i=0;i<maxDigitLength;i+=groupSize)
        if((i+groupSize)<maxDigitLength)
            numberOfIterations+=countingSort(arr,ascending,arraySize,groupSize,i);
        else
            numberOfIterations+=countingSort(arr,ascending,arraySize,maxDigitLength-i,i);
            
    return numberOfIterations;
}