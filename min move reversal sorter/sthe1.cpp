#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <tuple>
using std::vector;
using std::tuple;
using namespace std;

int heuristic(const vector<int> &input,int n)
{
    int t=0;
    int count=0;

    for(int i=1;i<n;i++){
        if(input[i] != t+1)
            count++;
        t=input[i];
    }
    if(t+1!=n)
        count++;
    return ceil(float(count)/3);
}

void shift(vector<int> &input,int i, int j,int s){
    
    rotate(input.begin()+i,input.begin()+i+s,input.begin()+j+1);
    
}
bool sorted(vector<int> &input,int n){
    for(int i=1;i<n;i++)
        if(input[i]!=i)
            return false;
    return true;
}
bool dfs(vector<tuple<int,int,int>> &chain, int limit, vector<int> &input,int n){
    if (sorted(input,n))
    {
        return true;
    }
    else if (int(chain.size()) + heuristic(input,n) > limit)
    {
        return false;
    }
    else
    {
        int i=1;
        while(i<n){
            if(input[i]!=i)
                break;
            i++;
        }
        //first unordered index
        
        int last=n;
        while(last>i){
            if(input[last]!=last)
                break;
            last--;
        }
        //last unordered index
        vector<int> temp;
        tuple<int,int,int> newElement;
        for (; i <n-1; i++)
        {
            
            for (int j=i+1; j <last; j++)
            {
                for(int s=1;s<=j-i;s++){
                    get<0>(newElement)=i;
                    get<1>(newElement)=j;
                    get<2>(newElement)=s;
                    temp=input;
                        
                    chain.push_back(newElement);
                    shift(temp,i,j,s);
                    
                    
                    
                    if (dfs(chain, limit, temp,n))
                    {
                        return true;
                    }
                    chain.pop_back();
                }
            }
            
        }
        return false;
    }
}
vector<tuple<int,int,int>> findShiftOrder(vector<int> &input)
{
    vector<tuple<int,int,int>> chain = {};
    int n=input.size();
    for(int limit=heuristic(input,n); dfs(chain, limit, input,n)==false; limit++);

    return chain;
}


int main(void){
    vector<int> input;
    int x;
    input.push_back(0);
    while(cin>>x){
        if(x=='\0')
            break;
        input.push_back(x);
    }
    
    vector<tuple<int,int,int>> chain=findShiftOrder(input);
    
    x=chain.size();
    
    cout<<x<<"\n";
    
    for(int i=0;i<x;i++){
        cout<<get<0>(chain[i])<<" "<<get<1>(chain[i])<<" "<<get<2>(chain[i])<<"\n";
    }
    

    
}
