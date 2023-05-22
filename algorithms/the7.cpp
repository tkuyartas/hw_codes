#include "the7.h"
/*struct Road {
    std::pair<int, int> buildings;
    int time;
    Road(std::pair<int, int> bs, int t) : buildings(bs), time(t) {}
};*/
// Implement your solution after this line
void fw_sp(std::vector<std::vector<int>>& adj_m , std::vector<std::vector<int>>& path_m, int n) {
    for(int k=0;k<n;k++){
    
        for(int i=0;i<n;i++){
            
            for(int j=0;j<n;j++){
                
                if(adj_m[i][k]==INT_MAX || adj_m[k][j]==INT_MAX)
                    continue;
                if(adj_m[i][j]>adj_m[i][k]+adj_m[k][j]){
                    adj_m[i][j]=adj_m[i][k]+adj_m[k][j];
                    path_m[i][j]=path_m[i][k];
                }
            }
        }
    }

}
std::vector<int> make_path(std::vector<std::vector<int>> path_m,int s,int d){
    std::vector<int> path;
    path.push_back(s);
    int start=s;
    int end=d; 
    while(start!=end){
        start=path_m[start][end];
        path.push_back(start);
    }
    return path;
}
void add(std::vector<int>& a,std::vector<int> b){
    for(unsigned i=1;i<b.size();i++)
        a.push_back(b[i]);
}
void CanCatch(int n, std::vector<Road> roads, int s, int d, int x, int y, int l, int printPath) {
    std::vector<int> path;
    std::vector<std::vector<int>> adj_m( n , std::vector<int> (n, INT_MAX));
    std::vector<std::vector<int>> path_m( n , std::vector<int> (n, -1));
    for(int i=0;i<n;i++){
        adj_m[i][i]=0;
        path_m[i][i]=i;
    }
    for(unsigned i=0;i<roads.size();i++){
        path_m[roads[i].buildings.first][roads[i].buildings.second]=roads[i].buildings.second;
        path_m[roads[i].buildings.second][roads[i].buildings.first]=roads[i].buildings.first;
        
        adj_m[roads[i].buildings.first][roads[i].buildings.second]=roads[i].time;
        adj_m[roads[i].buildings.second][roads[i].buildings.first]=roads[i].time;
        
    }
    int min;
    int imp_flag=1;
    fw_sp(adj_m,path_m,n);
    
    int s_d=adj_m[s][d];
    int s_x=adj_m[s][x];
    int s_y=adj_m[s][y];
    int x_y=adj_m[x][y];
    int y_x=adj_m[y][x];
    int y_d=adj_m[y][d];
    int x_d=adj_m[x][d];
    
    if(s_y+y_x+x_d<=l && s_y+y_x+x_d<s_x+x_y+y_d){ //s-y-x-d
    
        min=s_y+y_x+x_d;
        path=make_path(path_m,s,y);
        add(path,make_path(path_m,y,x));
        add(path,make_path(path_m,x,d));
        std::cout<<"YES BOTH "<<min<<"MINUTES\n";
        
        
    }else if(s_x+x_y+y_d<=l){    //s-x-y-d
        min=s_x+x_y+y_d;
        path=make_path(path_m,s,x);
        add(path,make_path(path_m,x,y));
        add(path,make_path(path_m,y,d));
        std::cout<<"YES BOTH "<<min<<"MINUTES\n";
    }else if(s_y+y_d<=l && s_x+x_d>s_y+y_d){    //s-y-d
        min=s_y+y_d;
        path=make_path(path_m,s,y);
        add(path,make_path(path_m,y,d));
        std::cout<<"YES DORM "<<min<<"MINUTES\n";
    }else if(s_x+x_d<=l){    //s-x-d
        min=s_x+x_d;
        path=make_path(path_m,s,x);
        add(path,make_path(path_m,x,d));
        std::cout<<"YES PRINTER "<<min<<"MINUTES\n";
    }else if(s_d<=l){        //s-d
        min=s_d;
        path=make_path(path_m,s,d);
        std::cout<<"YES DIRECTLY "<<min<<"MINUTES\n";
        
    }else{
        imp_flag=0;
        std::cout<<"IMPOSSIBLE\n";
    }
    if(printPath && imp_flag)
        PrintRange(path.begin(), path.end());
    
   
    
    // You can use the PrintRange function to print the path.
    // It helps you print elements of containers with iterators (e.g., std::vector).
    // Usage: PrintRange(path.begin(), path.end());
    
}


