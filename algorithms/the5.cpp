#include "the5.h"
bool breaker=true;
/* 
    in the5.h "struct Room" is defined as below:
    
    struct Room {
        int id;
        char content;
        vector<Room*> neighbors;
    };

*/
void maze_visit(Room* room, int i, vector<int>* path, unsigned* color, int* p,vector<Room*> maze){
    if(breaker){
        color[i]=1;
        if(room->content=='*'){
                    breaker=false;
                    path->push_back(p[i]);
                    color[i]=2;
        }else{
            
            for(int j=0;j<room->neighbors.size();j++){
                if(color[room->neighbors[j]->id]==0){
                    path->push_back(room->neighbors[j]->id);
                    if(room->neighbors[j]->content=='*')
                        for(int k=0;k<maze.size();k++){
                            color[k]=2;
                        }
                    if(p[room->neighbors[j]->id]==-1)
                        p[room->neighbors[j]->id]=i;
                    maze_visit(room->neighbors[j],room->neighbors[j]->id,path,color,p,maze);
                }
            }
            if(p[i]!=-1)
                path->push_back(p[i]);
            color[i]=2;
        }
    }
}


vector<int> maze_trace(vector<Room*> maze) { 
    vector<int> path;
    
    // 0 white
    // 1 gray
    // 2 black
    unsigned color[maze.size()];
    
    int p[maze.size()];
    
    for(int i=0;i<maze.size();i++){
        color[i]=0;
        p[i]=-1;
    }
    
    path.push_back(0);
    maze_visit(maze[0],0,&path,color,p,maze);
    
    

    return path; // this is a dummy return value. YOU SHOULD CHANGE THIS!
}



