#ifndef   MY_H_FILE       //如果没有定义这个宏  
#define   MY_H_FILE       //定义这个宏  


#include <cstdio>
#include <iostream>
#include <malloc.h> 
#include <vector>
#include <math.h>
#include <unordered_set>
#include <fstream>
using namespace std;
#define maxNum 100 //定义邻接举证的最大定点数 
#define INF 32767				//定义 ∞

using namespace std;
//图的邻接矩阵表示结构  
typedef struct
{
	char v[maxNum];//图的顶点信息  
	int e[maxNum][maxNum];//图的顶点信息  
	int vNum;//顶点个数  
	int eNum;//边的个数  
}graph;
typedef struct edge
{
	int x;
	int y;
	edge(int a, int b) :x(a), y(b) {}
}edge;
typedef struct
{
	int vnum;
	int eNum;
	edge v[maxNum];
	int e[maxNum][maxNum];
}DirectedGraph;
typedef struct
{
	size_t operator()(const edge& e1) const
	{
		return hash<int>()(e1.x) ^ hash<int>()(e1.y);
	}
}edge_Hash;
typedef struct
{
	bool operator()(const edge& e1, const edge& e2) const noexcept {
		return (e1.x == e2.x && e1.y == e2.y) || (e1.y == e2.x && e1.x == e2.y);
	}
}edge_equal;
DirectedGraph* getDmGraph(vector<graph> graphs, vector<vector<edge>> allocation);//获取交换有向图 
int computeRank(vector<edge> Ai, graph gi);//计算rank
void intialGraphicMatroids3(vector<graph>& graphicMatroid);//初始化图拟阵
void deleteUndirectedEdge(graph* g, edge e);//删除无向图的某条边 
void dfs1(DirectedGraph* g, int i, int x, int* visited1);//dfs遍历图 
void DFS(graph* g);//遍历图 
bool hasPath(DirectedGraph DmGraph, int x, int y);//判断点(x,y)之间是否有路径 
void initialGraph(graph* g1, graph* g2);//初始化无向图 
bool isIndependent(graph* g);//判断图是否有环
void generateEdgeSet(int m, vector<edge>* edgeSet);//生产边集
void initialDmGraph(DirectedGraph* DmGraph);//初始化有向图 
void paintGraph(graph* graph, vector<edge> Ain); //使用边集作图
bool isEdgeTheSame(edge e1, edge e2);//判断两条边是否相同
bool isAlreadyInGraph(edge e, vector<edge> Ai, edge edgeThatGoneOut);//判断边是否已经在图中
bool FisAlreadyInGraph(edge e, vector<edge> Ai);//判断边是否已经在图中
void findFi(graph graphi, vector<edge> Ai, vector<edge>& Fi);//查找Fi，Fi中的元素为{g :I+g∈Ii}
void getF(vector<graph> graphs, vector<vector<edge>> A, vector<vector<edge>>* F);//总的F
void addDirectEdge(DirectedGraph* DmGraph, edge e1, edge e2);//为有向图添加边
int getIndex(DirectedGraph* DmGraph, edge e1);// 得到某条边在交换图中对应的编号
void deleteUndirectedEdge(graph* g, edge e);//删除无向图中的某一条边 
void findDirectEdge(graph graphi, vector<edge> Ai, DirectedGraph* DmGraph);//为交换图添加有向边
void getS(vector<vector<edge>> A, vector<edge>* s);//获取s s代表还未被分配到玩家手中的物品
void printGraph(DirectedGraph DmGraph);//在屏幕中输出图的信息 
void getFirstAllocation(vector<graph> graphs, vector<vector<edge>>* allocation);//得到最初始的分配
bool isSInF(edge s, vector<edge> Fi);//判断物品s，是不是已经存在于F之中
bool hasPathFToS(vector<edge> Fi, edge s, DirectedGraph* DmGraph);//判断F-s之间是否存在路径
void printAllocations(vector<vector<edge>> allocations);//在屏幕上输出分配结果
bool transferGoods(vector<edge>* Aj, vector<edge>* Ak, graph graphi);//将物品在Ai,Aj之间交换
void Balancedallocations(vector<vector<edge>>& allocations, graph graphi);//使分配满足|Ai-Aj|<=1 for all i,j
void getParetoOptimal(vector<graph> graphs, vector<vector<edge>>* allocations);//获取满足帕累托最优的解
int computeMMSi(graph graphi, int n);//计算某个玩家i的MaxMin Share
void computeMMS(vector<graph> graphs, vector<int>& MMS);//计算所有玩家的 MaxMin Share
void MMS_Allocations(vector<graph> graphs);//获取MMS allocations 
void EdgeSetToIntSet(DirectedGraph* DmGraph, vector<edge>& EdgeSet, vector<int>& IntSet);//将边集转换为点集 
void IntSetToEdgeSet(DirectedGraph* DmGraph, vector<edge>& EdgeSet, vector<int>& IntSet);//将点集转换为边集
void MMS_Allocations(vector<graph> graphs);//获取MMS allocations 
void PMMS_Allocations(vector<graph> graphs);//获取PMMS allocations 
void ppath(int path[][maxNum], int i, int j, vector<int>& aimdPath);//递归floyd算法中的路径 得出最短路径
vector<int> Floyd(DirectedGraph g, vector<int> Fi, vector<int> Aj);	//弗洛伊德算法从每对顶点之间的最短路径
int isInPath(vector<int> path, int g);//判断某个物品是否存在于路径之中
void SymmetricDiffrence(vector<int> path, vector<int>& Ai); // 完成 AiΔP  P={gi,gi+1 :g1∈shortestPath} 
int computePMMS(graph graphi, vector<edge> allocationsi, vector<edge> allocationsj);//计算任意Pairwise MaxMin Share 
void printRawAllocations(vector<vector<edge>> allocations);//打印分配结果所包含的边

#endif  
