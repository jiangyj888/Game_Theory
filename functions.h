#ifndef   MY_H_FILE       //���û�ж��������  
#define   MY_H_FILE       //���������  


#include <cstdio>
#include <iostream>
#include <malloc.h> 
#include <vector>
#include <math.h>
#include <unordered_set>
#include <fstream>
using namespace std;
#define maxNum 100 //�����ڽӾ�֤����󶨵��� 
#define INF 32767				//���� ��

using namespace std;
//ͼ���ڽӾ����ʾ�ṹ  
typedef struct
{
	char v[maxNum];//ͼ�Ķ�����Ϣ  
	int e[maxNum][maxNum];//ͼ�Ķ�����Ϣ  
	int vNum;//�������  
	int eNum;//�ߵĸ���  
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
DirectedGraph* getDmGraph(vector<graph> graphs, vector<vector<edge>> allocation);//��ȡ��������ͼ 
int computeRank(vector<edge> Ai, graph gi);//����rank
void intialGraphicMatroids3(vector<graph>& graphicMatroid);//��ʼ��ͼ����
void deleteUndirectedEdge(graph* g, edge e);//ɾ������ͼ��ĳ���� 
void dfs1(DirectedGraph* g, int i, int x, int* visited1);//dfs����ͼ 
void DFS(graph* g);//����ͼ 
bool hasPath(DirectedGraph DmGraph, int x, int y);//�жϵ�(x,y)֮���Ƿ���·�� 
void initialGraph(graph* g1, graph* g2);//��ʼ������ͼ 
bool isIndependent(graph* g);//�ж�ͼ�Ƿ��л�
void generateEdgeSet(int m, vector<edge>* edgeSet);//�����߼�
void initialDmGraph(DirectedGraph* DmGraph);//��ʼ������ͼ 
void paintGraph(graph* graph, vector<edge> Ain); //ʹ�ñ߼���ͼ
bool isEdgeTheSame(edge e1, edge e2);//�ж��������Ƿ���ͬ
bool isAlreadyInGraph(edge e, vector<edge> Ai, edge edgeThatGoneOut);//�жϱ��Ƿ��Ѿ���ͼ��
bool FisAlreadyInGraph(edge e, vector<edge> Ai);//�жϱ��Ƿ��Ѿ���ͼ��
void findFi(graph graphi, vector<edge> Ai, vector<edge>& Fi);//����Fi��Fi�е�Ԫ��Ϊ{g :I+g��Ii}
void getF(vector<graph> graphs, vector<vector<edge>> A, vector<vector<edge>>* F);//�ܵ�F
void addDirectEdge(DirectedGraph* DmGraph, edge e1, edge e2);//Ϊ����ͼ��ӱ�
int getIndex(DirectedGraph* DmGraph, edge e1);// �õ�ĳ�����ڽ���ͼ�ж�Ӧ�ı��
void deleteUndirectedEdge(graph* g, edge e);//ɾ������ͼ�е�ĳһ���� 
void findDirectEdge(graph graphi, vector<edge> Ai, DirectedGraph* DmGraph);//Ϊ����ͼ��������
void getS(vector<vector<edge>> A, vector<edge>* s);//��ȡs s����δ�����䵽������е���Ʒ
void printGraph(DirectedGraph DmGraph);//����Ļ�����ͼ����Ϣ 
void getFirstAllocation(vector<graph> graphs, vector<vector<edge>>* allocation);//�õ����ʼ�ķ���
bool isSInF(edge s, vector<edge> Fi);//�ж���Ʒs���ǲ����Ѿ�������F֮��
bool hasPathFToS(vector<edge> Fi, edge s, DirectedGraph* DmGraph);//�ж�F-s֮���Ƿ����·��
void printAllocations(vector<vector<edge>> allocations);//����Ļ�����������
bool transferGoods(vector<edge>* Aj, vector<edge>* Ak, graph graphi);//����Ʒ��Ai,Aj֮�佻��
void Balancedallocations(vector<vector<edge>>& allocations, graph graphi);//ʹ��������|Ai-Aj|<=1 for all i,j
void getParetoOptimal(vector<graph> graphs, vector<vector<edge>>* allocations);//��ȡ�������������ŵĽ�
int computeMMSi(graph graphi, int n);//����ĳ�����i��MaxMin Share
void computeMMS(vector<graph> graphs, vector<int>& MMS);//����������ҵ� MaxMin Share
void MMS_Allocations(vector<graph> graphs);//��ȡMMS allocations 
void EdgeSetToIntSet(DirectedGraph* DmGraph, vector<edge>& EdgeSet, vector<int>& IntSet);//���߼�ת��Ϊ�㼯 
void IntSetToEdgeSet(DirectedGraph* DmGraph, vector<edge>& EdgeSet, vector<int>& IntSet);//���㼯ת��Ϊ�߼�
void MMS_Allocations(vector<graph> graphs);//��ȡMMS allocations 
void PMMS_Allocations(vector<graph> graphs);//��ȡPMMS allocations 
void ppath(int path[][maxNum], int i, int j, vector<int>& aimdPath);//�ݹ�floyd�㷨�е�·�� �ó����·��
vector<int> Floyd(DirectedGraph g, vector<int> Fi, vector<int> Aj);	//���������㷨��ÿ�Զ���֮������·��
int isInPath(vector<int> path, int g);//�ж�ĳ����Ʒ�Ƿ������·��֮��
void SymmetricDiffrence(vector<int> path, vector<int>& Ai); // ��� Ai��P  P={gi,gi+1 :g1��shortestPath} 
int computePMMS(graph graphi, vector<edge> allocationsi, vector<edge> allocationsj);//��������Pairwise MaxMin Share 
void printRawAllocations(vector<vector<edge>> allocations);//��ӡ�������������ı�

#endif  
