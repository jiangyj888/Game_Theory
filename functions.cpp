#include "functions.h"
int visited[maxNum];//通过visited数组来标记这个顶点是否被访问过，0表示未被访问，1表示被访问  
int DFS_Count;//连通部件个数，用于测试无向图是否连通，DFS_Count=1表示只有一个连通部件，所以整个无向图是连通的  
int vN;
int isInPath(vector<int> path, int g) {
	int nextGood = -1;
	for (int i = 0; i < path.size() - 1; i++) {
		if (g == path[i]) {
			return path[i + 1];
		}
	}
	return nextGood;
}
//返回新的Ai  
void SymmetricDiffrence(vector<int> path, vector<int>& Ai) {
	int temp[100];
	for (int i = 0; i < Ai.size(); i++) {
		temp[i] = -1;
	}
	for (int i = 0; i < Ai.size(); i++) {
		temp[i] = Ai[i];
	}
	//	int flag=Ai.size();
	for (int i = 0; i < Ai.size(); i++) {
		int nextGood = isInPath(path, Ai[i]);
		int flag = 0;
		if (nextGood != -1) {
			for (int j = 0; j < Ai.size(); j++) {
				if (temp[j] == Ai[i]) {
					temp[j] = -1;
					flag++;
				}
				if (temp[j] == nextGood) {
					temp[j] = -1;
					flag++;
				}
			}
			if (flag == 0) {
				for (int k = 0; k < Ai.size(); k++) {
					if (temp[k] == -1) {
						temp[k] = Ai[i];
						break;
					}
				}
				for (int k = 0; k < Ai.size(); k++) {
					if (temp[k] == -1) {
						temp[k] = nextGood;
						break;
					}

				}

			}
			if (flag == 1) {
				for (int k = 0; k < Ai.size(); k++) {
					if (temp[k] == -1) {
						temp[k] = nextGood;
						break;
					}
				}
			}
		}
	}
	vector<int> newAi;
	for (int i = 0; i < Ai.size(); i++) {
		if (temp[i] != -1) {
			newAi.push_back(temp[i]);
		}
	}
	Ai = newAi;
}
vector<int> Floyd(DirectedGraph g, vector<int> Fi, vector<int> Aj)
{
	vector<int> aimdPath;
	int A[maxNum][maxNum], path[maxNum][maxNum];
	int i, j, k, n = g.vnum;
	for (i = 1; i <= n; i++)						//给A数组置初值
		for (j = 1; j <= n; j++)
		{
			A[i][j] = g.e[i][j];
			path[i][j] = -1;
		}
	for (k = 1; k <= n; k++)						//计算Ak
	{
		for (i = 1; i <= n; i++)
			for (j = 1; j <= n; j++)
				if (A[i][j] > (A[i][k] + A[k][j]))
				{
					A[i][j] = A[i][k] + A[k][j];
					path[i][j] = k;
				}
	}
	//	printf("\n输出最短路径:\n");
		//	DisPath(A,path,n);			//输出最短路径
	int min = INF, g1, gt;
	for (int i = 0; i < Fi.size(); i++) {
		for (int j = 0; j < Aj.size(); j++) {
			if (A[Fi[i]][Aj[j]] != INF && Fi[i] != Aj[j])
			{
				if (A[Fi[i]][Aj[j]] < min) {
					min = A[Fi[i]][Aj[j]];
					g1 = Fi[i];
					gt = Aj[j];
				}
			}
		}
	}
	if (min != INF) {
		//		printf("  从%d到%d路径为:", g1, gt);
		aimdPath.push_back(g1);
		//		printf("%d,", g1);
		ppath(path, g1, gt, aimdPath);
		//		printf("%d", gt);
		aimdPath.push_back(gt);
		//		printf("\t路径长度为:%d\n", min);
	}
	return aimdPath;
}
void ppath(int path[][maxNum], int i, int j, vector<int>& aimdPath)
{
	int k;
	k = path[i][j];
	if (k == -1)  return;
	ppath(path, i, k, aimdPath);
	printf("%d,", k);
	aimdPath.push_back(k);
	ppath(path, k, j, aimdPath);
}
void dfs1(DirectedGraph* g, int i, int x, int* visited1) {
	//cout<<"顶点"<<g->v[i]<<"已经被访问"<<"\n";  
	//cout<<"顶点"<<i<<"已经被访问"<<"\n";  
	visited1[i] = 1;//标记顶点i被访问  
	//pre[i]=++point;  
	for (int j = 1; j <= g->vnum; j++)
	{
		if (g->e[i][j] == 1 && g->e[i][j] != 0 && visited1[j] == 0)
			dfs1(g, j, x, visited1);
	}
	// post[i]=++point;  

}
bool hasPath(DirectedGraph DmGraph, int x, int y) {
	int visited1[100];
	for (int i = 1; i < 100; i++) {
		visited1[i] = 0;
	}
	dfs1(&DmGraph, x, y, visited1);
	if (visited1[y] == 1) return true;
	return false;
}
void initialGraph(graph* g1, graph* g2) {
	g1->eNum = 0;
	g1->vNum = g2->vNum;
	//	cout<<g1->eNum<<","<<g1->vNum;
	for (int i = 1; i <= g1->vNum; i++) {
		for (int j = 1; j <= g1->vNum; j++) {
			if (i == j) {
				g1->e[i][j] = 0;
			}
			else {
				g1->e[i][j] = INF;
			}
		}

	}
	//	for(int i=1;i<=g1->vNum;i++){
	//		for(int j=1;j<=g1->vNum;j++){
	//			cout<<g1->e[i][j]<<" ";
	//		}
	//		cout<<"\n";
	//		
	//	}
	//	cout<<"------------------------------------";

}
void dfs(graph* g, int i)
{
	//cout<<"顶点"<<g->v[i]<<"已经被访问"<<"\n";  
	//cout<<"顶点"<<i<<"已经被访问"<<"\n";  
	visited[i] = 1;//标记顶点i被访问  

	for (int j = 1; j <= g->vNum; j++)
	{
		if (g->e[i][j] == 1 && g->e[i][j] != 0 && visited[j] == 0)
			dfs(g, j);
	}

}

void DFS(graph* g)
{
	int i;
	//初始化visited数组，表示一开始所有顶点都未被访问过  
	for (i = 1; i <= g->vNum; i++)
	{
		visited[i] = 0;

	}
	//初始化pre和post  

	//初始化连通部件数为0  
	DFS_Count = 0;
	//深度优先搜索  
	for (i = 1; i <= g->vNum; i++)
	{
		if (visited[i] == 0)//如果这个顶点为被访问过，则从i顶点出发进行深度优先遍历  
		{
			DFS_Count++;//统计调用void dfs(graph *g,int i);的次数  
			dfs(g, i);
		}
	}
}

bool isIndependent(graph* g) {
	DFS(g);
	if (g->eNum + DFS_Count > g->vNum)
		return 0;
	else
		return 1;
}
void generateEdgeSet(int m, vector<edge>* edgeSet) {
	for (int i = 1; i <= m; i++) {
		for (int j = i + 1; j <= m; j++) {
			edge edge = { i,j };
			edgeSet->push_back(edge);
			//	cout<<i<<","<<j<<"\n";
		}
	}
}
void initialDmGraph(DirectedGraph* DmGraph) {
	vector<edge> edgeSet;
	generateEdgeSet(DmGraph->vnum, &edgeSet);
	int n = DmGraph->vnum;
	for (int i = 0; i < edgeSet.size(); i++) {
		DmGraph->v[i + 1] = edgeSet[i];
		//	cout<<DmGraph->v[i+1].x<<","<<DmGraph->v[i+1].y<<"\n";
	}
	for (int i = 1; i <= n * (n - 1) / 2; i++) {
		for (int j = 1; j <= n * (n - 1) / 2; j++) {
			if (i == j) {
				DmGraph->e[i][j] = 0;
			}
			else {
				DmGraph->e[i][j] = INF;
			}
		}
		DmGraph->vnum = n * (n - 1) / 2;
		DmGraph->eNum = 0;

	}

}
void paintGraph(graph* graph, vector<edge> Ain) {
	if (Ain.size() == 0) {
		//cout<<"没有添加任何边！"<<"\n"; 
		return;
	}
	for (int i = 0; i < Ain.size(); i++) {
		graph->e[Ain[i].x][Ain[i].y] = 1;
		graph->e[Ain[i].y][Ain[i].x] = 1;
		graph->eNum++;
	}
	//		for(int i=1;i<=graph->vNum;i++){
	//			for(int j=1;j<=graph->vNum;i++){
	//				cout<<graph->e[i][j]<<" ";
	//			}
	//			cout<<"???";
	//		}
}
bool isEdgeTheSame(edge e1, edge e2) {
	return (e1.x == e2.x && e2.y == e1.y) || (e1.y == e2.x && e2.y == e1.x);
}
bool isAlreadyInGraph(edge e, vector<edge> Ai, edge edgeThatGoneOut) {

	if (isEdgeTheSame(e, edgeThatGoneOut) == true) return true;
	for (int i = 0; i < Ai.size(); i++) {
		if (isEdgeTheSame(e, Ai[i]) == true) return true;
	}
	return false;
}
bool FisAlreadyInGraph(edge e, vector<edge> Ai) {
	for (int i = 0; i < Ai.size(); i++) {
		if (isEdgeTheSame(e, Ai[i]) == true) return true;
	}
	return false;
}
void addUndirectedEdge(graph* g, edge e) {
	g->e[e.x][e.y] = 1;
	g->e[e.y][e.x] = 1;
	g->eNum = g->eNum + 1;
}
void findFi(graph graphi, vector<edge> Ai, vector<edge>& Fi) {

	graph g;
	initialGraph(&g, &graphi);
	paintGraph(&g, Ai);
	for (int j = 1; j <= g.vNum; j++) {
		for (int k = j; k <= g.vNum; k++) {
			if (graphi.e[j][k] == 1) {
				edge edge = { j,k };
				if (!FisAlreadyInGraph(edge, Ai)) {

					addUndirectedEdge(&g, edge);
					if (isIndependent(&g)) {
						Fi.push_back(edge);
					}
					deleteUndirectedEdge(&g, edge);
				}
			}

		}
	}
	//return Fi;

}
void getS(vector<vector<edge>> A, vector<edge>* s) {
	vector<edge> edgeSet;
	generateEdgeSet(vN, &edgeSet);
	vector<edge> allocatedGoods;
	int n = A.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < A[i].size(); j++) {
			allocatedGoods.push_back(A[i][j]);
			//	cout<<"allocated goods:"<<A[i][j].x<<","<<A[i][j].y<<" ";
		}
	}

	for (int i = 0; i < edgeSet.size(); i++) {
		int flag = 0;
		for (int j = 0; j < allocatedGoods.size(); j++) {
			if (isEdgeTheSame(edgeSet[i], allocatedGoods[j])) {
				flag = 1;
			}
		}
		if (flag == 0) {
			s->push_back(edgeSet[i]);
		}
	}
}
void getF(vector<graph> graphs, vector<vector<edge>> A, vector<vector<edge>>* F) {
	for (int i = 0; i < graphs.size(); i++) {
		vector<edge> Fi;
		findFi(graphs[i], A[i], Fi);
		F->push_back(Fi);
	}
}
void addDirectEdge(DirectedGraph* DmGraph, edge e1, edge e2) {
	int x1 = 0, x2 = 0;
	for (int i = 1; i <= DmGraph->vnum; i++) {
		if (isEdgeTheSame(DmGraph->v[i], e1)) {
			x1 = i;
		}

	}
	for (int j = 1; j <= DmGraph->vnum; j++) {
		if (isEdgeTheSame(DmGraph->v[j], e2)) {
			x2 = j;
		}

	}

	//	cout<<x1<<","<<x2<<"\n";
	//	cout<<"4---------------"<<"\n";
	DmGraph->e[x1][x2] = 1;
	DmGraph->eNum++;
}
int getIndex(DirectedGraph* DmGraph, edge e1) {
	int x1 = 0;
	for (int i = 1; i <= DmGraph->vnum; i++) {
		if (isEdgeTheSame(DmGraph->v[i], e1)) {
			x1 = i;
		}

	}
	return x1;
}
void deleteUndirectedEdge(graph* g, edge e) {
	g->e[e.x][e.y] = INF;
	g->e[e.y][e.x] = INF;
	g->eNum = g->eNum - 1;
}
void findDirectEdge(graph graphi, vector<edge> Ai, DirectedGraph* DmGraph) {
	vector<vector<edge>> Ain;
	vector<edge> goods;
	//vector<edge> edgeSet;
	for (int i = 0; i < Ai.size(); i++) {
		vector<edge> A;
		for (int j = 0; j < Ai.size(); j++) {
			if (i == j) {
				goods.push_back(Ai[i]);

				continue;
			}
			else {
				A.push_back(Ai[j]);
			}
		}
		Ain.push_back(A);
	}
	//	cout<<goods[0].x<<","<<goods[0].y<<"\n";
	//	cout<<Ain.size()<<"\n";
	for (int i = 0; i < Ai.size(); i++) {

		graph g1;
		initialGraph(&g1, &graphi);

		paintGraph(&g1, Ain[i]);//用去除某个边的边集作图
//		cout<<g1.vNum<<"\n";
//		for(int f=1;f<=g1.vNum;f++){
//			for(int g=1;g<=g1.vNum;g++){
//				cout<<g1.e[f][g]<<" ";
//			}
//			cout<<"\n";
//		}
		for (int j = 1; j <= g1.vNum; j++) {
			for (int k = 1; k <= g1.vNum; k++) {
				if (graphi.e[j][k] == 1) {
					edge edge = { j,k };
					//					cout<<edge.x<<","<<edge.y<<"\n";
					//					cout<<"1---------------"<<"\n";
					if (!isAlreadyInGraph(edge, Ain[i], goods[i])) {
						//						cout<<"entered:"<<edge.x<<","<<edge.y<<"\n";
						//						cout<<"2---------------"<<"\n";
						addUndirectedEdge(&g1, edge);
						if (isIndependent(&g1)) {
							//							cout<<"independent!"<<"\n";
							//							cout<<"3---------------"<<"\n";
							addDirectEdge(DmGraph, goods[i], edge);
						}
						deleteUndirectedEdge(&g1, edge);
					}
				}

			}
		}

	}
}
DirectedGraph* getDmGraph(vector<graph> graphs, vector<vector<edge>> allocation) {
	DirectedGraph* DmGraph = (DirectedGraph*)malloc(sizeof(DirectedGraph));
	DmGraph->vnum = vN;
	initialDmGraph(DmGraph);
	for (int i = 0; i < allocation.size(); i++) {
		findDirectEdge(graphs[i], allocation[i], DmGraph);
	}
	return DmGraph;
}
void printGraph(DirectedGraph DmGraph) {
	for (int i = 1; i <= DmGraph.vnum; i++) {
		for (int j = 1; j <= DmGraph.vnum; j++) {
			cout << DmGraph.e[i][j] << " ";
		}
		cout << "\n";
	}
}
void getFirstAllocation(vector<graph> graphs, vector<vector<edge>>* allocation) {
	int maxSize = graphs.size();
	unordered_set<edge, edge_Hash, edge_equal> allocations;
	for (int i = 0; i < maxSize; i++) {
		graph g = graphs[i];//取出图 
		for (int j = 1; j <= g.vNum; j++) {
			bool flag = 0;
			for (int k = j; k <= g.vNum; k++) {
				if (g.e[j][k] == 1) {
					edge e = { j,k };
					int current_size = allocations.size();
					allocations.insert(e);
					if (current_size == allocations.size()) {
						continue;
					}
					else {
						vector<edge> v;
						//					 	cout<<"被分配的边"<<e.x<<","<<e.y<<"\n"; 
						v.push_back(e);
						allocation->push_back(v);
						flag = 1;
						break;
					}
				}
			}
			if (flag == 1) {
				flag = 0;
				break;
			}


		}
	}
}
bool isSInF(edge s, vector<edge> Fi) {
	for (int i = 0; i < Fi.size(); i++) {
		//		cout<<"s:"<<s.x<<","<<s.y<<"f:"<<Fi[i].x<<","<<Fi[i].y<<"\n";
		if (isEdgeTheSame(s, Fi[i])) {
			return true;
		}
	}
	return false;
}
bool hasPathFToS(vector<edge> Fi, edge s, DirectedGraph* DmGraph) {
	for (int i = 0; i < Fi.size(); i++) {
		int x = getIndex(DmGraph, Fi[i]);
		int y = getIndex(DmGraph, s);
		if (hasPath(*DmGraph, x, y)) {
			return true;
		}
	}
	return false;
}
void printAllocations(vector<vector<edge>> allocations) {
	vector<edge> edgeSet;
	generateEdgeSet(vN, &edgeSet);
	int count = 0;
	for (int i = 0; i < allocations.size(); i++) {
		cout << "第" << i + 1 << "个玩家分得 ";
		count += allocations[i].size();
		cout << allocations[i].size() << "个物品" << "\n";
		for (int j = 0; j < allocations[i].size(); j++) {
			for (int k = 0; k < edgeSet.size(); k++) {
				if (isEdgeTheSame(allocations[i][j], edgeSet[k])) {
					cout << "物品编号: " << k + 1 << " ";
				}

			}
		}
		cout << "\n";
	}
	cout << "共有" << count << "个物品被分配出去\n";
}
bool transferGoods(vector<edge>* Aj, vector<edge>* Ak, graph graphi) {

	for (int i = 0; i < Ak->size(); i++) {
		graph graph;
		initialGraph(&graph, &graphi);
		paintGraph(&graph, *Aj);
		edge e = (*Ak)[i];
		addUndirectedEdge(&graph, e);
		if (isIndependent(&graph)) {
			Aj->push_back(e);
			Ak->erase(Ak->begin() + i);

			return true;

		}
		deleteUndirectedEdge(&graph, e);

	}
	return false;
}
void Balancedallocations(vector<vector<edge>>& allocations, graph graphi) {
	while (1) {
		int count = 0;
		for (int i = 0; i < allocations.size(); i++) {
			for (int j = 0; j < allocations.size(); j++) {
				if (i != j) {
					int rank1 = computeRank(allocations[i], graphi);
					int rank2 = computeRank(allocations[j], graphi);
					if (rank1 > rank2&& rank1 - rank2 >= 2) {
						bool isTransfered = transferGoods(&allocations[j], &allocations[i], graphi);
						if (isTransfered) count++;

					}
					else if (rank1 < rank2 && rank2 - rank1 >= 2) {
						bool isTransfered = transferGoods(&allocations[i], &allocations[j], graphi);
						if (isTransfered) count++;
					}
				}
			}
		}

		if (count == 0) {
			//		printAllocations(allocations);
			//    	cout<<"------------------\n";
			break;
		}

	}//
//	printAllocations(allocations);

}
void getParetoOptimal(vector<graph> graphs, vector<vector<edge>>* allocations) {
	getFirstAllocation(graphs, allocations);
	DirectedGraph* DmGraph = getDmGraph(graphs, *allocations);
	vector<vector<edge>> F;
	getF(graphs, *allocations, &F);
	vector<edge> s;
	getS(*allocations, &s);

	for (int i = 0; i < s.size(); i++) {
		//		cout<<"第"<<i<<"次循环"; 
		for (int j = 0; j < F.size(); j++) {
			//			cout<<"第"<<j<<"个Fi\n"<<"--------------------\n";
			int t1 = isSInF(s[i], F[j]);
			int t2 = hasPathFToS(F[j], s[i], DmGraph);
			//			cout<<"t1:"<<t1<<"t2:"<<t2<<"\n";
			if (t1 || t2) {
				((*allocations)[j]).push_back(s[i]);
				//				cout<<"第"<<i<<"个s-->"<<"第"<<j<<"个"<<"玩家\n"; 
				break;
			}
		}
	}

	//	printAllocations(*allocations);
}
int computeMMSi(graph graphi, int n) {
	vector<graph> graphs;
	for (int i = 0; i < n; i++) {
		graphs.push_back(graphi);
	}
	vector<vector<edge>> allocations;
	getParetoOptimal(graphs, &allocations);
	Balancedallocations(allocations, graphi);
	int min = INF;
	for (int i = 0; i < allocations.size(); i++) {
		if (computeRank(allocations[i], graphi) < min) {
			min = computeRank(allocations[i], graphi);
		}
	}
	return min;
}
void computeMMS(vector<graph> graphs, vector<int>& MMS) {
	for (int i = 0; i < graphs.size(); i++) {
		MMS.push_back(computeMMSi(graphs[i], graphs.size()));
	}
}
void EdgeSetToIntSet(DirectedGraph* DmGraph, vector<edge>& EdgeSet, vector<int>& IntSet) {
	for (int i = 0; i < EdgeSet.size(); i++) {
		int x = getIndex(DmGraph, EdgeSet[i]);
		IntSet.push_back(x);
	}
}

void IntSetToEdgeSet(DirectedGraph* DmGraph, vector<edge>& EdgeSet, vector<int>& IntSet) {
	vector<edge> edgeSet;
	for (int i = 0; i < IntSet.size(); i++) {
		edge e = DmGraph->v[IntSet[i]];
		edgeSet.push_back(e);
	}
	EdgeSet = edgeSet;
}
int computeRank(vector<edge> Ai, graph gi) {
	vector<edge> newAi;
	for (int i = 0; i < Ai.size(); i++) {
		if (gi.e[Ai[i].x][Ai[i].y] != INF) {
			newAi.push_back(Ai[i]);
		}
	}
	graph newGi;
	initialGraph(&newGi, &gi);
	unsigned int count = 0;
	for (int i = 0; i < newAi.size(); i++) {
		newGi.e[newAi[i].x][newAi[i].y] = 1;
		newGi.e[newAi[i].y][newAi[i].x] = 1;
		newGi.eNum++;
		int flag = isIndependent(&newGi);
		if (flag == 0) {
			newGi.e[newAi[i].x][newAi[i].y] = INF;
			newGi.e[newAi[i].y][newAi[i].x] = INF;
			newGi.eNum--;
			count++;
		}
	}
	return newAi.size() - count;

}

int computePMMS(graph graphi, vector<edge> allocationsi, vector<edge> allocationsj) {
	graph gi;
	//	vector<graph> graphs;
	initialGraph(&gi, &graphi);
	paintGraph(&gi, allocationsi);
	paintGraph(&gi, allocationsj);

	return computeMMSi(gi, 2);
}

int  computeNums(int array[maxNum][maxNum]) {
	int count = 0;
	for (int i = 0; i < vN; i++) {
		for (int j = i + 1; j < vN; j++) {
			if (array[i][j] == 1) {
				count++;
			}
		}
	}
	return count;
}
void printRawAllocations(vector<vector<edge>> allocations) {
	for (int i = 0; i < allocations.size(); i++) {
		cout << "玩家" << i + 1 << "获得的边有：\n";
		for (int j = 0; j < allocations[i].size(); j++) {
			cout << "（" << allocations[i][j].x << "," << allocations[i][j].y << ") ";
		}
		cout << endl;
	}
}
void intialGraphicMatroids3(vector<graph>& graphicMatroid) {
	int array[100][100];
	ifstream myfile("out.txt", ios::in);

	if (!myfile.is_open())
	{
		cout << "can not open this file" << endl;
		return;
	}
	int s, t;
	char buf1[1024];
	char buf2[1024];
	myfile >> buf1 >> s >> buf2 >> t;
	vN = s;
	//	cout << s << " " << t << endl;
	while (t > 0) {
		graph graphi;
		for (int i = 0; i < s; i++) {
			for (int j = 0; j < s; j++) {
				myfile >> array[i][j];
			}
		}
		int eNum = computeNums(array);
		graphi.eNum = eNum;
		graphi.vNum = vN;
		for (int k = 0; k < s; k++) {
			for (int q = 0; q < s; q++) {
				//    	cout << array[i][j] << " ";
				graphi.e[k + 1][q + 1] = array[k][q];
				//cout << array[k][q] << " ";
			}

		}
		graphicMatroid.push_back(graphi);

		t--;
	}

}