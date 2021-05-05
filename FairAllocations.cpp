// FairAllocations.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//


#include "functions.h"
using namespace std;
extern int vN;
void MMS_Allocations(vector<graph> graphs, vector<vector<edge>>& allocations) {
	//	vector<vector<edge>> allocations;
	getParetoOptimal(graphs, &allocations);
	//	cout<<"paretoOptimal:\n";
	//	printAllocations(allocations);
	vector<int> MMS;
	computeMMS(graphs, MMS);
	cout << "各个玩家的maxmin share的值为:\n";
	for (int i = 0; i < MMS.size(); i++) {
		cout << "玩家" << i + 1 << ":" << MMS[i] << " ";
	}
	cout << "\n";
	cout << endl;
	//	cout << "----------------------------------------------------------------------------------------";
	DirectedGraph* DmGraph;
	DmGraph = getDmGraph(graphs, allocations);
	//printGraph(*DmGraph);
	vector<int> Sminus;
	vector<int> Splus;
	for (int i = 0; i < allocations.size(); i++) {
		if (computeRank(allocations[i], graphs[i]) < MMS[i]) {
			Sminus.push_back(i);
		}
		else if (computeRank(allocations[i], graphs[i]) > MMS[i]) {
			Splus.push_back(i);

		}
	}
	while (Sminus.size() > 0) {
		for (int i = 0; i < Sminus.size(); i++) {
			//count++;
			vector<edge> Fi;
			findFi(graphs[Sminus[i]], allocations[Sminus[i]], Fi);
			vector<int> newFi;
			EdgeSetToIntSet(DmGraph, Fi, newFi);
			int index;
			if (Splus.size() != 0) {
				index = Splus[0];
			}
			Splus.erase(Splus.begin());
			vector<int> newAj;
			EdgeSetToIntSet(DmGraph, allocations[index], newAj);
			vector<int> shortestPath = Floyd(*DmGraph, newFi, newAj);
			for (int j = 0; j < allocations.size(); j++) {
				if (j != Sminus[i] && j != index) {
					vector<int> IntSet;
					EdgeSetToIntSet(DmGraph, allocations[j], IntSet);
					SymmetricDiffrence(shortestPath, IntSet);
					IntSetToEdgeSet(DmGraph, allocations[j], IntSet);
				}
			}

			vector<int> IntSet1;
			EdgeSetToIntSet(DmGraph, allocations[Sminus[i]], IntSet1);
			SymmetricDiffrence(shortestPath, IntSet1);
			IntSetToEdgeSet(DmGraph, allocations[Sminus[i]], IntSet1);
			int gt = shortestPath[shortestPath.size() - 1];
			edge et = DmGraph->v[gt];
			int g1 = shortestPath[0];
			edge e1 = DmGraph->v[g1];
			allocations[Sminus[i]].push_back(e1);
			int ap;
			for (ap = 0; ap < allocations[index].size(); ap++) {
				if (isEdgeTheSame(et, allocations[index][ap])) {
					break;
				}
			}
			allocations[index].erase(allocations[index].begin() + ap);
			if (computeRank(allocations[Sminus[i]], graphs[Sminus[i]]) == MMS[Sminus[i]]) {
				Sminus.erase(Sminus.begin() + i);
			}
			if (computeRank(allocations[index], graphs[index]) > MMS[index]) {
				Splus.push_back(index);

			}
			DmGraph = getDmGraph(graphs, allocations);
			//printGraph(*DmGraph);
		}

	}



}
void PMMS_Allocations(vector<graph> graphs, vector<vector<edge>>& allocations) {

	getParetoOptimal(graphs, &allocations);
	int n = graphs.size();
	//	int count = 10;
	while (true) {
		int flag = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					int ranki = computeRank(allocations[i], graphs[i]);
					int PMMSij = computePMMS(graphs[i], allocations[i], allocations[j]);
					if (ranki < PMMSij) {
						flag = 1;
						vector<edge> Fi;
						findFi(graphs[i], allocations[i], Fi);
						//if (Fi.size() == 0) flag = 0;
						for (int m = 0; m < Fi.size(); m++) {
							int xxx = 0;
							for (int n = 0; n < allocations[j].size(); n++) {
								if (isEdgeTheSame(Fi[m], allocations[j][n])) {
									allocations[i].push_back(allocations[j][n]);
									allocations[j].erase(allocations[j].begin() + n);
									xxx = 1;
									break;
								}
							}
							if (xxx == 1) {
								break;
							}
						}
					}
				}
			}
		}
		if (flag == 0) {
			break;
		}

	}

}
int main() {


	vector<vector<edge>> allocations;
	vector<vector<edge>> P_allocations;
	vector<graph> graphicMatroids;
	vector<edge> edgeSet;
	intialGraphicMatroids3(graphicMatroids);
	generateEdgeSet(vN, &edgeSet);
	
//	int MSS1 =computeMMSi(graphicMatroids[0], 9);
	MMS_Allocations(graphicMatroids, allocations);
	//将剩余的物品加给第一个玩家
	for (int i = 0; i < vN * (vN - 1) / 2; i++) {
		int flag = 0;
		for (int j = 0; j < allocations.size(); j++) {
			for (int k = 0; k < allocations[j].size(); k++) {
				if (isEdgeTheSame(edgeSet[i], allocations[j][k])) {
					flag = 1;
				}
			}
		}
		if (flag == 0) {
			allocations[0].push_back(edgeSet[i]);
		}
	}
	PMMS_Allocations(graphicMatroids, P_allocations);
	cout << "----------------------------------------------------------------------------------------\n";
	cout << "MMS Allocations!\n";
	printRawAllocations(allocations);
	cout << endl;
	cout << "----------------------------------------------------------------------------------------\n";
	cout << "PMMS Allocations!\n";
	printRawAllocations(P_allocations);
	cout << endl;
	cout << "----------------------------------------------------------------------------------------\n";
	cout << "MMS Allocations!\n";
	printAllocations(allocations);
	cout << endl;
//	PMMS_Allocations(graphicMatroids, P_allocations);
	cout << "----------------------------------------------------------------------------------------\n";
	cout << "PMMS Allocations!\n";
	printAllocations(P_allocations);
	cout << endl;
	cout << "----------------------------------------------------------------------------------------\n";
	cout << "MMS Values!\n";
	for(int i=0;i<allocations.size();i++)
		cout << computeRank(allocations[i], graphicMatroids[i])<<" ";
	cout << "\n";
	cout << "----------------------------------------------------------------------------------------\n";
	cout << "MMS Values!\n";
	for (int i = 0; i < P_allocations.size(); i++)
		cout << computeRank(P_allocations[i], graphicMatroids[i]) << " ";
	cout << "\n";

	
}


