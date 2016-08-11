/*
 * InfluenceSetApprox.hpp
 *
 *  Created on: Dec 31, 2015
 *      Author: Rohit
 */

#ifndef INFLUENCESETAPPROX_HPP_
#define INFLUENCESETAPPROX_HPP_

#include <iostream>
#include <map>
#include <utility>
#include "hyperloglog.hpp"
#include <stack>
#include "modifiedhyperloglog.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <algorithm>
#include <set>
#include "Timer.h"
using namespace std;
using namespace hll;
using namespace mhll;

namespace isa {
struct edge {
	string fromVertex;
	string toVertex;
	long time;
};
inline bool operator<(const edge& lhs, const edge& rhs) {
	return lhs.fromVertex < rhs.fromVertex;
}
struct node {
	string nodeid;
	HyperLogLog nodeset;
	double estimate;
};

bool sortByEstimate(const node &lhs, const node &rhs) {
	return lhs.estimate > rhs.estimate;
}
class InfluenceSetApprox {
public:
	map<string, HyperLogLog> nodesummary;
	int windowpercent;
	uint8_t numberofbuckets;
	string datafile;
	map<string, ModifiedHyperLogLog> nodes;
	string outfile;
	long long window;
	InfluenceSetApprox(int wp = 1, uint8_t b = 7, string file = "",
			string ofile = "", long long wind = 0) throw (invalid_argument) {
		windowpercent = wp;
		numberofbuckets = b;
		datafile = file;
		outfile = ofile;
		window = wind;
	}
	void compute(bool write) {
		//	cout << "window percent:: " << windowpercent << endl;
		stack<edge> edges;

		long dstart = 0;
		long dend = 0;
		int count = 0;
		string fromVertex, toVertex;
		long etime;
		vector<string> temp;
		ifstream infile(datafile.c_str());

		//	std::cout << "line:" << datafile << std::endl;
		double starttime = time(0);
		edge tempEdge;
		std::map<string, ModifiedHyperLogLog>::iterator it;
		while (infile >> fromVertex >> toVertex >> etime) {
			if (count == 0) {
				dstart = etime;
			}
			dend = etime;

			tempEdge.fromVertex = fromVertex;
			tempEdge.toVertex = toVertex;
			tempEdge.time = etime;
			edges.push(tempEdge);
			it = nodes.find(fromVertex);

			if (it == nodes.end()) {
				ModifiedHyperLogLog mhll1(numberofbuckets);
				nodes[fromVertex] = mhll1;
			}

			it = nodes.find(toVertex);

			if (it == nodes.end()) {
				ModifiedHyperLogLog mhll2(numberofbuckets);
				nodes[toVertex] = mhll2;
			}
			count++;
		}
		dend = dend / 100;
		dstart = dstart / 100;
		//cout << "dend ::" << dend << endl;
		//cout << "dstart ::" << dstart << endl;
		if (window == 0)
			window = (dend - dstart) * windowpercent;
		double endtime = time(0);
		int edgelength = edges.size();
		double timeadd = 0;
		double t1 = 0;
		double t2 = 0;
		double timeunion = 0;

		//	std::cout << "#nodes :" << nodes.size() << std::endl;
		//	std::cout << "#edges :" << edgelength << std::endl;
		//std::cout << "window :" << window << std::endl;
		//	std::cout << "time taken to read :" << (endtime - starttime)
		//			<< std::endl;
		starttime = time(0);
		while (edgelength > 0) {
			tempEdge = edges.top();
			t1 = time(0);
			nodes[tempEdge.fromVertex].add(tempEdge.toVertex.c_str(),
					tempEdge.toVertex.size(), tempEdge.time);
			t2 = time(0);
			timeadd += t2 - t1;
			nodes[tempEdge.fromVertex].merge(nodes[tempEdge.toVertex],
					tempEdge.time, window);
			t1 = time(0);
			timeunion += t1 - t2;
			edges.pop();
			edgelength--;
		}
		endtime = time(0);

		//std::cout << "time ," << (endtime - starttime) << ", window , "
		//		<< windowpercent << "," << endl;
		//	std::cout << "time taken add :" << timeadd << std::endl;
		//	std::cout << "time taken union :" << timeunion << std::endl;
		//	typedef std::map<string, HyperLogLog>::iterator nodeit;
		if (write) {
			typedef std::map<string, ModifiedHyperLogLog>::iterator nodeit;
			ofstream result;
			result.open(outfile.c_str());
			//cout<<outfile<<endl;
			for (nodeit iterator = nodes.begin(); iterator != nodes.end();
					iterator++) {
				result << iterator->first << "," << iterator->second.estimate()
						<< "\n";
				//nodesummary[iterator->first]=iterator->second.convertToHLL();

			}

			result.close();
		}

	}
	void findseed(string keyfile, vector<int> seedcount) {

		node nd;
		compute(false);
		//	std::cout << "finished compute" << std::endl;
		Platform::Timer timer;
		for (int seedc : seedcount) {

			//	std::cout << "staring" << std::endl;
			timer.Start();
			vector<node> nodelist(seedc);
			vector<string> seedlist(seedc);
			set<string> selectednodes;
			typedef std::map<string, ModifiedHyperLogLog>::iterator nodeit;
			for (nodeit iterator = nodes.begin(); iterator != nodes.end();
					iterator++) {
				nd.nodeid = iterator->first;
				nd.nodeset = iterator->second.convertToHLL();
				nd.estimate = nd.nodeset.estimate();
				if (nd.estimate > 2)
					nodelist.push_back(nd);
			}
			sort(nodelist.begin(), nodelist.end(), sortByEstimate);
			seedlist[0] = nodelist[0].nodeid;
			HyperLogLog seed(numberofbuckets);

			double max = -1.0, tempsize = 0.0;

			string newseed = "";
			HyperLogLog is(numberofbuckets), temp(numberofbuckets), maxinf(
					numberofbuckets), newinf(numberofbuckets);
			is = nodelist[0].nodeset;
			selectednodes.insert(nodelist[0].nodeid);

			for (int i = 1; i < seedc; i++) {

				max = 0.0;
				newseed = "";
				for (node &n : nodelist) {
					if (selectednodes.find(n.nodeid) == selectednodes.end()) {
						if (n.nodeset.estimate() < max) {
							break;
						}
						temp = is;
						temp.merge(n.nodeset);
						tempsize = temp.estimate();
						if (tempsize > max) {
							max = tempsize;
							maxinf = temp;
							newseed = n.nodeid;
							newinf = n.nodeset;
						}
					}
				}
				seedlist[i] = newseed;
				selectednodes.insert(newseed);
				is.merge(newinf);
			}

			//std::cout << windowpercent << "," << seedc << ","
			//		<< timer.LiveElapsedMilliseconds() << std::endl;
			ofstream result;
			stringstream resultfile;
			resultfile << keyfile << "_" << windowpercent << "_100" << ".keys";
			std::cout << resultfile.str() << endl;
			result.open(resultfile.str().c_str());
			for (string &n : seedlist) {
				result << n << "\n";

			}

			result.close();
		}

	}

	void testQuery(string keyfile, int seedcount) {
		map<string, HyperLogLog> nodelist;

		double timeused(0);
		compute(false);
		vector<string> seedlist(seedcount);
		Platform::Timer timer;
		int nodecount = 0;
		typedef std::map<string, ModifiedHyperLogLog>::iterator nodeit;
		for (nodeit iterator = nodes.begin(); iterator != nodes.end();
				iterator++) {
			nodecount++;
			nodelist[iterator->first] = iterator->second.convertToHLL();

		}

		for (int m = 1; m < 21; m++) {
			int count = seedcount * m;
			timeused = 0;
			for (int j = 0; j < 1000; j++) {
				HyperLogLog seed(numberofbuckets);
				timer.Start();
				int idx = 0;
				HyperLogLog temp(numberofbuckets);
				//	std::string s = std::to_string(42);
				for (int i = 1; i < count; i++) {
					idx = rand() % nodecount;
					//		s = std::to_string(idx);

					seed.merge(nodelist[""]);
				}
				timeused += timer.LiveElapsedMilliseconds();
			}
			timeused = timeused / 1000;
			std::cout << count << "," << timeused << std::endl;
		}

	}

};
}

#endif /* INFLUENCESETAPPROX_HPP_ */
