/*
 * RunTest.cpp
 *
 *  Created on: Jan 1, 2016
 *      Author: Rohit
 */
#include <iostream>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include "hyperloglog.hpp"
#include "modifiedhyperloglog.hpp"
#include "influenceSetApprox.hpp"
#include "Properties.h"
#include "Split.h"

using namespace std;
using namespace hll;
using namespace mhll;
using namespace isa;
int main() {
	Properties props;

	props.Read("config.properties");

	string graphFolder, graphFiles, windows, windowtype, seed, outputFolder;
	string buckets;

	props.GetValue("graphFolder", graphFolder);
	props.GetValue("outputFolder", outputFolder);
	props.GetValue("graphFiles", graphFiles);
	props.GetValue("windows", windows);
	props.GetValue("windowtype", windowtype);
	props.GetValue("seed", seed);

	props.GetValue("l", buckets);
	uint8_t l = atoi(buckets.c_str());
	vector<string> swindow;
	vector<int> iwindow;
	swindow = Tools::Split(windows, ',');
	for (int i = 0; i < swindow.size(); i++) {
		iwindow.push_back(atoi(swindow[i].c_str()));
	}
	vector<string> sseed;
	vector<int> iseed;
	sseed = Tools::Split(seed, ',');
	for (int i = 0; i < sseed.size(); i++) {
		iseed.push_back(atoi(sseed[i].c_str()));
	}
	//for (int i = 10; i < 1000;) {

	//	iseed.push_back(i);
	//	i = i + 20;
	//}
	//iseed.push_back(1000);
	vector<string> files = Tools::Split(graphFiles, ',');
	long long abswindow = 0;
	for (string file : files) {
		std::cout << "*****************" << std::endl;
		std::cout << file << std::endl;
		for (int wind : iwindow) {
			//for (int wind=81;wind<100;) {
			string input = graphFolder + file + ".txt";
			stringstream w;
			w << wind;
			string output = outputFolder + file + "_" + w.str() + ".csv";
			if (windowtype.compare("h")) {
				abswindow = wind * 60 * 60;
			}
			InfluenceSetApprox my(wind, l, input, output, abswindow);
				//my.compute(true);
			//	std::cout << wind << std::endl;
			//my.testQuery(outputFolder + file, iseed);
			//wind += 5;
			my.findseed(outputFolder + file, iseed);
			//input = graphFolder + file + ".txt";

		}

	}

	//InfluenceSetApprox
	//my(0.1, 7, "C:\\phd\\testdata\\inputC\\dblp_coauthor.txt",
	//"C:\\phd\\testdata\\outputC\\dblp_coauthor.txt"
//	);
	//double t1 = time(0);
	//my.compute(true);
//	my.findseed("C:\\phd\\testdata\\outputC\\dblp_coauthor", 50);
//	std::cout << "time :" << (time(0) - t1) << std::endl;
	/*
	 string fromVertex, toVertex;
	 long etime;
	 edge tempEdge;
	 map<string, ModifiedHyperLogLog> nodes;
	 stack<edge> edges;
	 ifstream infile("C:\\phd\\testdata\\inputC\\slashdot-threads.txt");
	 double t1 = time(0);
	 std::map<string, ModifiedHyperLogLog>::iterator it;

	 while (infile >> fromVertex >> toVertex >> etime) {

	 }
	 std::cout << "time :" << (time(0) - t1) << std::endl;
	 */
	return 0;
}

void testHLL() {

	HyperLogLog hll(7);

	ModifiedHyperLogLog mll(7);
	vector<string> somedata(100);
	for (int i = 0; i < 100; ++i) {
		somedata[i] = i;
	}
	vector<std::string>::iterator iter = somedata.begin();
	vector<std::string>::iterator iter_end = somedata.end();

	;
	vector<string> somenewdata(100);
	for (int i = 0; i < 100; ++i) {
		somenewdata[i] = i;
	}
	iter = somenewdata.begin();
	iter_end = somenewdata.end();
	for (; iter != iter_end; ++iter) {
		mll.add(iter->c_str(), iter->size(), (long) 1);
	}

	double cardinality = mll.convertToHLL().estimate();
	std::cout << "Cardinality:" << cardinality << std::endl;
}
