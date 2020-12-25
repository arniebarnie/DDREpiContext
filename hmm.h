//
// Created by arnav on 12/30/18.
//

#ifndef DDREPICONTEXT_HMM_H
#define DDREPICONTEXT_HMM_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <random>
#include <unordered_map>

int DEFAULT_NUMSPLITBINS = 5000;

struct RecIntString
{
    int nindex;
    std::string sz;

    RecIntString(int n, std::string& s)
    {
        nindex = n;
        sz = s;
    }

    RecIntString() = default;
};

int RecIntStringCompare(RecIntString& a, RecIntString& b)
{
    return a.sz.compare(b.sz);
}

struct ObservedRec
{
    int nobserved;
    std::vector<bool> flagA;

public:
    ObservedRec(int n, std::vector<bool>& f) {
        nobserved = n;
        flagA = f;
    }

    ObservedRec(const ObservedRec& o)
    {
        nobserved = o.nobserved;
        flagA = o.flagA;
    }

    ObservedRec() = default;
};

struct HMM {

    double lglike;

    std::string outdir = "";

    std::string indir = "";

    std::vector<double> init;

    std::vector<std::vector<double> > transmatrix;

    std::vector<int> tmnum;

    std::vector<std::vector<int> > tmindex;

    std::vector<std::vector<bool> > elim;

    std::vector<int> tmnumCol;

    std::vector<std::vector<int> > tmindexCol;

    std::vector<std::vector<std::vector<double> > > emisprobs;

    int states;

    std::vector<std::vector<int> > dataObservedIndex;

    std::vector<std::vector<bool> > dataObservedValues;

    std::vector<std::vector<bool> > dataNotMissing;

    std::vector<std::vector<bool> > dataObservedSeqFlags;

    std::vector<std::string> datasets;

    int datasetnum;

    std::vector<std::string> cellSeq;

    std::vector<std::string> chromSeq;


    std::vector<std::string> chromfiles;

    std::vector<int> stateordering;

    std::vector<int> colordering;
};


#endif //PORT_HMM_H
