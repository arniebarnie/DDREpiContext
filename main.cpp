#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <Python.h>
#include <boost/python.hpp>

#include "util.h"
#include "hmm.h"

namespace fs = boost::filesystem;
namespace py = boost::python;

const double EPSILON = std::pow(10, -300);
const std::string SZSEGMENTEXTENSION = "_segments.bed";
const double SPARSECUTOFFRATIO = 0.7;
const double SPARSECUTOFFLOOSERRATIO = 0.8;
const int BINSIZE = 200;
const int MAXITER = 1000;
const double CONVERGEDIFF = .001;
const double INFOSMOOTH = .02;
const int BUCKETS = 2;


void load(HMM& hmm)
{
    std::vector<std::string> chromfilesall;

    for (auto i = fs::directory_iterator(hmm.indir); i != fs::directory_iterator(); i++)
        if (!fs::is_directory(i->path()))
            chromfilesall.push_back(i->path().filename().string());


    for (std::string& f : chromfilesall)
        if (f.find("_binary") != std::string::npos && f[0] != '.')
            hmm.chromfiles.push_back(f);

    hmm.cellSeq = std::vector<std::string>(hmm.chromfiles.size());
    hmm.chromSeq = std::vector<std::string>(hmm.chromfiles.size());

    hmm.dataObservedIndex = std::vector<std::vector<int> >(hmm.chromfiles.size());

    std::unordered_map<int, ObservedRec> hmObserved;

    int nobserved = 0;

    std::string szLine;
    for (int nfile = 0; nfile < hmm.chromfiles.size(); nfile++)
    {
        std::cout << "reading\t" + hmm.indir + " " + hmm.chromfiles[nfile] << std::endl;

        //std::istream br = util::getReader();
        std::ifstream br(hmm.indir + "/" + hmm.chromfiles[nfile]);

        std::getline(br, szLine);

        std::vector<std::string> tokens;
        boost::split(tokens, szLine, boost::is_any_of("\t"));

        hmm.cellSeq[nfile] = tokens[0];
        hmm.chromSeq[nfile] = tokens[1];
        tokens.clear();


        std::getline(br, szLine);
        boost::split(tokens, szLine, boost::is_any_of("\t"));
        if (!nfile)
            hmm.datasets = std::vector<std::string>(tokens);
        tokens.clear();

        hmm.datasetnum = hmm.datasets.size();

        std::vector<std::string> aldata;
        int cnter = 0;
        while (std::getline(br, szLine))
        {
            szLine.erase(std::remove(szLine.begin(), szLine.end(), '\t'), szLine.end());

            boost::split(tokens, szLine, boost::is_any_of("\t"));

            aldata.push_back(szLine);
        }

        int nsize = aldata.size();
        hmm.dataObservedIndex[nfile] = std::vector<int>(nsize);

        for (int nrow = 0; nrow < nsize; nrow++)
        {
            int theBigInteger = std::stoi(aldata[nrow]);

            if (hmObserved.count(theBigInteger))
            {
                hmObserved[theBigInteger].flagA[nfile] = true;
                hmm.dataObservedIndex[nfile][nrow] = hmObserved[theBigInteger].nobserved;
            }
            else
            {
                std::vector<bool> flagA(hmm.chromfiles.size());
                flagA[nfile] = true;

                hmObserved[theBigInteger] = ObservedRec(nobserved, flagA);
                hmm.dataObservedIndex[nfile][nrow] = nobserved;

                nobserved++;
            }
        }

        br.close();
    }

    hmm.dataObservedValues = std::vector<std::vector<bool> >(nobserved, std::vector<bool>(hmm.datasetnum));
    hmm.dataNotMissing = std::vector<std::vector<bool> >(nobserved, std::vector<bool>(hmm.datasetnum));
    hmm.dataObservedSeqFlags = std::vector<std::vector<bool> >(hmm.chromfiles.size(), std::vector<bool>(nobserved));

    for (auto& pairs : hmObserved)
    {
        std::string szmapping = std::to_string(pairs.first);
        ObservedRec theObservedRec = pairs.second;

        int ncurrindex = theObservedRec.nobserved;

        int numch = szmapping.size();
        int numleading0 = hmm.datasetnum - numch;

        for (int nj = 0; nj < numleading0; nj++)
        {
            hmm.dataObservedValues[ncurrindex][nj] = false;
            hmm.dataNotMissing[ncurrindex][nj] = true;
        }

        int nmappedindex = numleading0;
        for (int nj = 0; nj < numch; nj++)
        {
            char ch = szmapping[nj];

            hmm.dataObservedValues[ncurrindex][nmappedindex] = (ch == '1');
            hmm.dataNotMissing[ncurrindex][nmappedindex] = (ch == '0' || ch == '1');
            nmappedindex++;
        }

        for (int nj = 0; nj < hmm.chromfiles.size(); nj++)
            hmm.dataObservedSeqFlags[nj][ncurrindex] = theObservedRec.flagA[nj];
    }
}

void initialize(HMM& hmm)
{
    hmm.init = std::vector<double>(hmm.states);
    hmm.emisprobs = std::vector<std::vector<std::vector<double> > >(hmm.states, std::vector<std::vector<double> >(hmm.datasetnum, std::vector<double>(BUCKETS)));
    hmm.elim = std::vector<std::vector<bool> >(hmm.states, std::vector<bool>(hmm.states));
    hmm.transmatrix = std::vector<std::vector<double> >(hmm.states, std::vector<double>(hmm.states));
    hmm.tmindex = std::vector<std::vector<int> >(hmm.states, std::vector<int>(hmm.states));
    hmm.tmnum = std::vector<int>(hmm.states);
    hmm.tmindexCol = std::vector<std::vector<int> >(hmm.states, std::vector<int>(hmm.states));
    hmm.tmnumCol = std::vector<int>(hmm.states);

    std::vector<std::vector<int> > traindataObservedIndexPair(hmm.dataObservedIndex.size());

    std::vector<std::vector<bool> > alobservedpairflags;

    std::unordered_map<std::string, int> hmObserved;

    int nobserved = 0;

    for (int nseq = 0; nseq < hmm.dataObservedIndex.size(); nseq++)
    {
        traindataObservedIndexPair[nseq] = std::vector<int>(hmm.dataObservedIndex[nseq].size() - 1);
        std::vector<bool> currvals = hmm.dataObservedValues[hmm.dataObservedIndex[nseq][0]];

        for (int nindex = 0; nindex < hmm.dataObservedIndex[nseq].size() - 1; nindex++)
        {
            std::vector<bool> nextvals = hmm.dataObservedValues[hmm.dataObservedIndex[nseq][nindex + 1]];
            std::string theBigInteger;
            for (int nmark = 0; nmark < hmm.datasetnum; nmark++)
                theBigInteger += (currvals[nmark] && nextvals[nmark] ? "1" : "0");

            int ncurrobserved = 0;
            if (!hmObserved.count(theBigInteger))
            {
                hmObserved[theBigInteger] = nobserved;
                ncurrobserved = nobserved++;

                std::vector<bool> pairflags(hmm.datasetnum);
                for (int nmark = 0; nmark < hmm.datasetnum; nmark++)
                    pairflags[nmark] = (currvals[nmark] && nextvals[nmark]);
                alobservedpairflags.push_back(pairflags);
            }
            else
                ncurrobserved = hmObserved[theBigInteger];

            traindataObservedIndexPair[nseq][nindex] = ncurrobserved;
            currvals = nextvals;

        }
    }

    int numels = alobservedpairflags.size();

    std::vector<int> tallys(numels);

    int ntotaltally = 0;

    for (int nseq = 0; nseq < traindataObservedIndexPair.size(); nseq++)
    {
        for (int nindex = 0; nindex < traindataObservedIndexPair[nseq].size(); nindex++)
            tallys[traindataObservedIndexPair[nseq][nindex]]++;
        ntotaltally += traindataObservedIndexPair[nseq].size();
    }

    for (int nj = 0; nj < hmm.datasetnum; nj++)
    {
        hmm.emisprobs[0][nj][1] = INFOSMOOTH * 1.0 / BUCKETS;
        hmm.emisprobs[0][nj][0] = 1 - hmm.emisprobs[0][nj][1];
    }

    std::vector<int> partitionTally(hmm.states);
    partitionTally[0] = ntotaltally;

    std::vector<int> initStateAssign(numels);

    std::vector<int> backptr(hmm.states);

    std::vector<std::vector<int> > nextpartitionTally(hmm.states - 1, std::vector<int>(hmm.datasetnum));

    for (int niteration = 1; niteration < hmm.states; niteration++)
    {
        for (int ni = 0; ni < niteration - 1; ni++)
            for (int nmark = 0; nmark < hmm.datasetnum; nmark++)
                nextpartitionTally[ni][nmark] = 0;

        for (int nel = 0; nel < tallys.size(); nel++)
            for (int nsplitmark = 0; nsplitmark < hmm.datasetnum; nsplitmark++)
                if (alobservedpairflags[nel][nsplitmark])
                    nextpartitionTally[initStateAssign[nel]][nsplitmark] += tallys[nel];

        double dbestinformationchange = 0;
        int nbestsplit = -1;
        int nbestsplitmark = -1;

        for (int nsplit = 0; nsplit < niteration; nsplit++)
        {
            double dprobfull = partitionTally[nsplit] / (double) ntotaltally;
            double dprobfullterm = (dprobfull > 0 ? dprobfull * std::log(dprobfull) : 0);

            for (int nsplitmark = 0; nsplitmark < hmm.datasetnum; nsplitmark++)
            {
                double dprob1 = (partitionTally[nsplit] - nextpartitionTally[nsplit][nsplitmark]) / (double) ntotaltally;
                double dprob2 = nextpartitionTally[nsplit][nsplitmark] / (double) ntotaltally;
                double dinformationchange = dprobfullterm;
                if (dprob1 > 0)
                    dinformationchange -= dprob1 * std::log(dprob1);
                if (dprob2 > 0)
                    dinformationchange -= dprob2 * std::log(dprob2);

                if (dinformationchange > dbestinformationchange)
                {
                    dbestinformationchange = dinformationchange;
                    nbestsplit = nsplit;
                    nbestsplitmark = nsplitmark;
                }
            }
        }

        int numnewsplit = nextpartitionTally[nbestsplit][nbestsplitmark];
        partitionTally[niteration] = numnewsplit;
        partitionTally[nbestsplit] -= numnewsplit;
        backptr[niteration] = nbestsplit;

        for (int nel = 0; nel < tallys.size(); nel++)
            if (initStateAssign[nel] == nbestsplit && alobservedpairflags[nel][nbestsplitmark])
                initStateAssign[nel] = niteration;
    }

    std::vector<std::vector<int> > postally(hmm.states, std::vector<int>(hmm.datasetnum));

    for (int nel = 0; nel < tallys.size(); nel++)
        for (int nmark = 0; nmark < hmm.datasetnum; nmark++)
            if (alobservedpairflags[nel][nmark])
            {
                int ncurrstate = initStateAssign[nel];
                do
                {
                    postally[ncurrstate][nmark] += tallys[nel];
                    ncurrstate = backptr[ncurrstate];
                } while (ncurrstate != 0);
            }

    std::vector<int> partitionTallySum(partitionTally.size());

    for (int nstate = 1; nstate < hmm.states; nstate++)
    {
        int ncurrstate = nstate;
        int ntotsum = 0;
        int ncurrval = partitionTally[ncurrstate];
        do
        {
            partitionTallySum[ncurrstate] += ncurrval;
            ncurrstate = backptr[ncurrstate];
        } while (ncurrstate != 0);
    }

    for (int nstate = 1; nstate < hmm.states; nstate++)
        for (int nj = 0; nj < hmm.datasetnum; nj++)
        {
            hmm.emisprobs[nstate][nj][1] = INFOSMOOTH * 1.0 / BUCKETS + (1 - INFOSMOOTH) * (postally[nstate][nj] / (double) partitionTallySum[nstate]);
            hmm.emisprobs[nstate][nj][0] = 1 - hmm.emisprobs[nstate][nj][1];
        }

    std::vector<int> numstarts(hmm.states);
    int nvalidseq = 0;
    for (int nseq = 0; nseq < traindataObservedIndexPair.size(); nseq++)
        if(!traindataObservedIndexPair[nseq].empty())
        {
            numstarts[initStateAssign[traindataObservedIndexPair[nseq][0]]]++;
            nvalidseq++;
        }

    for (int ni = 0; ni < hmm.init.size(); ni++)
        hmm.init[ni] = INFOSMOOTH * 1.0 / hmm.states + (1 - INFOSMOOTH) * numstarts[ni] / (double) nvalidseq;

    std::vector<std::vector<int> > transitiontally(hmm.states, std::vector<int>(hmm.states));
    int nnextstate;
    for (int nseq = 0; nseq < traindataObservedIndexPair.size(); nseq++)
        if (!traindataObservedIndexPair[nseq].empty())
        {
            int nprevstate = initStateAssign[traindataObservedIndexPair[nseq][0]];
            for (int nindex = 1; nindex < traindataObservedIndexPair[nseq].size(); nindex++)
            {
                nnextstate = initStateAssign[traindataObservedIndexPair[nseq][nindex]];
                transitiontally[nprevstate][nnextstate]++;
                nprevstate = nnextstate;
            }
        }

    for (int ni = 0; ni < hmm.states; ni++)
    {
        double dnumfromi = std::accumulate(transitiontally[ni].begin(), transitiontally[ni].end(), 0.0);

        for (int nj = 0; nj < hmm.states; nj++)
        {
            hmm.elim[ni][nj] = false;
            hmm.transmatrix[ni][nj] = INFOSMOOTH * 1.0 / hmm.states + (1 - INFOSMOOTH) * (transitiontally[ni][nj] / (double) dnumfromi);
            hmm.tmindex[ni][nj] = hmm.tmindexCol[ni][nj] = nj;
        }
        hmm.tmnum[ni] = hmm.tmnumCol[ni] = hmm.states;
    }
}

void orderCorrel(HMM& hmm, std::vector<std::vector<double> >& data, std::vector<int>& ordering)
{
    std::vector<bool> assignedrow(ordering.size());
    std::vector<int> temproworder(ordering.size());
    double dmintotsum = std::numeric_limits<double>::infinity();

    std::vector<std::vector<double> > correlationdistance(ordering.size(), std::vector<double>(ordering.size()));
    for (int ni = 0; ni < ordering.size(); ni++)
        for (int nj = 0; nj < ordering.size(); nj++)
            correlationdistance[ni][nj] = std::sqrt(1 - util::correlation(data[ni], data[nj]));

    for (int ninitrow = 0; ninitrow < ordering.size(); ninitrow++)
    {
        temproworder[0] = ninitrow;

        std::fill(assignedrow.begin(), assignedrow.end(), false);
        assignedrow[ninitrow] = true;

        double dtotsum = 0;
        int nminrow;
        int nprevminrow = ninitrow;

        for (int ncurrow = 1; ncurrow < ordering.size(); ncurrow++)
        {
            double dmindist = std::numeric_limits<double>::infinity();
            nminrow = 0;

            for (int nrow = 0; nrow < ordering.size(); nrow++)
                if (!assignedrow[nrow])
                    if (correlationdistance[nrow][nprevminrow] < dmindist)
                    {
                        dmindist = correlationdistance[nrow][nprevminrow];
                        nminrow = nrow;
                    }

            dtotsum += dmindist;
            temproworder[ncurrow] = nminrow;
            assignedrow[nminrow] = true;
            nprevminrow = nminrow;
        }
        if (dtotsum < dmintotsum)
        {
            dmintotsum = dtotsum;
            std::copy(temproworder.begin(), temproworder.end(), ordering.begin());
        }
    }
}

void orderStates(HMM& hmm)
{
    std::vector<std::vector<double> > emissionprobspos(hmm.states, std::vector<double>(hmm.datasetnum));
    for (int ni = 0; ni < hmm.states; ni++)
        for (int nj = 0; nj < hmm.datasetnum; nj++)
            emissionprobspos[ni][nj] = hmm.emisprobs[ni][nj][1];
    orderCorrel(hmm, emissionprobspos, hmm.stateordering);
}

void orderCol(HMM& hmm)
{
    std::vector<std::vector<double> > emissionprobspostranspose(hmm.datasetnum, std::vector<double>(hmm.states));
    for (int ni = 0; ni < hmm.datasetnum; ni++)
        for (int nj = 0; nj < hmm.states; nj++)
            emissionprobspostranspose[ni][nj] = hmm.emisprobs[nj][ni][1];
    orderCorrel(hmm, emissionprobspostranspose, hmm.colordering);
}

void printE(HMM& hmm, int niteration)
{
    std::ofstream pw(hmm.outdir + "/emissions_" + std::to_string(hmm.states) + ".txt");

    if (niteration <= 1) {
        std::cout << "Writing to file " << hmm.outdir + "/emissions_" + std::to_string(hmm.states) + ".txt" << std::endl;
//        for (int ni = 0; ni < hmm.datasets.size(); ni++)
//            std::cout << "\t" << hmm.datasets[ni];
//        std::cout << std::endl;
    }

    pw << "state (Emission order)";
    for (int ni = 0; ni < hmm.datasets.size(); ni++) {
        pw << "\t" << hmm.datasets[hmm.colordering[ni]];
        //std::cout << "\t" << hmm.datasets[hmm.colordering[ni]];
    }
    pw << std::endl;
    //std::cout << std::endl;

    for (int ni = 0; ni < hmm.states; ni++)
    {
        pw << ni + 1;
        for (int nj = 0; nj < hmm.emisprobs[ni].size(); nj++)
            pw << "\t" << hmm.emisprobs[hmm.stateordering[ni]][hmm.colordering[nj]][1];
        pw << std::endl;
    }
    pw << std::endl;
    pw.close();
}

void printP(HMM& hmm, int niteration)
{
    std::ofstream pw(hmm.outdir + "/model_" + std::to_string(hmm.states) + ".txt");

    if (niteration <= 1)
        std::cout << "Writing to file " << hmm.outdir + "/model_" + std::to_string(hmm.states) + ".txt" << std::endl;

    pw << hmm.states << "\t" << hmm.datasetnum << "\tE\t" << hmm.lglike << "\t" << niteration << std::endl;

    for (int ni = 0; ni < hmm.states; ni++)
        pw << "init\t" << ni + 1 << "\t" << hmm.init[hmm.stateordering[ni]] << std::endl;

    for (int ni = 0; ni < hmm.transmatrix.size(); ni++)
        for (int nj = 0; nj < hmm.transmatrix[hmm.stateordering[ni]].size(); nj++)
            pw << "transmatrix\t" << ni + 1 << "\t" << nj + 1 << "\t" << hmm.transmatrix[hmm.stateordering[ni]][hmm.stateordering[nj]] << std::endl;

    for (int ni = 0; ni < hmm.states; ni++)
        for (int nj = 0; nj < hmm.emisprobs[hmm.stateordering[ni]].size(); nj++)
            for (int nk = 0; nk < hmm.emisprobs[hmm.stateordering[ni]][nj].size(); nk++)
                pw << "emisprobs\t" << ni + 1 << "\t" << nj << "\t" << hmm.datasets[nj] << "\t" << nk << "\t" << hmm.emisprobs[hmm.stateordering[ni]][nj][nk] << std::endl;

    pw.close();
}

void writeT(HMM& hmm, int niteration)
{
    std::ofstream pw(hmm.outdir + "/transitions_" + std::to_string(hmm.states) + ".txt");

    if (niteration <= 1)
        std::cout << "Writing to file " << hmm.outdir + "/transitions_" + std::to_string(hmm.states) + ".txt" << std::endl;

    pw << "state (from\\to) (Emission order)";
    for (int ni = 0; ni < hmm.states; ni++)
        pw << "\t" << (ni + 1);
    pw << std::endl;

    for (int ni = 0; ni < hmm.states; ni++)
    {
        pw << (ni + 1);
        for (int nj = 0; nj < hmm.states; nj++)
            pw << "\t" << hmm.transmatrix[hmm.stateordering[ni]][hmm.stateordering[nj]];
        pw << std::endl;
    }
    pw.close();
}

void train(HMM& hmm)
{
    int niteration = 1;
    bool bconverged = false;
    double dzerotransitioncutoff = std::pow(10, -8);
    int nsparsecutoff = (int) (hmm.states * SPARSECUTOFFRATIO);
    int nsparsecutofflooser = (int) (hmm.states * SPARSECUTOFFLOOSERRATIO);

    std::vector<int> numtime(hmm.dataObservedIndex.size());
    std::transform(hmm.dataObservedIndex.begin(), hmm.dataObservedIndex.end(), numtime.begin(), [](auto& a) { return a.size(); });
    int nmaxtime = *std::max_element(numtime.begin(), numtime.end());

    //std::cout << "Maximum number of locations\t" << nmaxtime << std::endl;

    std::vector<std::vector<double> > emissionproducts(hmm.dataObservedValues.size(), std::vector<double>(hmm.states));
    std::vector<double> emissionproducts_scale;
    std::vector<double> tempproductbetaemiss(hmm.states);
    std::vector<std::vector<double> > alpha(nmaxtime, std::vector<double>(hmm.states));
    std::vector<double> gamma_nt(hmm.states);
    std::vector<double> beta_nt(hmm.states);
    std::vector<double> beta_ntp1(hmm.states);
    std::vector<double> scale(nmaxtime);
    std::vector<std::vector<double> > coltransitionprobs(hmm.states, std::vector<double>(hmm.states));
    std::vector<std::vector<double> > gammainitstore(hmm.dataObservedIndex.size(), std::vector<double>(hmm.states));
    std::vector<std::vector<std::vector<double> > > sxistore(hmm.dataObservedIndex.size(), std::vector<std::vector<double> >(hmm.states, std::vector<double>(hmm.states)));
    std::vector<std::vector<std::vector<std::vector<double> > > > gammaksumstore(hmm.dataObservedIndex.size(), std::vector<std::vector<std::vector<double> > >(hmm.states, std::vector<std::vector<double> >(hmm.datasetnum, std::vector<double>(BUCKETS))));
    std::vector<std::vector<double> > sumforsxi(hmm.states, std::vector<double>(hmm.states));
    std::vector<std::vector<double> > gammaObservedSum(hmm.dataObservedValues.size(), std::vector<double>(hmm.states));

    int nelim = 0;
    long long ltimeitr = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    double dprevloglike = -std::numeric_limits<double>::infinity();

    do
    {
        hmm.lglike = 0;

        for (int nseq = 0; nseq < hmm.dataObservedIndex.size(); nseq++)
        {
            for (int ns = 0; ns < gammaksumstore[nseq].size(); ns++)
                for (int nmark = 0; nmark < gammaksumstore[nseq][ns].size(); nmark++)
                    for (int nbucket = 0; nbucket < gammaksumstore[nseq][ns][nmark].size(); nbucket++)
                        gammaksumstore[nseq][ns][nmark][nbucket] = 0;

            for (int ni = 0; ni < sxistore[nseq].size(); ni++)
                for (int nj = 0; nj < sxistore[nseq][ni].size(); nj++)
                    sxistore[nseq][ni][nj] = 0;

            for (int ncombo = 0; ncombo < gammaObservedSum.size(); ncombo++)
                for (int ns = 0; ns < gammaObservedSum[ncombo].size(); ns++)
                    gammaObservedSum[ncombo][ns] = 0;

            for (int ni = 0; ni < emissionproducts.size(); ni++)
                if (hmm.dataObservedSeqFlags[nseq][ni])
                {
                    bool ballzero = true;

                    for (int ns = 0; ns < hmm.states; ns++)
                    {
                        emissionproducts[ni][ns] = 1;

                        for (int nmod = 0; nmod < hmm.datasetnum; nmod++)
                            if (hmm.dataNotMissing[ni][nmod])
                                emissionproducts[ni][ns] *= (hmm.dataObservedValues[ni][nmod] ? hmm.emisprobs[ns][nmod][1] : hmm.emisprobs[ns][nmod][0]);

                         if (emissionproducts[ni][ns] >= EPSILON)
                             ballzero = false;
                    }

                    if (ballzero)
                        for (int ns = 0; ns < hmm.states; ns++)
                            emissionproducts[ni][ns] = EPSILON;
                }

            for (int ns = 0; ns < hmm.states; ns++)
                alpha[0][ns] = hmm.init[ns] * emissionproducts[hmm.dataObservedIndex[nseq][0]][ns];

            scale[0] = std::accumulate(alpha[0].begin(), alpha[0].end(), 0.0);

            std::transform(alpha[0].begin(), alpha[0].end(), alpha[0].begin(), [&](auto d) {return d / scale[0]; });

            hmm.lglike += std::log(scale[0]);

            for (int ni = 0; ni < hmm.states; ni++)
                for (int nj = 0; nj < hmm.states; nj++)
                    coltransitionprobs[ni][nj] = hmm.transmatrix[nj][ni];

            for (int nt = 1; nt < numtime[nseq]; nt++)
            {
                scale[nt] = 0;
                for (int ns = 0; ns < hmm.states; ns++)
                {
                    double dtempsum = 0;
                    if (hmm.tmnumCol[ns] < nsparsecutoff)
                        for (int nj = 0; nj < hmm.tmnumCol[ns]; nj++)
                            dtempsum += coltransitionprobs[ns][hmm.tmindexCol[ns][nj]] * alpha[nt - 1][hmm.tmindexCol[ns][nj]];
                    else
                        for (int nj = 0; nj < hmm.states; nj++)
                            dtempsum += coltransitionprobs[ns][nj] * alpha[nt - 1][nj];

                    alpha[nt][ns] = dtempsum * emissionproducts[hmm.dataObservedIndex[nseq][nt]][ns];
                    scale[nt] += alpha[nt][ns];
                }

                std::transform(alpha[nt].begin(), alpha[nt].end(), alpha[nt].begin(), [&](auto d) {return d / scale[nt]; });
                hmm.lglike += std::log(scale[nt]);
            }

            int nlastindex = numtime[nseq] - 1;
            double dinitval = 1.0 / scale[nlastindex];

            std::fill(beta_ntp1.begin(), beta_ntp1.end(), dinitval);

            double ddenom = 0;

            for (int ns = 0; ns < gamma_nt.size(); ns++)
            {
                ddenom += alpha[nlastindex][ns] * beta_ntp1[ns];
                gamma_nt[ns] = alpha[nlastindex][ns] * beta_ntp1[ns];
            }

            std::transform(gamma_nt.begin(), gamma_nt.end(), gamma_nt.begin(), [&](auto d) {return d / ddenom; });

            std::transform(gammaObservedSum[hmm.dataObservedIndex[nseq][nlastindex]].begin(), gammaObservedSum[hmm.dataObservedIndex[nseq][nlastindex]].end(), gamma_nt.begin(), gammaObservedSum[hmm.dataObservedIndex[nseq][nlastindex]].begin(), std::plus<>());

            for (int nt = nlastindex - 1; nt >= 0; nt--)
            {
                for (int ns = 0; ns < hmm.states; ns++)
                    tempproductbetaemiss[ns] = beta_ntp1[ns] * emissionproducts[hmm.dataObservedIndex[nseq][nt + 1]][ns];

                double dsumbeta = 0;

                for (int ni = 0; ni < hmm.states; ni++)
                {
                    double dtempsum = 0;

                    if (hmm.tmnum[ni] < nsparsecutoff)
                        for (int nj = 0; nj < hmm.tmnum[ni]; nj++)
                            dtempsum += hmm.transmatrix[ni][hmm.tmindex[ni][nj]] * tempproductbetaemiss[hmm.tmindex[ni][nj]];
                    else
                        for (int nj = 0; nj < hmm.states; nj++)
                            dtempsum += hmm.transmatrix[ni][nj] * tempproductbetaemiss[nj];

                    if (dtempsum / scale[nt] > std::numeric_limits<double>::infinity())
                        beta_nt[ni] = std::numeric_limits<double>::infinity();
                    else
                        beta_nt[ni] = dtempsum / scale[nt];
                }

                ddenom = 0;

                for (int ns = 0; ns < gamma_nt.size(); ns++)
                {
                    ddenom += alpha[nt][ns] * beta_nt[ns];
                    gamma_nt[ns] = alpha[nt][ns] * beta_nt[ns];
                }

                std::transform(gamma_nt.begin(), gamma_nt.end(), gamma_nt.begin(), [&](auto d) {return d / ddenom; });

                for (int ns = 0; ns < hmm.states; ns++)
                    gammaObservedSum[hmm.dataObservedIndex[nseq][nt]][ns] += gamma_nt[ns];

                double dsum = 0;

                for (int ni = 0; ni < hmm.states; ni++)
                {
                    if (hmm.tmnum[ni] < nsparsecutofflooser)
                        for (int nj = 0; nj < hmm.tmnum[ni]; nj++)
                        {
                            dsum += hmm.transmatrix[ni][hmm.tmindex[ni][nj]] * alpha[nt][ni] * tempproductbetaemiss[hmm.tmindex[ni][nj]];
                            sumforsxi[ni][hmm.tmindex[ni][nj]] = hmm.transmatrix[ni][hmm.tmindex[ni][nj]] * alpha[nt][ni] * tempproductbetaemiss[hmm.tmindex[ni][nj]];
                        }
                    else
                        for (int nj = 0; nj < hmm.states; nj++)
                        {
                            dsum += hmm.transmatrix[ni][nj] * alpha[nt][ni] * tempproductbetaemiss[nj];
                            sumforsxi[ni][nj] = hmm.transmatrix[ni][nj] * alpha[nt][ni] * tempproductbetaemiss[nj];
                        }
                }

                for (int ni = 0; ni < hmm.states; ni++)
                {
                    if (hmm.tmnum[ni] < nsparsecutoff)
                        for (int nj = 0; nj < hmm.tmnum[ni]; nj++)
                            sxistore[nseq][ni][hmm.tmindex[ni][nj]] += sumforsxi[ni][hmm.tmindex[ni][nj]] / dsum;
                    else
                        for (int nj = 0; nj < hmm.states; nj++)
                            sxistore[nseq][ni][nj] += sumforsxi[ni][nj] / dsum;
                }
                beta_ntp1 = beta_nt;
            }

            for (int ns = 0; ns < hmm.states; ns++)
                gammainitstore[nseq][ns] = gamma_nt[ns];

            for (int nindex = 0; nindex < gammaObservedSum.size(); nindex++)
                if (hmm.dataObservedSeqFlags[nseq][nindex])
                    for (int ns = 0; ns < hmm.states; ns++)
                        for (int nmark = 0; nmark < hmm.states; nmark++)
                            if (hmm.dataNotMissing[nindex][nmark])
                            {
                                if (hmm.dataObservedValues[nindex][nmark])
                                    gammaksumstore[nseq][ns][nmark][1] += gammaObservedSum[nindex][ns];
                                else
                                    gammaksumstore[nseq][ns][nmark][0] += gammaObservedSum[nindex][ns];
                            }

            if (niteration > 1 || nseq == hmm.dataObservedIndex.size() - 1)
            {
                double dsum = 0;

                for (int ni = 0; ni < hmm.states; ni++)
                {
                    double dgammainitsum = 0;
                    for (int nitr = 0; nitr < hmm.dataObservedIndex.size(); nitr++)
                        dgammainitsum += gammainitstore[nitr][ni];

                    hmm.init[ni] = dgammainitsum;
                    dsum += dgammainitsum;
                }

                std::transform(hmm.init.begin(), hmm.init.end(), hmm.init.begin(), [&](auto d) {return d / dsum; });

                bool bchange = false;
                for (int ni = 0; ni < hmm.transmatrix.size(); ni++)
                {
                    dsum = 0;

                    for (int nj = 0; nj < hmm.tmnum[ni]; nj++)
                    {
                        hmm.transmatrix[ni][hmm.tmindex[ni][nj]] = 0;
                        for (int nitr = 0; nitr < hmm.dataObservedIndex.size(); nitr++)
                            hmm.transmatrix[ni][hmm.tmindex[ni][nj]] += sxistore[nitr][ni][hmm.tmindex[ni][nj]];

                        dsum += hmm.transmatrix[ni][hmm.tmindex[ni][nj]];
                    }

                    for (int nj = 0; nj < hmm.tmnum[ni]; nj++)
                    {
                        hmm.transmatrix[ni][hmm.tmindex[ni][nj]] /= dsum;

                        if (hmm.transmatrix[ni][hmm.tmindex[ni][nj]] < dzerotransitioncutoff && ni != hmm.tmindex[ni][nj])
                        {
                            hmm.elim[ni][hmm.tmindex[ni][nj]] = true;
                            bchange = true;
                            nelim++;
                            hmm.transmatrix[ni][hmm.tmindex[ni][nj]] = 0;
                        }
                    }
                }

                if (bchange)
                {
                    for (int ni = 0; ni < hmm.transmatrix.size(); ni++)
                    {
                        int nindex = 0;
                        ddenom = 0;

                        for (int nj = 0; nj < hmm.transmatrix[ni].size(); nj++)
                            if (!hmm.elim[ni][nj])
                            {
                                hmm.tmindex[ni][nindex++] = nj;
                                ddenom += hmm.transmatrix[ni][nj];
                            }

                        std::transform(hmm.transmatrix[ni].begin(), hmm.transmatrix[ni].end(), hmm.transmatrix[ni].begin(), [&](auto d) {return d / ddenom; });

                        hmm.tmnum[ni] = nindex;
                    }

                    for (int ni = 0; ni < hmm.transmatrix.size(); ni++)
                    {
                        int nindex = 0;
                        for (int nj = 0; nj < hmm.transmatrix[ni].size(); nj++)
                            if (!hmm.elim[nj][ni])
                                hmm.tmindexCol[ni][nindex++] = nj;

                        hmm.tmnumCol[ni] = nindex;
                    }
                }

                for (int ns = 0; ns < hmm.states; ns++)
                    for (int nmark = 0; nmark < hmm.emisprobs[ns].size(); nmark++)
                    {
                        double dgammadenom = 0;

                        for (int nbucket = 0; nbucket < BUCKETS; nbucket++)
                        {
                            hmm.emisprobs[ns][nmark][nbucket] = 0;

                            for (int nitr = 0; nitr < hmm.dataObservedIndex.size(); nitr++)
                                hmm.emisprobs[ns][nmark][nbucket] += gammaksumstore[nitr][ns][nmark][nbucket];

                            dgammadenom += hmm.emisprobs[ns][nmark][nbucket];
                        }

                        if (dgammadenom > 0)
                            std::transform(hmm.emisprobs[ns][nmark].begin(), hmm.emisprobs[ns][nmark].end(), hmm.emisprobs[ns][nmark].begin(), [&](auto d) {return d / dgammadenom; });
                    }
            }

        }

        double ddiff = hmm.lglike - dprevloglike;

        dprevloglike = hmm.lglike;

        orderStates(hmm);
        orderCol(hmm);

        writeT(hmm, niteration);
        printE(hmm, niteration);
        printP(hmm, niteration);

        long long ltimefinal = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        double dtimechange = (ltimefinal - ltimeitr) / 1000.0;
        bconverged = (((niteration >= MAXITER) || ((ddiff < CONVERGEDIFF) && (CONVERGEDIFF >= 0))));

        if (niteration == 1)
        {
            //std::cout << "Iteration: " + std::to_string(niteration) + "\tTime: " + std::to_string((long long) std::round(dtimechange)) + "s" << std::endl;
            printf("%10s %25s %10s %20s\n", "Iteration", "Estimated Log Likelihood", "Change", "Total Time (secs)");
            printf("%10s %25s %10s %20s\n", std::to_string(niteration).c_str(), std::to_string(hmm.lglike).c_str(), "-", std::to_string(dtimechange).c_str());
        }
        else
        {
            //std::cout << "Iteration: " + std::to_string(niteration) + "\tTime: " + std::to_string((long long) std::round(dtimechange)) + "s" << std::endl;
            printf("%10s %25s %10s %20s\n", std::to_string(niteration).c_str(), std::to_string(hmm.lglike).c_str(), std::to_string(ddiff).c_str(), std::to_string(dtimechange).c_str());
        }
        niteration++;
    } while (!bconverged);
    std::cout << std::endl;
}

void makeSegmentation(HMM& hmm)
{
    int nsparsecutoff = (int) (hmm.states * SPARSECUTOFFRATIO);

    std::vector<int> numtime(hmm.dataObservedIndex.size());

    std::transform(hmm.dataObservedIndex.begin(), hmm.dataObservedIndex.end(), numtime.begin(), [&](auto d) { return d.size(); });
    int nmaxtime = *std::max_element(numtime.begin(), numtime.end());

    std::vector<std::vector<double> > emissionproducts(hmm.dataObservedValues.size(), std::vector<double>(hmm.states));
    std::vector<double> tempproductbetaemiss(hmm.states);
    std::vector<std::vector<double> > alpha(nmaxtime, std::vector<double>(hmm.states));
    std::vector<std::vector<double> > gamma(nmaxtime, std::vector<double>(hmm.states));
    std::vector<double> beta_nt(hmm.states);
    std::vector<double> beta_ntp1(hmm.states);
    std::vector<double> scale(nmaxtime);
    std::vector<std::vector<double> > coltransitionprobs(hmm.states, std::vector<double>(hmm.states));

    std::unordered_map<std::string, std::ofstream> hmcellToFourColPW;
    std::vector<RecIntString> ordered;
    for (int nindex = 0; nindex < hmm.chromfiles.size(); nindex++)
        ordered.emplace_back(nindex, hmm.chromfiles[nindex]);
    std::sort(ordered.begin(), ordered.end(), RecIntStringCompare);

    for (int nseq = 0; nseq < hmm.dataObservedIndex.size(); nseq++)
    {
        int nordered_nseq = ordered[nseq].nindex;

        std::string szprefix;
        if (!hmm.cellSeq[nordered_nseq].empty()) szprefix += hmm.cellSeq[nordered_nseq] + "_";
        szprefix += std::to_string(hmm.states);

        if (!hmcellToFourColPW.count(hmm.cellSeq[nordered_nseq]))
        {
            std::string szsegmentoutfilename = hmm.outdir + "/" + szprefix + SZSEGMENTEXTENSION;

            std::cout << "Writing to file " << szsegmentoutfilename << std::endl;

            hmcellToFourColPW[hmm.cellSeq[nordered_nseq]] = std::ofstream(szsegmentoutfilename);
        }
        std::ofstream& pwbed = hmcellToFourColPW[hmm.cellSeq[nordered_nseq]];

        for (int ni = 0; ni < emissionproducts.size(); ni++)
        {
                if (hmm.dataObservedSeqFlags[nordered_nseq][ni])
                {

                    bool ballzero = true;

                    for (int ns = 0; ns < hmm.states; ns++)
                    {
                        double dproduct = 1;

                        for (int nmod = 0; nmod < hmm.datasetnum; nmod++)
                            if (hmm.dataNotMissing[ni][nmod])
                            {
                                if (hmm.dataObservedValues[ni][nmod])
                                    dproduct *= hmm.emisprobs[ns][nmod][1];
                                else
                                    dproduct *= hmm.emisprobs[ns][nmod][0];
                            }
                        emissionproducts[ni][ns] = dproduct;

                        if (dproduct >= EPSILON)
                            ballzero = false;
                    }

                    if (ballzero)
                        for (int ns = 0; ns < hmm.states; ns++)
                            emissionproducts[ni][ns] = EPSILON;
                }
        }

        for (int ns = 0; ns < hmm.states; ns++)
            alpha[0][ns] = hmm.init[ns] * emissionproducts[hmm.dataObservedIndex[nordered_nseq][0]][ns];
        scale[0] = std::accumulate(alpha[0].begin(), alpha[0].end(), 0.0);

        std::transform(alpha[0].begin(), alpha[0].end(), alpha[0].begin(), [&](auto d){ return d/scale[0]; });

        for (int ni = 0; ni < hmm.states; ni++)
            for (int nj = 0; nj < hmm.states; nj++)
                coltransitionprobs[ni][nj] = hmm.transmatrix[nj][ni];

        for (int nt = 1; nt < numtime[nordered_nseq]; nt++)
        {
            double dscale = 0;

            for (int ns = 0; ns < hmm.states; ns++)
            {
                double dtempsum = 0;
                if (hmm.tmnumCol[ns] < nsparsecutoff)
                    for (int nj = 0; nj < hmm.tmnumCol[ns]; nj++)
                    {
                        int nmappedindex = hmm.tmindexCol[ns][nj];
                        dtempsum += coltransitionprobs[ns][nmappedindex] * alpha[nt - 1][nmappedindex];
                    }
                else
                    for (int nj = 0; nj < hmm.states; nj++)
                        dtempsum += coltransitionprobs[ns][nj] * alpha[nt - 1][nj];

                alpha[nt][ns] = dtempsum * emissionproducts[hmm.dataObservedIndex[nordered_nseq][nt]][ns];
                dscale += alpha[nt][ns];
            }

            scale[nt] = dscale;
            std::transform(alpha[nt].begin(), alpha[nt].end(), alpha[nt].begin(), [&](auto d) { return d / scale[nt]; });

        }

        int nlastindex = numtime[nordered_nseq] - 1;
        double dinitval = 1.0 / scale[nlastindex];
        std::fill(beta_ntp1.begin(), beta_ntp1.end(), dinitval);

        int nmappedindexouter;
        double ddenom = 0;

        for (int ns = 0; ns < gamma[nlastindex].size(); ns++)
        {
            gamma[nlastindex][ns] = alpha[nlastindex][ns] * beta_ntp1[ns];
            ddenom += gamma[nlastindex][ns];
        }

        std::transform(gamma[nlastindex].begin(), gamma[nlastindex].end(), gamma[nlastindex].begin(), [&](auto d) { return d / ddenom; });

        for (int nt = nlastindex - 1; nt >= 0; nt--)
        {
            int ntp1 = nt + 1;

            std::transform(beta_ntp1.begin(), beta_ntp1.end(), emissionproducts[hmm.dataObservedIndex[nordered_nseq][ntp1]].begin(), tempproductbetaemiss.begin(), std::multiplies<>());

            for (int ni = 0; ni < hmm.states; ni++)
            {
                double dtempsum = 0;

                if (hmm.tmnum[ni] < nsparsecutoff)
                    for (int nj = 0; nj < hmm.tmnum[ni]; nj++)
                    {
                        nmappedindexouter = hmm.tmindex[ni][nj];
                        dtempsum += hmm.transmatrix[ni][nmappedindexouter] * tempproductbetaemiss[nmappedindexouter];
                    }
                else
                    for (int nj = 0; nj < hmm.states; nj++)
                        dtempsum += hmm.transmatrix[ni][nj] * tempproductbetaemiss[nj];

                double dratio = dtempsum / scale[nt];
                beta_nt[ni] = std::min(dratio, std::numeric_limits<double>::infinity());
            }

            for (int ns = 0; ns < gamma[nlastindex].size(); ns++)
                gamma[nt][ns] = alpha[nt][ns] * beta_nt[ns];

            ddenom = std::accumulate(gamma[nt].begin(), gamma[nt].end(), 0.0);
            std::transform(gamma[nt].begin(), gamma[nt].end(), gamma[nt].begin(), [&](auto d) { return d / ddenom; });

            beta_ntp1 = beta_nt;
        }

        int nstart = 0;

        double dmaxval = 0;
        int nmaxstate = 0;

        for (int ns = 0; ns < gamma[0].size(); ns++)
        {
            int nmappedstate = hmm.stateordering[ns];
            double dprob = gamma[0][nmappedstate];

            if (dprob > dmaxval)
            {
                dmaxval = dprob;
                nmaxstate = ns;
            }
        }

        int nmaxstateprev = nmaxstate;

        for (int nt = 1; nt < numtime[nordered_nseq]; nt++)
        {
            dmaxval = 0;
            nmaxstate = 0;
            for (int ns = 0; ns < gamma[nt].size(); ns++)
            {
                double dprob = gamma[nt][hmm.stateordering[ns]];

                if (dprob > dmaxval)
                {
                    dmaxval = dprob;
                    nmaxstate = ns;
                }
            }

            if (nmaxstateprev != nmaxstate)
            {
                pwbed << hmm.chromSeq[nordered_nseq] << "\t" << (nstart * BINSIZE) << "\t" << (nt * BINSIZE) << "\t" << "E" << (nmaxstateprev + 1) << std::endl;
                nstart = nt;
                nmaxstateprev = nmaxstate;
            }
        }

        int nlastcoordinate = numtime[nordered_nseq] * BINSIZE;
        pwbed << hmm.chromSeq[nordered_nseq] << "\t" << (nstart * BINSIZE) << "\t" << nlastcoordinate << "\t" << "E" << (nmaxstateprev + 1) << std::endl;
    }

    for (auto& p : hmcellToFourColPW) p.second.close();

}

void build(HMM& hmm)
{
    initialize(hmm);
    train(hmm);
}

int intersect(std::vector<std::pair<int, int> >& a, std::vector<std::pair<int, int> >& b)
{
    std::vector<std::pair<int, int> > overlap;

    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());

    int ni = 0, nj = 0;

    while (ni < a.size() && nj < b.size())
    {
        int lo = std::max(a[ni].first, b[nj].first);
        int hi = std::min(a[ni].second, b[nj].second);
        if (lo <= hi)
            overlap.emplace_back(lo, hi);

        if (a[ni].second < b[nj].second)
            ni++;
        else
            nj++;
    }

    int sum = 0;
    for (auto& pair : overlap) sum += pair.second - pair.first;

    return sum;
}

void enrichmap(std::unordered_map<std::string, std::vector<double> >& overlaptable)
{
    Py_Initialize();
    py::object mmodule = py::import("__main__");
    py::object mnamespace = mmodule.attr("__dict__");
    py::object future = py::import("__future__");
    py::object builtins = py::import("__builtin__");
    py::object sns = py::import("seaborn");
    py::object plt = py::import("matplotlib.pyplot");
    py::object pd = py::import("pandas");
    py::object np = py::import("numpy");
    py::object print = builtins.attr("print");
    py::object zip = builtins.attr("zip");
    py::object dict = builtins.attr("dict");
    py::object list = builtins.attr("list");
    py::object DataFrame = pd.attr("DataFrame");
    py::object heatmap = sns.attr("heatmap");
    py::object show = plt.attr("show");
    py::object figure = plt.attr("figure");

    py::dict table;
    for (auto& pair : overlaptable)
    {
        table[pair.first] = py::dict();

        for (int nstate = 1; nstate < pair.second.size(); nstate++)
            table[pair.first][nstate] = overlaptable[pair.first][nstate];
    }

    py::object df = DataFrame(table).attr("transpose")();

    figure(py::object(), py::make_tuple(15, 8));
    py::object graph = heatmap(df, py::object(), py::object(), "Reds");
    graph.attr("set_xlabel")("MRE States");
    graph.attr("set_ylabel")("Features");
    show();
}

void enrichment(std::string& segdir, std::string& coorddir, std::string& outdir, int numstates)
{
    std::vector<std::string> segfiles;

    for (auto i = fs::directory_iterator(segdir); i != fs::directory_iterator(); i++)
        if (!fs::is_directory(i->path()))
        {
            std::string f = i->path().filename().string();
            if (f.find(SZSEGMENTEXTENSION) != std::string::npos && f[0] != '.')
                segfiles.push_back(f);
        }

    std::unordered_map<std::string, std::vector<std::pair<int, int> > > segintervals[numstates + 1][segfiles.size()];
    std::vector<double> stateprop(numstates + 1);
    int cellsizesum = 0;

    std::string line;
    std::vector<std::string> tokens;
    for (int nfile = 0; nfile < segfiles.size(); nfile++)
    {
        std::cout << "reading\t" + segdir + " " + segfiles[nfile] << std::endl;

        //std::istream br = util::getReader();
        std::ifstream br(segdir + "/" + segfiles[nfile]);

        int mx = 0;
        while (std::getline(br, line))
        {
            boost::split(tokens, line, boost::is_any_of("\t"));

            int state = std::stoi(tokens[3].substr(1));
            std::string chr = tokens[0];
            int start = std::stoi(tokens[1]), end = std::stoi(tokens[2]);

            cellsizesum += end - start;
            stateprop[state] += end - start;

            segintervals[state][nfile][chr].push_back({start, end});

            tokens.clear();
        }
    }

    std::transform(stateprop.begin(), stateprop.end(), stateprop.begin(), [&] (auto d) { return d / cellsizesum; });

    //std::cout << "Cell Size Sum\t" << cellsizesum << std::endl;

    std::vector<std::string> coordfiles;
    for (auto i = fs::directory_iterator(coorddir); i != fs::directory_iterator(); i++)
        if (!fs::is_directory(i->path()))
        {
            std::string f = i->path().filename().string();
            if (f.find(".bed") != std::string::npos && f[0] != '.')
                coordfiles.push_back(f);
        }

    std::vector<std::string> features(coordfiles.size());
    std::transform(coordfiles.begin(), coordfiles.end(), features.begin(), [&](auto d) { return d.substr(0, d.find('.')); });

    for (auto& s : features) std::cout << s << "\t";
    std::cout << std::endl;

    std::unordered_map<std::string, std::vector<std::pair<int, int> > > featintervals;
    std::unordered_map<std::string, double> featprop;
    std::unordered_map<std::string, std::vector<double> > overlaptable;

    for (int nfile = 0; nfile < coordfiles.size(); nfile++)
    {
        std::cout << "reading\t" + coorddir + " " + coordfiles[nfile] << std::endl;

        std::ifstream br(coorddir + "/" + coordfiles[nfile]);

        while (std::getline(br, line))
        {
            boost::split(tokens, line, boost::is_any_of("\t"));

            std::string chr = tokens[0];
            int start = std::stoi(tokens[1]), end = std::stoi(tokens[2]);

            featintervals[chr].push_back({start, end});
            featprop[features[nfile]] += end - start;

            tokens.clear();
        }

        //std::cout << features[nfile] << "\t" << featprop[features[nfile]] << std::endl;
        featprop[features[nfile]] /= cellsizesum;

        for (int nstate = 1; nstate <= numstates; nstate++)
            for (int ncell = 0; ncell < segfiles.size(); ncell++)
                if (!segintervals[nstate][ncell].empty())
                {
                    int overlap = 0;

                    for (auto &pair : segintervals[nstate][ncell]) {
                        std::string chr = pair.first;

                        overlap += intersect(pair.second, featintervals[chr]);
                    }

                    if (overlaptable[features[nfile]].size() != numstates + 1)
                        overlaptable[features[nfile]].resize(numstates + 1);

                    overlaptable[features[nfile]][nstate] += overlap;
                }
        featintervals.clear();
    }

    std::ofstream pw(outdir + "/overlap_" + std::to_string(numstates) + ".txt");

    pw << "Features";
    for (int ni = 1; ni <= numstates; ni++) pw << "\t" << ni;
    pw << std::endl;

    for (auto& feat : features)
    {
        for (int ni = 1; ni <= numstates; ni++)
        {
            overlaptable[feat][ni] /= cellsizesum;
            overlaptable[feat][ni] /= stateprop[ni] * featprop[feat];
        }

        pw << feat;
        for (int i = 1; i <= numstates; i++) pw << "\t" << overlaptable[feat][i];
        pw << std::endl;
    }

    pw.close();

    enrichmap(overlaptable);
}

int main(int argc, char* argv[]) {
    std::string szprefixpath = fs::current_path().string() + "/";

    std::cout << szprefixpath << std::endl;

    if (std::string(argv[1]) == "GenParam" && argc == 5)
    {
        std::string szinputdir = argv[2];
        std::string szoutputdir = argv[3];

        if (!fs::exists(szinputdir))
        {
            std::cerr << szinputdir + " does not exist!" << std::endl;
            exit(1);
        }

        if (!fs::exists(szoutputdir))
            if (!fs::create_directory(szoutputdir))
            {
                std::cerr << szoutputdir + " does not exist and could not be created!" << std::endl;
                exit(1);
            }

        int numstates = 0;
        try {
            numstates = std::stoi(argv[4]);
        }
        catch (std::invalid_argument& e)
        {
            std::cerr << "'states' argument must be an integer!" << std::endl;
            exit(1);
        }

        HMM hmm;
        hmm.indir = szinputdir;
        hmm.outdir = szoutputdir;
        hmm.states = numstates;

        load(hmm);

        hmm.stateordering = std::vector<int>(hmm.states);
        std::iota(hmm.stateordering.begin(), hmm.stateordering.end(), 0);
        hmm.colordering = std::vector<int>(hmm.datasetnum);
        std::iota(hmm.colordering.begin(), hmm.colordering.end(), 0);

        build(hmm);
        makeSegmentation(hmm);
    }
    else if (std::string(argv[1]) == "EnrichmentScore" && argc == 6)
    {
        std::string segdir = argv[2];
        std::string coorddir = argv[3];
        std::string outdir = argv[4];
        int numstates = std::stoi(argv[5]);

        enrichment(segdir, coorddir, outdir, numstates);
    }
    else
    {
        std::cout << "ERROR: Command not found OR Incorrect # of arguments" << std::endl;
    }
}