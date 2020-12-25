//
// Created by arnav on 12/28/18.
//

#ifndef DDREPICONTEXT_UTIL_H
#define DDREPICONTEXT_UTIL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace util {
    /*
    std::istream getReader(std::string szFile)
    {
        if (boost::algorithm::ends_with(szFile, ".gz")) {
            std::ifstream file(szFile, std::ios_base::in | std::ios_base::binary);
            boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
            inbuf.push(boost::iostreams::gzip_decompressor());
            inbuf.push(file);
            return std::istream(&inbuf);
        } else {
            std::filebuf fb;
            fb.open(szFile, std::ios::in);
            return std::istream(&fb);
        }
    }
     */

    void convert10tob(int N, int b, std::string& s)
    {
        if (N == 0)
            return;
        int x = N % b;
        N /= b;
        if (x < 0)
            N += 1;
        convert10tob(N, b, s);
        s += std::to_string(x < 0 ? x + (b * -1) : x);
    }

    double euclid(std::vector<double>& xvalues, std::vector<double>& yvalues)
    {
        double ddist = 0;
        for (int ni = 0; ni < xvalues.size(); ni++)
            ddist += std::pow(xvalues[ni] - yvalues[ni], 2);
        return sqrt(ddist);
    }

    double correlation(std::vector<double>& xvalues, std::vector<double>& yvalues)
    {
        double dsumx = 0, dsumy = 0, dsumxsq = 0, dsumysq = 0, dsumxy = 0;
        double dvarx, dvary;
        int numvalues = xvalues.size();

        for (int nindex = 0; nindex < xvalues.size(); nindex++)
        {
            dsumx += xvalues[nindex];
            dsumy += yvalues[nindex];
            dsumxsq += xvalues[nindex] * xvalues[nindex];
            dsumysq += yvalues[nindex] * yvalues[nindex];
            dsumxy += xvalues[nindex] * yvalues[nindex];
        }

        if (!numvalues)
            return 0;

        dvarx = dsumxsq - dsumx * dsumx / numvalues;
        dvary = dsumysq - dsumy * dsumy / numvalues;
        double dvarxdvary = dvarx * dvary;

        if (dvarxdvary <= 0)
            return 0;
        return (dsumxy - dsumx * dsumy / numvalues) / sqrt(dvarxdvary);
    }
};


#endif //PORT_UTIL_H
