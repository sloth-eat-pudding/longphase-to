#ifndef PHASINGGRAPH_H
#define PHASINGGRAPH_H

#include "Util.h"
#include "ParsingBam.h"
#include "PhasingProcess.h"


typedef std::pair<int, int> PosAllele;
typedef std::map<std::string, int> ReadBaseMap;

class SubEdge{
    
    private:
        int readCount;
        // Edge information. The vector store next pos
        // < next position, read name >
        std::map<int, std::vector<std::string> > *refRead;
        std::map<int, std::vector<std::string> > *altRead;
        // sum of edge pair quality, pos1 quality + pos2 quality
        // < next position, quality sum >
        std::map<int, int> *refQuality;
        std::map<int, int> *altQuality;
        // < next position, read count >
        std::map<int, float> *refReadCount;
        std::map<int, float> *altReadCount;

    public:
        
        SubEdge();
        ~SubEdge();
        
        void destroy();
        
        void addSubEdge(Variant &currentNode, Variant &connectNode, std::string readName, int baseQuality, double edgeWeight);
        std::pair<float,float> BestPair(int targetPos);
        float getRefReadCount(int targetPos);
        float getAltReadCount(int targetPos);        
        
        std::vector<std::string> showEdge(std::string message);
        std::vector<std::pair<int,int>> getConnectPos();
 
        int getQuality(PosAllele targetPos);
        int getAvgQuality(PosAllele targetPos);

};


struct VariantEdge{
    int currPos;
    SubEdge* alt;
    SubEdge* ref;

    VariantEdge(int currPos);
    // node pair 
    std::pair<PosAllele,PosAllele> findBestEdgePair(int targetPos, bool isONT, double diffRatioThreshold, bool debug);
    // number of read of two node. AA and AB combination
    std::pair<int,int> findNumberOfRead(int targetPos);
};


struct BlockRead{
    std::map<std::string,int> readVec;
    
    void recordRead(std::string readName);
};

struct EdgeResult{
    int rr;
    int ra;
    int ar;
    int aa;
};

enum Haplotype {
    HAPLOTYPE1 = 0,
    HAPLOTYPE2 = 1
};

class VairiantGraph{
    
    private:
        PhasingParameters *params;
        std::string *ref;
        std::vector<std::string> dotResult;
        std::vector<ReadVariant> *readVariant;
        
        // By default, a Map in C++ is sorted in increasing order based on its key.
        // position, edge
        std::map<int,VariantEdge*> *edgeList;
        // position, type
        std::map<int,int> *variantPosType;

        // position phasing result
        PosPhasingResult *posPhasingResult;
        // store phased read and read's haplotype
        std::map<std::string,int> *readHpMap;

        // produce PS tag and determine phased GT tag
        void storeResultPath();
        
        void readCorrection();

        void edgeConnectResult();
        

    public:
    
        VairiantGraph(std::string &ref, PhasingParameters &params);
        ~VairiantGraph();
    
        void addEdge(std::vector<ReadVariant> &in_readVariant);
        
        void phasingProcess(PosPhasingResult &posPhasingResult);
        void writingDotFile(std::string dotPrefix);
        std::map<std::string,int>* getReadHP();
        int totalNode();

        void destroy();
        
};




#endif
