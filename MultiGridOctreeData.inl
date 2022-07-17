
#include "Octree.h"

// VertexData
long long VertexData::CenterIndex(const OctNode *node, const int &maxDepth, int idx[DIMENSION]) {
    int d,o[3];
    node->depthAndOffset(d,o);
    for(int i=0;i<DIMENSION;++i){
        idx[i]=BinaryNode<float>::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);
    }
    return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}

// SortedTreeNodes
SortedTreeNodes::SortedTreeNodes(void) {
    nodeCount=NULL;
    treeNodes=NULL;
    maxDepth=0;
}

SortedTreeNodes::~SortedTreeNodes(void) {
    if(nodeCount) delete [] nodeCount;
    if(treeNodes) delete [] treeNodes;
}

void SortedTreeNodes::set(OctNode& root, const int& setIndex) {
    if(nodeCount) delete [] nodeCount;
    if(treeNodes) delete [] treeNodes;

    maxDepth=root.maxDepth()+1;
    // record the starting address of each depth, save in [0,maxDepth()+1]
    nodeCount=new int[maxDepth+1];
    // create an array contains pointer to all nodes
    treeNodes=new OctNode*[root.nodes()];

    OctNode* temp=root.nextNode();
    int i,cnt=0;
    while(temp){
        treeNodes[cnt++]=temp;
        temp=root.nextNode(temp);
    }
    // small depth node at front
    qsort(treeNodes,cnt,sizeof(const OctNode*),OctNode::CompareForwardPointerDepths);
    for(i=0;i<maxDepth;++i){nodeCount[i]=0;}
    for(i=0;i<cnt;++i){
        if(setIndex) treeNodes[i]->nodeData.nodeIndex=i;
        nodeCount[treeNodes[i]->depth()+1]++;
    }
    for(i=1;i<=maxDepth;++i)
        nodeCount[i]+=nodeCount[i-1];

}
