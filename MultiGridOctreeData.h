//
// Created by 徐京伟 on 2022/7/15.
//

#ifndef PRACTICE_MULTIGRIDOCTREEDATA_H
#define PRACTICE_MULTIGRIDOCTREEDATA_H


#include <unordered_map>
using std::unordered_map;

class VertexData {
    /**     Assume the maxDepth is set,
     *      then the resolution of valueTable is (1<<maxDepth+1).
     *      Return the index of $node's center in fData's valueTable.
     *      Value of this center is valueTable[node.off*res2+returnValue],
     *      equals fData.baseFunctions[off](returnValue).
     */
    static long long CenterIndex(const OctNode* node, const int& maxDepth, int index[DIMENSION]);
};

class SortedTreeNodes {
    /**     the array to save the pointer of OctNode                */
    OctNode** treeNodes;
    /**     record the starting address offset of each depth,
     *      [0, root->maxDepth()+1] is valid.                       */
    int* nodeCount;
    /**     root->maxDepth()+1                                      */
    int maxDepth;
    SortedTreeNodes(void);
    ~SortedTreeNodes(void);
    /**     save all nodes in subtree from $root in &this object    */
    void set(OctNode& root, const int& setIndex);
};


#include "MultiGridOctreeData.inl"
#endif //PRACTICE_MULTIGRIDOCTREEDATA_H
