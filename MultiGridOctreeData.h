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

template<int Degree>
class Octree {
    OctNode::NeighborKey neighborKey;

    float radius;

    /**     $index should be fixed to (res2*o+o')    */
    float GetLaplacian(const int index[DIMENSION]) const;
    float GetDivergence(const int index[DIMENSION], const Point3D<float>& normal) const;

    class DivergenceFunction{
    public:
        Point3D<float> normal;
        Octree<Degree>* ot;
        int index[DIMENSION], scratch[DIMENSION];
        /**     $node2 doesn't matter   */
        void Function(OctNode* node1, OctNode* node2);
    };

    class LaplacianProjectionFunction{
    public:
        double value;
        Octree<Degree>* ot;
        int index[DIMENSION], scratch[DIMENSION];
        /**     $node2 doesn't matter   */
        void Function(OctNode* node1, OctNode* node2);
    };

    /*
    class LaplacianMatrixFunction{
    public:
        int x2,y2,z2,d2;
        Octree<Degree>* ot;
        int index[DIMENSION], scratch[DIMENSION];
        int elementCount,offset;
        MatrixEntry<float>* rowElements;    // an array record Nth element's value
        int Function(OctNode* node1, OctNode* node2);
    };
    */


    bool NonLinearUpdateWeightContribution(OctNode* node, const Point3D<float>& position);

public:
    std::vector<Point3D<float> >*normals;
    float postNormalSmooth;
    OctNode tree;
    FunctionData<Degree,double> fData;

    Octree(void);

    /**     set the fData of Octree     */
    void setFunctionData(const PPolynomial<Degree>& ReconstructionFunction,
                         const int& maxDepth,
                         const int& normalize,
                         const float& normalSmooth=-1);

};

#include "MultiGridOctreeData.inl"
#endif //PRACTICE_MULTIGRIDOCTREEDATA_H