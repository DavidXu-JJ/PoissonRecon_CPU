//
// Created by 徐京伟 on 2022/7/15.
//

#ifndef PRACTICE_MULTIGRIDOCTREEDATA_H
#define PRACTICE_MULTIGRIDOCTREEDATA_H


#include <unordered_map>
using std::unordered_map;

class VertexData {
public:

    static long long EdgeIndex(const OctNode* node,const int& eIndex,const int& maxDepth,int index[DIMENSION]);

    static long long FaceIndex(const OctNode* node,const int& fIndex,const int& maxDepth,int index[DIMENSION]);

    /**     Assume the maxDepth is set,
     *      then the resolution of valueTable is (1<<maxDepth+1).
     *      Assign $index with the index of $node's center in fData's valueTable.
     *      Value of this center is valueTable[node.off*res2+index[i]],
     *      equals fData.baseFunctions[off](index[i]).           */
    static long long CenterIndex(const OctNode* node, const int& maxDepth, int index[DIMENSION]);

    /**     Similar to CenterIndex()    */
    static long long CornerIndex(const OctNode* node, const int& cIndex, const int& maxDepth, int index[DIMENSION]);
};

class SortedTreeNodes {
public:
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
    void set(OctNode& root, const bool& setIndex);
};

/**     Similar to class SortedTreeNodes, sort in descending depth order     */
class SortedTreeLeaves{
public:
    OctNode** treeLeaves;
    int leafCount;
    SortedTreeLeaves(void);
    ~SortedTreeLeaves(void);
    void set(OctNode& root);
    void set(OctNode& root,const int& maxDepth);
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

    /**     Use to project fixed-depth solution to neighbor nodes' children */
    class LaplacianProjectionFunction{
    public:
        double value;
        Octree<Degree>* ot;
        int index[DIMENSION], scratch[DIMENSION];
        /**     $node2 doesn't matter   */
        void Function(OctNode* node1, OctNode* node2);
    };

    /**     Use to get Laplacian entries    */
    class LaplacianMatrixFunction{
    public:
        int x2,y2,z2,d2;
        /**     use to call GetLaplacian()   */
        Octree<Degree>* ot;
        int index[DIMENSION], scratch[DIMENSION];
        int elementCount,offset;
        /**     an array record Nth element's value     */
        MatrixEntry<float>* rowElements;
        int Function(OctNode* node1, OctNode* node2);
    };

    class RestrictedLaplacianMatrixFunction{
    public:
        int depth,offset[3];
        /**     use to call GetLaplacian()   */
        Octree<Degree>* ot;
        float radius;
        int index[DIMENSION], scratch[DIMENSION];
        int elementCount;
        MatrixEntry<float>* rowElements;
        int Function(const OctNode* node1, const OctNode* node2);
    };

    class PointIndexValueFunction{
    public:
        int res2;
        double* valueTables;
        int index[DIMENSION];
        float value;
//        int cnt=0;
        /**     calculate $node in which iso-value surface   */
        void Function(const OctNode* node);
    };

    class AdjacencyCountFunction{
    public:
        int adjacencyCount;
        /**     $node2 doesn't matter   */
        void Function(const OctNode* node1, const OctNode* node2);
    };

    class AdjacencySetFunction{
    public:
        int* adjacencies,adjacencyCount;
        /**     $node2 doesn't matter   */
        /**     Add index of node1 to adjacencies   */
        void Function(const OctNode* node1, const OctNode* node2);
    };

    /**     Use to create children node     */
    class RefineFunction{
    public:
        int depth;
//        int cnt=0;
        /**     Call node1->initChildren if node1->depth()<this->depth,
         *      $node2 doesn't matter   */
        void Function(OctNode* node1,const OctNode* node2);
    };

    /**     Assign matrix to be the lower triangle Laplacian matrix of all nodes with $depth    */
    int GetFixedDepthLaplacian(SparseSymmetricMatrix<float>& matrix, const int& depth, const SortedTreeNodes& sNodes);

    /**     The entries is restricted in the sub-tree of restricted Node $rNode,
     *      this will be call multiple times in SolveFixedDepthMatrix()
     *      Get their Laplacian matrix entries.                         */
    int GetRestrictedFixedDepthLaplacian(SparseSymmetricMatrix<float>& matrix,
                                         const int* entries,const int& entryCount,
                                         const OctNode* rNode,const float& radius,
                                         const SortedTreeNodes& sNodes);


    /**     Use the Lx=v Solution at $depth to project fixed-depth solution back onto deeper nodes' residual    */
    int SolveFixedDepthMatrix(const int& depth, const SortedTreeNodes& sNodes);

    /**     Use multiple small sub-trees to solve the Lx=v      */
    int SolveFixedDepthMatrix(const int& depth, const int& startingDepth, const SortedTreeNodes& sNodes);

    /**     Split the cube based on cornerValues,
     *      make the complicated node splits into small cubes.  */
    void SetIsoSurfaceCorners(const float& isoValue,const int& subdivisionDepth,const int& fullDepthIso);

    /**     check if the face of node is on the boundary face at subdivideDepth */
    static int IsBoundaryFace(const OctNode* node,const int& faceIndex,const int& subdivideDepth);

    void PreValidate(OctNode* node,const float& isoValue,const int& maxDepth,const int& subdivideDepth);

    void PreValidate(const float& isoValue,const int& maxDepth,const int& subdivideDepth);

    void Validate(OctNode* node,const float& isoValue,const int& maxDepth,const int& fullDepthIso,const int& subdivideDepth);

    /**     Subdivide the tree by multiple checks   */
    void Validate(OctNode* node,const float& isoValue,const int& maxDepth,const int& fullDepthIso);

    /**     Subdivide and delete $node.nodeData.isoData.cornerValues,
     *      eight children are assigned their own cornerValues. */
    void Subdivide(OctNode* node,const float& isoValue,const int& maxDepth);

    static int InteriorFaceRootCount(const OctNode* node,const int& faceIndex,const float& isoValue,const int maxDepth);

    /**     Count the edge root on neighbor cube        */
    static int EdgeRootCount(const OctNode* node,const int& edgeIndex,const float& isoValue,const int& maxDepth);

    /**     calculate the point $position contribution to node's neighbors   */
    int NonLinearUpdateWeightContribution(OctNode* node, const Point3D<float>& position);

    /**     get 1 / (sum of weight in $node's neighbors) of $position point  */
    float NonLinearGetSampleWeight(OctNode* node, const Point3D<float>& position);

    void NonLinearGetSampleDepthAndWeight(OctNode* node, const Point3D<float>& position,
                                          const float& samplesPerNode,
                                          float& depth,
                                          float& weight);

    /**     Update this->normals member,
     *      only update $node's neighbor node.
     *      normals[neighbors->nodeData.nodeIndex] += $normal * weight      */
    int NonLinearSplatOrientedPoint(OctNode* node, const Point3D<float>& position, const Point3D<float>& normal);

    /**     Reach the kernelDepth node in &this tree which is nearest to $point.
     *      No new node is created.
     *      The node depth updated decided by $samplePerNode, $minDepth, $maxDepth,
     *      update neighbor node of the node at decided depth.
     *      Update this->normals member.                                    */
    void NonLinearSplatOrientedPoint(const Point3D<float>& point,
                                     const Point3D<float>& normal,
                                     const int& kernelDepth,
                                     const float& samplesPerNode,
                                     const int& minDepth,
                                     const int& maxDepth);

    /**     Find whether there is valid normal in $node and its subnode     */
    int HasNormals(OctNode* node,const float& epsilon);

public:
    /**     normals[$node->nodeData.nodeIndex] denote the normal in $node   */
    std::vector<Point3D<float> >*normals;
    float postNormalSmooth;
    OctNode tree;
    FunctionData<Degree,double> fData;

    Octree(void);

    /**     set the this->fData of Octree,
     *      set the this->radius = abs(fData.polys[0].start) = 1.5  */
    void setFunctionData(const PPolynomial<Degree>& ReconstructionFunction,
                         const int& maxDepth,
                         const int& normalize,
                         const float& normalSmooth=-1);

    /**     Initialized the weight contribution of the node at each depth.
     *      Finally use the normalized oriented point in file
     *      to update this->normals of the node at kernelDepth.     */
    int setTree(char* fileName,
                const int& maxDepth,
                const int& binary,
                const int& kernelDepth,
                const float& samplesPerNode,
                const float& scaleFactor,
                Point3D<float>& center,
                float& scale,
                const int& resetSampleDepth=1);

    /**     Make the node pointer with no normal points to NULL
     *      (the memory hasn't been released?)                      */
    void ClipTree(void);

    /**     For the node with valid normal,
     *      call initChildren() to the node close to valid node
     *      until depth >= ( valid node's depth - $refineNeighbors )     */
    void finalize1(const int& refineNeighbors=-1);

    /**     Use this->normals to update nodeData.value,
     *      nodaData.value now is <divergence, Fo>.
     *      Use Length of this->normals[$index] to replace old $index node.nodeData.centerWeightContribution */
    void SetLaplacianWeights(void);

    /**     Similar to finalize1(), judge if the node is valid with <divergence, Fo>  */
    void finalize2(const int& refineNeighbors=-1);

    /**     nodeData.value: <divergence, Fo> ->
     *      solution of surface function = sum(Fo * x), o in every node   */
    int LaplacianMatrixIteration(const int& subdivideDepth);

    /**     Get iso-value from all nodes */
    float GetIsoValue(void);

//    void GetMCIsoTriangles(const float& isoValue,CoredMeshData* mesh,const int& fullDepthIso=0);
};

#include "MultiGridOctreeData.inl"
#endif //PRACTICE_MULTIGRIDOCTREEDATA_H
