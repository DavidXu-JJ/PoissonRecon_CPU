//
// Created by 徐京伟 on 2022/7/15.
//

#ifndef PRACTICE_OCTREE_H
#define PRACTICE_OCTREE_H

#include "BinaryNode.h"
#include "MarchingCubes.h"

#define DIMENSION 3

class IsoNodeData{
public:
    union{
        float cornerValues[Cube::CORNERS];
        struct{
            int mcIndex;    //marching cubes index
            /**     $eSegmentCount encode the number of eSegment on the face,
             *      every face has two bits to use.
             *      0010 means faceIndex=0 has 2 edges,
             *      1101 means faceIndex=1 has 3 edges and faceIndex=0 has 1 edge
             */
            int eSegmentCount;
            /**     $eSegments[0] encode the edgeIndex add on each face,
             *      every 4 bits denote one edge(12 edges need 4 bits to recognize),
             *      every face can save at most 1 edges segment
             */
            long long eSegments[2];
        };
    };

    IsoNodeData(void);
    ~IsoNodeData(void);

    /**     get the times that addEdgeSegment() called on face $faceIndex   */
    inline int edgeCount(const int& faceIndex) const;
    /**     get two edges' indices added on face $faceIndex       */
    inline int edgeIndex(const int& faceIndex, const int& e1, const int& e2) const;
    int addEdgeSegment(const int& edgeIndex1, const int& edgeIndex2);
};

class NodeData{
public:
    static int UseIndex;
    union{
        IsoNodeData isoNode;
        struct{
            int nodeIndex;
            float centerWeightContribution;
        };
    };
    float value;

    NodeData(void);
    ~NodeData(void);
};

class OctNode
{
private:

    /**     implement F(children,$node) to all children,
     *      except &this itself                             */
    template<class NodeAdjacencyFunction>
    void __processNodeNodes(OctNode* node,NodeAdjacencyFunction* F);

    /**     this function is recursive call F->Function() depends on whether node has children  */
    template<class TerminatingNodeAdjacencyFunction>
    static void __ProcessNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                           OctNode* node1, const float& radius1,
                                           OctNode* node2, const float& radius2, const float& cWidth2,
                                           TerminatingNodeAdjacencyFunction* F);

    /**     Only call F->Function() at fixed depth.
     *      If depth < $depth, it will go deeper.    */
    template<class NodeAdjacencyFunction>
    static void __ProcessFixedDepthNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                     OctNode* node1, const float& radius1,
                                                     OctNode* node2, const float& radius2, const float& cWidth2,
                                                     const int& depth,
                                                     NodeAdjacencyFunction* F);


    /**     Used for refine node, count adjacent nodes.
     *      F->Function() will be called when depth <= $depth in this function.
     *      If depth < $depth, it will go deeper.                   */
    template<class NodeAdjacencyFunction>
    static void __ProcessMaxDepthNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                   OctNode* node1, const float& radius1,
                                                   OctNode* node2, const float& radius2, const float& cWidth2,
                                                   const int& depth,
                                                   NodeAdjacencyFunction* F);

    /**     Used to get Laplacian matrix's entries.
     *      Only successfully call F->Function() when two nodes is close enough and at same depth.
     *      node1 keeps the same, node2 will dig deeper.    */
    template<class NodeAdjacencyFunction>
    static void __ProcessTerminatingNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                      OctNode* node1, const float& radius1,
                                                      OctNode* node2, const float& radius2, const float& cWidth2,
                                                      NodeAdjacencyFunction* F);

    /**     $dir denotes the direction axis, $off denotes bigger or smaller face.
     *      no new node will be created,
     *      return the nearest node along $dir and $off direction.
     *      Depth of the return node may not be the same as &this node.             */
    const OctNode* __faceNeighbor(const int& dir,const int& off) const;
    OctNode* __faceNeighbor(const int&dir, const int& off,const int& forceChildren);

    /**     $o denotes the orientation, $i denotes the coordinate in other two dims,
     *      $idx denotes the index of other two dims/direction,
     *      return the node has the same edge denote by the params with &this node.
     *      Depth of the return node may not be the same as &this node.             */
    const OctNode* __edgeNeighbor(const int& o,const int i[2],const int idx[2]) const;
    OctNode* __edgeNeighbor(const int& o,const int i[2],const int idx[2],const int& forceChildren);

    /**     if distance between $center1 and $center2 on every dimension is smaller than $dWidth,
     *      they are determined as overlap each other.                              */
    static inline bool Overlap(const Point3D<float>& center1, const Point3D<float>& center2,const float& dWidth);
    /**     another point can be seen as (0,0,0)                                    */
    static inline bool Overlap(const float& c1, const float& c2, const float& c3, const float& dWidth);

    /**     return an encoded int,
     *      contains information about point (dx, dy, dz) with cRadius2 in [-d, d]  */
    inline static int ChildOverlap(const float& dx, const float& dy, const float& dz, const float& d, const float& cRadius2);


public:
    // use to encode the offset and depth
    static const int DepthShift, OffsetShift, OffsetShift1, OffsetShift2, OffsetShift3;
    static const int DepthMask, OffsetMask;

    static inline void DepthAndOffset(const long long& index, int& depth,int offset[3]);
    static inline void CenterAndWidth(const long long& index, Point3D<float>& center, float& width);
    static inline int Depth(const long long& index);

    OctNode* parent;
    OctNode* children;
    /**     in Function initChildren():
     *      if the node at depth n,
     *      the $off will range in [2^n-1, 2^(n+1)-2].
     *      $off value denotes the index in FunctionData's baseFunction.
     *      $off is global and offset get by depthAndOffset() is regional on its depth  */
    int d,off[3];
    NodeData nodeData;

    OctNode(void);
    ~OctNode(void);

    bool initChildren(void);
    inline int depth(void) const;
    inline Point3D<float> center(void) const;
    inline float width(void) const;

    inline void depthAndOffset(int& depth,int offset[3]) const;
    inline void centerAndWidth(Point3D<float>& center, float& width) const;

    static inline void Index(const int& depth,const int offset[3],int& d,int off[3]);

    /**     return the number of leaves under this node     */
    int leaves(void) const;

    /**     return the number of leaves under this node with depth <= $maxDepth     */
    int maxDepthLeaves(const int& maxDepth) const;

    /**     return the number of nodes under this node,
     *      including this node itself                      */
    int nodes(void) const;

    /**     return the depth of the subtree from &this node,
     *      leaves return 0                                 */
    int maxDepth(void) const;

    /**     if give a current node,
     *      return the next leaf in the subtree of &this node.
     *      if not,
     *      return the first leaf in this subtree           */
    const OctNode* nextLeaf(const OctNode* current=NULL) const;
    OctNode* nextLeaf(OctNode* current=NULL);

    const OctNode* nextNode(const OctNode* currentNode=NULL) const;
    OctNode* nextNode(OctNode* currentNode=NULL);


    /**     give a current node,
     *      return the next branch in the subtree of &this node  */
    const OctNode* nextBranch(const OctNode* current) const;
    OctNode* nextBranch(OctNode* current);

    /**     initialize a full subtree from &this node with depth $maxDepth  */
    void setFullDepth(const int& maxDepth);

    /**     implement F(children,$node) to all children,
     *      processCurrent decides where &this node is processed            */
    template<class NodeAdjacencyFunction>
    void processNodeNodes(OctNode* node, NodeAdjacencyFunction* F,const int& processCurrent=1);

    const OctNode* root(void) const;

    bool write(const char* fileName) const;
    bool write(FILE* fp) const;
    bool read(const char* fileName);
    bool read(FILE* fp);

    OctNode& operator = (const OctNode& node);

    /**     Similar to ProcessMaxDepthNodeAdjacentNodes(),
     *      but don't take depth into account/consideration.    */
    template<class NodeAdjacencyFunction>
    static void ProcessNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                         OctNode* node1, const float& radius1,
                                         OctNode* node2, const float& radius2, const float& width2,
                                         NodeAdjacencyFunction* F,const int& processCurrent=1);
    /**     Similar to ProcessMaxDepthNodeAdjacentNodes(),
     *      but don't take depth into account/consideration.    */
    template<class NodeAdjacencyFunction>
    static void ProcessNodeAdjacentNodes(OctNode* node1, const float& radius1,
                                         OctNode* node2, const float& radius2,
                                         NodeAdjacencyFunction* F,const int& processCurrent=1);

    /**     Very similar to ProcessMaxDepthNodeAdjacentNodes()  */
    template<class NodeAdjacencyFunction>
    static void ProcessFixedDepthNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                   OctNode* node1, const float& radius1,
                                                   OctNode* node2, const float& radius2, const float& width2,
                                                   const int& depth,
                                                   NodeAdjacencyFunction* F,
                                                   const int& processCurrent=1);

    template<class NodeAdjacencyFunction>
    static void ProcessFixedDepthNodeAdjacentNodes(OctNode* node1, const float& radius1,
                                                   OctNode* node2, const float& radius2,
                                                   const int& depth,
                                                   NodeAdjacencyFunction* F,
                                                   const int& processCurrent=1);

    template<class NodeAdjacencyFunction>
    static void ProcessMaxDepthNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                 OctNode* node1, const float& radius1,
                                                 OctNode* node2, const float& radius2, const float& width2,
                                                 const int& depth,
                                                 NodeAdjacencyFunction* F,
                                                 const int& processCurrent=1);

    template<class NodeAdjacencyFunction>
    static void ProcessMaxDepthNodeAdjacentNodes(OctNode* node1, const float& radius1,
                                                 OctNode* node2, const float& radius2,
                                                 const int& depth,
                                                 NodeAdjacencyFunction* F,
                                                 const int& processCurrent=1);

    template<class TerminatingNodeAdjacencyFunction>
    static void ProcessTerminatingNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                    OctNode* node1, const float& radius1,
                                                    OctNode* node2, const float& radius2, const float& width2,
                                                    TerminatingNodeAdjacencyFunction* F,
                                                    const int& processCurrent=1);

    template<class TerminatingNodeAdjacencyFunction>
    static void ProcessTerminatingNodeAdjacentNodes(OctNode* node1, const float& radius1,
                                                 OctNode* node2, const float& radius2,
                                                 TerminatingNodeAdjacencyFunction* F,
                                                 const int& processCurrent=1);


    /**     return the nearest corner index in a node center at $center */
    static int CornerIndex(const Point3D<float>& center,const Point3D<float>& p);

    const OctNode* faceNeighbor(const int& faceIndex) const;
    OctNode* faceNeighbor(const int& faceIndex, const int& forceChildren=0);

    const OctNode* edgeNeighbor(const int& edgeIndex) const;
    OctNode* edgeNeighbor(const int& edgeIndex,const int& forceChildren=0);

    const OctNode* cornerNeighbor(const int& cornerIndex) const;
    OctNode* cornerNeighbor(const int& cornerIndex,const int& forceChildren=0);

    /**     return the nearest leaf node address to $p in &this subtree */
    const OctNode* getNearestLeaf(const Point3D<float>& p) const;
    OctNode* getNearestLeaf(const Point3D<float>& p) ;

    /**     check whether the $eIndex1 edge in $node1 and the $eIndex2 edge in $node2 are on the same line
     *      (don't need to intersect)                                   */
    static bool CommonEdge(const OctNode* node1, const int& eIndex1, const OctNode* node2, const int& eIndex2);

    /**     pointer v1 and v2 is actually the pointer of pointer of OctNode
     *      (OctNode ** v1, OctNode ** v2).
     *      smaller depth will be sorted at front,
     *      greater depth will be sorted at back.                       */
    static int CompareForwardPointerDepths(const void* v1, const void* v2);
    /**     pointer v1 and v2 is actually the pointer of pointer of OctNode
     *      (OctNode ** v1, OctNode ** v2).
     *      greater depth will be sorted at front,
     *      smaller depth will be sorted at back.                       */
    static int CompareBackwardPointerDepths(const void* v1, const void* v2);

    /**     Check whether these two node is close enough      */
    static inline int Overlap2(const int& depth1, const int offSet1[DIMENSION], const float& multiplier1,
                               const int& depth2, const int offSet2[DIMENSION], const float& multiplier2);

    class Neighbors{
    public:
        /**     neighbors[1][1][1] is node itself,
         *      the offset of its neighbors in xyz axis denote its index    */
        OctNode* neighbors[3][3][3];
        Neighbors(void);
        void clear(void);
    };

    class NeighborKey{
        Neighbors* neighbors;
    public:
        NeighborKey(void);
        ~NeighborKey(void);

        /**     create $depth neighbors element in NeighborKey              */
        void set(const int& depth);
        /**     neighbors[node->depth()] will be set to be the neighbors of $node and finally return,
         *      all element of neighbors[3][3][3] will be valid.            */
        Neighbors& setNeighbors(OctNode* node);
        /**     neighbors[node->depth()] will be set to be the neighbors of $node and finally return,
         *      not all element of neighbors[3][3][3] is valid,
         *      this function won't init new children node                  */
        Neighbors& getNeighbors(OctNode* node);
    };


};



#include "Octree.inl"
#endif //PRACTICE_OCTREE_H
