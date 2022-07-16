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

    /**     get the times that *addEdgeSegment* called on face $faceIndex   */
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
public:
    // use to encode the offset and depth
    static const int DepthShift, OffsetShift, OffsetShift1, OffsetShift2, OffsetShift3;
    static const int DepthMask, OffsetMask;

    static inline void DepthAndOffset(const long long& index, int& depth,int offset[3]);
    static inline void CenterAndWidth(const long long& index, Point3D<float>& center, float& width);
    static inline int Depth(const long long& index);

    OctNode* parent;
    OctNode* children;
    /**     in Function *initChildren*:
     *      if the node at depth n,
     *      the $off will range in [2^n-1, 2^(n+1)-2].
     *      $off value denotes the index in FunctionData's baseFunction
     */
    int d,off[3];
    NodeData nodeData;

    OctNode(void);
    ~OctNode(void);

    int initChildren(void);
    inline int depth(void) const;

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

    const OctNode* root(void) const;

    int write(const char* fileName) const;
    int write(FILE* fp) const;
    int read(const char* fileName);
    int read(FILE* fp);

    OctNode& operator = (const OctNode& node);


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

        /**     create $depth neighbors     */
        void set(const int& depth);
        /**     neighbors[node->depth()] will be set to be the neighbors of $node and finally return,
         *      all element of neighbors[3][3][3] will be valid.             */
        Neighbors& setNeighbors(OctNode* node);
        /**     neighbors[node->depth()] will be set to be the neighbors of $node and finally return,
         *      not all element of neighbors[3][3][3] is valid,
         *      this function won't init new children node                   */
        Neighbors& getNeighbors(OctNode* node);
    };


};



#include "Octree.inl"
#endif //PRACTICE_OCTREE_H
