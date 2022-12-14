

#include <bitset>
#include <iostream>

template<class real>
void dbg(real a){
    std::cout<< "\033[33;1m"<<std::bitset<64>(a)<<"\033[39;0m"<<std::endl;
}

// IsoNodeData
IsoNodeData::IsoNodeData(void)
{
    eSegmentCount=0;
    eSegments[0]=eSegments[1]=0;
}

IsoNodeData::~IsoNodeData() {;}

inline int IsoNodeData::edgeCount(const int& faceIndex) const
{
    return (eSegmentCount >> (faceIndex*2) ) & 3;
}

inline int IsoNodeData::edgeIndex(const int& faceIndex,
                                  const int& e1,
                                  const int& e2         )const
{
    return int((eSegments[e2]>>(4*(faceIndex*2+e1)))&15);
}

int IsoNodeData::addEdgeSegment(const int &edgeIndex1, const int &edgeIndex2)
{
    // find the face contains two edges
    int faceIndex=Cube::FaceAdjacentToEdges(edgeIndex1,edgeIndex2);

    if(faceIndex<0) return -1;
    int eCount=edgeCount(faceIndex);
    eSegments[0]|=(long long)(edgeIndex1)<<(4*(faceIndex*2+eCount));
    eSegments[1]|=(long long)(edgeIndex2)<<(4*(faceIndex*2+eCount));
    eCount++;
    eSegmentCount &= ~(3<<(faceIndex*2));
    eSegmentCount |= eCount<<(faceIndex*2);

    return faceIndex;
}

// NodeData
int NodeData::UseIndex=1;
NodeData::NodeData(void)
{
    if(UseIndex) {
        nodeIndex=-1;
        centerWeightContribution=0;
        isoNode=NULL;
    } else {
        isoNode=new IsoNodeData();
    }
    value=0;
}

NodeData::~NodeData() {
    if(isoNode) delete isoNode;
    isoNode=NULL;
}

// OctNode private member

template<class NodeAdjacencyFunction>
void OctNode::__processNodeNodes(OctNode* node,NodeAdjacencyFunction* F) {
    F->Function(&children[0],node);
    F->Function(&children[1],node);
    F->Function(&children[2],node);
    F->Function(&children[3],node);
    F->Function(&children[4],node);
    F->Function(&children[5],node);
    F->Function(&children[6],node);
    F->Function(&children[7],node);
    if(children[0].children){children[0].__processNodeNodes(node,F);}
    if(children[1].children){children[1].__processNodeNodes(node,F);}
    if(children[2].children){children[2].__processNodeNodes(node,F);}
    if(children[3].children){children[3].__processNodeNodes(node,F);}
    if(children[4].children){children[4].__processNodeNodes(node,F);}
    if(children[5].children){children[5].__processNodeNodes(node,F);}
    if(children[6].children){children[6].__processNodeNodes(node,F);}
    if(children[7].children){children[7].__processNodeNodes(node,F);}
}


template<class TerminatingNodeAdjacencyFunction>
void OctNode::__ProcessNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                         OctNode* node1, const float& radius1,
                                         OctNode* node2, const float& radius2, const float& cWidth2,
                                         TerminatingNodeAdjacencyFunction* F)
{
    /**     SetLaplacianWeights(): DivergenceFunction
     *      (dx, dy, dz) is the vector(node1->node2)
     *      radius1:    1.5 * width of node1
     *      radius2:    1.5 * width of node2
     *      cWidth2:    radius of node2 = width of node2 / 2    */

    /**     SolveFixedDepthMatrix(): LaplacianProjectionFunction
     *      radius1:    1.4449  * width
     *      radius2:    1.0001  * width
     *      width2:     radius = width / 2                      */

    /**     $cWidth takes half of node2's radius because it's convenient to
     *      move the center position to generate new (dx, dy, dz).
     *      Another half of node2's radius is offset by $radius.    */
    float cWidth=cWidth2/2;
    float radius=radius2/2;
    int o=ChildOverlap(dx,dy,dz,radius1+radius,cWidth);
    if(o) {
        float dx1=dx-cWidth;
        float dx2=dx+cWidth;
        float dy1=dy-cWidth;
        float dy2=dy+cWidth;
        float dz1=dz-cWidth;
        float dz2=dz+cWidth;
        if(o&  1){
            F->Function(&node2->children[0],node1);
            if(node2->children[0].children)
                __ProcessNodeAdjacentNodes(dx1,dy1,dz1,node1,radius1,&node2->children[0],radius,cWidth,F);
        }
        if(o&  2){
            F->Function(&node2->children[1],node1);
            if(node2->children[1].children)
                __ProcessNodeAdjacentNodes(dx2,dy1,dz1,node1,radius1,&node2->children[1],radius,cWidth,F);
        }
        if(o&  4){
            F->Function(&node2->children[2],node1);
            if(node2->children[2].children)
                __ProcessNodeAdjacentNodes(dx1,dy2,dz1,node1,radius1,&node2->children[2],radius,cWidth,F);
        }
        if(o&  8){
            F->Function(&node2->children[3],node1);
            if(node2->children[3].children)
                __ProcessNodeAdjacentNodes(dx2,dy2,dz1,node1,radius1,&node2->children[3],radius,cWidth,F);
        }
        if(o& 16){
            F->Function(&node2->children[4],node1);
            if(node2->children[4].children)
                __ProcessNodeAdjacentNodes(dx1,dy1,dz2,node1,radius1,&node2->children[4],radius,cWidth,F);
        }
        if(o& 32){
            F->Function(&node2->children[5],node1);
            if(node2->children[5].children)
                __ProcessNodeAdjacentNodes(dx2,dy1,dz2,node1,radius1,&node2->children[5],radius,cWidth,F);
        }
        if(o& 64){
            F->Function(&node2->children[6],node1);
            if(node2->children[6].children)
                __ProcessNodeAdjacentNodes(dx1,dy2,dz2,node1,radius1,&node2->children[6],radius,cWidth,F);
        }
        if(o&128){
            F->Function(&node2->children[7],node1);
            if(node2->children[7].children)
                __ProcessNodeAdjacentNodes(dx2,dy2,dz2,node1,radius1,&node2->children[7],radius,cWidth,F);
        }
    }
}

template<class NodeAdjacencyFunction>
void OctNode::__ProcessFixedDepthNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                   OctNode* node1, const float& radius1,
                                                   OctNode* node2, const float& radius2, const float& cWidth2,
                                                   const int& depth,
                                                   NodeAdjacencyFunction* F)
{
    /**     $cWidth takes half of node2's radius because it's convenient to
     *      move the center position to generate new (dx, dy, dz).
     *      Another half of node2's radius is offset by $radius.    */
    float cWidth=cWidth2/2;
    float radius=radius2/2;

    int o=ChildOverlap(dx,dy,dz,radius1+radius,cWidth);
    /**     Above function used to check which children of node2 is close enough to node1   */

    if(o){
        float dx1=dx-cWidth;
        float dx2=dx+cWidth;
        float dy1=dy-cWidth;
        float dy2=dy+cWidth;
        float dz1=dz-cWidth;
        float dz2=dz+cWidth;
        if(node2->depth()==depth){
            if(o&  1)   F->Function(&node2->children[0],node1);
            if(o&  2)   F->Function(&node2->children[1],node1);
            if(o&  4)   F->Function(&node2->children[2],node1);
            if(o&  8)   F->Function(&node2->children[3],node1);
            if(o& 16)   F->Function(&node2->children[4],node1);
            if(o& 32)   F->Function(&node2->children[5],node1);
            if(o& 64)   F->Function(&node2->children[6],node1);
            if(o&128)   F->Function(&node2->children[7],node1);
        } else{
            if(o&  1) if(node2->children[0].children)
                    __ProcessFixedDepthNodeAdjacentNodes(dx1,dy1,dz1,node1,radius1,
                                                       &node2->children[0],radius,cWidth,
                                                       depth,F);
            if(o&  2) if(node2->children[1].children)
                    __ProcessFixedDepthNodeAdjacentNodes(dx2,dy1,dz1,node1,radius1,
                                                       &node2->children[1],radius,cWidth,
                                                       depth,F);
            if(o&  4) if(node2->children[2].children)
                    __ProcessFixedDepthNodeAdjacentNodes(dx1,dy2,dz1,node1,radius1,
                                                       &node2->children[2],radius,cWidth,
                                                       depth,F);
            if(o&  8) if(node2->children[3].children)
                    __ProcessFixedDepthNodeAdjacentNodes(dx2,dy2,dz1,node1,radius1,
                                                       &node2->children[3],radius,cWidth,
                                                       depth,F);
            if(o& 16) if(node2->children[4].children)
                    __ProcessFixedDepthNodeAdjacentNodes(dx1,dy1,dz2,node1,radius1,
                                                       &node2->children[4],radius,cWidth,
                                                       depth,F);
            if(o& 32) if(node2->children[5].children)
                    __ProcessFixedDepthNodeAdjacentNodes(dx2,dy1,dz2,node1,radius1,
                                                       &node2->children[5],radius,cWidth,
                                                       depth,F);
            if(o& 64) if(node2->children[6].children)
                    __ProcessFixedDepthNodeAdjacentNodes(dx1,dy2,dz2,node1,radius1,
                                                       &node2->children[6],radius,cWidth,
                                                       depth,F);
            if(o&128) if(node2->children[7].children)
                    __ProcessFixedDepthNodeAdjacentNodes(dx2,dy2,dz2,node1,radius1,
                                                       &node2->children[7],radius,cWidth,
                                                       depth,F);
        }
    }
}

template<class NodeAdjacencyFunction>
void OctNode::__ProcessMaxDepthNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                               OctNode* node1, const float& radius1,
                                               OctNode* node2, const float& radius2, const float& cWidth2,
                                               const int& depth,
                                               NodeAdjacencyFunction* F)
{
    /**     $cWidth takes half of node2's radius because it's convenient to
     *      move the center position to generate new (dx, dy, dz).
     *      Another half of node2's radius is offset by $radius.    */
    float cWidth=cWidth2/2;
    float radius=radius2/2;

    /**     (dx, dy, dz) is the vector(center1->center2)
     *      radius1:    width of node1 * 3
     *      radius:     radius of node2 / 2 = width of node2 / 4
     *      cWidth:     radius of node2 / 2 = width of node2 / 4     */
    int o=ChildOverlap(dx,dy,dz,radius1+radius,cWidth);
    /**     Above function used to check which children of node2 is close enough to node1   */

    if(o){
        float dx1=dx-cWidth;
        float dx2=dx+cWidth;
        float dy1=dy-cWidth;
        float dy2=dy+cWidth;
        float dz1=dz-cWidth;
        float dz2=dz+cWidth;
        if(node2->depth()<=depth){
            if(o&  1)   F->Function(&node2->children[0],node1);
            if(o&  2)   F->Function(&node2->children[1],node1);
            if(o&  4)   F->Function(&node2->children[2],node1);
            if(o&  8)   F->Function(&node2->children[3],node1);
            if(o& 16)   F->Function(&node2->children[4],node1);
            if(o& 32)   F->Function(&node2->children[5],node1);
            if(o& 64)   F->Function(&node2->children[6],node1);
            if(o&128)   F->Function(&node2->children[7],node1);
        }
        if(node2->depth()<depth){
            if(o&  1) if(node2->children[0].children)
                    __ProcessMaxDepthNodeAdjacentNodes(dx1,dy1,dz1,node1,radius1,
                                                   &node2->children[0],radius,cWidth,
                                                       depth,F);
            if(o&  2) if(node2->children[1].children)
                    __ProcessMaxDepthNodeAdjacentNodes(dx2,dy1,dz1,node1,radius1,
                                                       &node2->children[1],radius,cWidth,
                                                       depth,F);
            if(o&  4) if(node2->children[2].children)
                    __ProcessMaxDepthNodeAdjacentNodes(dx1,dy2,dz1,node1,radius1,
                                                       &node2->children[2],radius,cWidth,
                                                       depth,F);
            if(o&  8) if(node2->children[3].children)
                    __ProcessMaxDepthNodeAdjacentNodes(dx2,dy2,dz1,node1,radius1,
                                                       &node2->children[3],radius,cWidth,
                                                       depth,F);
            if(o& 16) if(node2->children[4].children)
                    __ProcessMaxDepthNodeAdjacentNodes(dx1,dy1,dz2,node1,radius1,
                                                       &node2->children[4],radius,cWidth,
                                                       depth,F);
            if(o& 32) if(node2->children[5].children)
                    __ProcessMaxDepthNodeAdjacentNodes(dx2,dy1,dz2,node1,radius1,
                                                       &node2->children[5],radius,cWidth,
                                                       depth,F);
            if(o& 64) if(node2->children[6].children)
                    __ProcessMaxDepthNodeAdjacentNodes(dx1,dy2,dz2,node1,radius1,
                                                       &node2->children[6],radius,cWidth,
                                                       depth,F);
            if(o&128) if(node2->children[7].children)
                    __ProcessMaxDepthNodeAdjacentNodes(dx2,dy2,dz2,node1,radius1,
                                                       &node2->children[7],radius,cWidth,
                                                       depth,F);
        }
    }
}


template<class NodeAdjacencyFunction>
void OctNode::__ProcessTerminatingNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                    OctNode* node1, const float& radius1,
                                                    OctNode* node2, const float& radius2, const float& cWidth2,
                                                    NodeAdjacencyFunction* F)
{
    /**     $cWidth takes half of node2's radius because it's convenient to
     *      move the center position to generate new (dx, dy, dz).
     *      Another half of node2's radius is offset by $radius.    */
    float cWidth=cWidth2/2;
    float radius=radius2/2;

    /**     (dx, dy, dz) is the vector(center1->center2)
     *      radius1:    width of node1 * 1.500001
     *      radius:     radius of node2 / 2 = width of node2 / 4
     *      cWidth:     radius of node2 / 2 = width of node2 / 4     */
    int o=ChildOverlap(dx,dy,dz,radius1+radius,cWidth);
    if(o){
        float dx1=dx-cWidth;
        float dx2=dx+cWidth;
        float dy1=dy-cWidth;
        float dy2=dy+cWidth;
        float dz1=dz-cWidth;
        float dz2=dz+cWidth;
        if(o&  1){
            if(F->Function(&node2->children[0],node1) && node2->children[0].children )
                __ProcessTerminatingNodeAdjacentNodes(dx1,dy1,dz1,node1,radius1,
                                                      &node2->children[0],radius,
                                                      cWidth,F);
        }
        if(o&  2){
            if(F->Function(&node2->children[1],node1) && node2->children[1].children )
                __ProcessTerminatingNodeAdjacentNodes(dx2,dy1,dz1,node1,radius1,
                                                      &node2->children[1],radius,
                                                      cWidth,F);
        }
        if(o&  4){
            if(F->Function(&node2->children[2],node1) && node2->children[2].children )
                __ProcessTerminatingNodeAdjacentNodes(dx1,dy2,dz1,node1,radius1,
                                                      &node2->children[2],radius,
                                                      cWidth,F);
        }
        if(o&  8){
            if(F->Function(&node2->children[3],node1) && node2->children[3].children )
                __ProcessTerminatingNodeAdjacentNodes(dx2,dy2,dz1,node1,radius1,
                                                      &node2->children[3],radius,
                                                      cWidth,F);
        }
        if(o& 16){
            if(F->Function(&node2->children[4],node1) && node2->children[4].children )
                __ProcessTerminatingNodeAdjacentNodes(dx1,dy1,dz2,node1,radius1,
                                                      &node2->children[4],radius,
                                                      cWidth,F);
        }
        if(o& 32){
            if(F->Function(&node2->children[5],node1) && node2->children[5].children )
                __ProcessTerminatingNodeAdjacentNodes(dx2,dy1,dz2,node1,radius1,
                                                      &node2->children[5],radius,
                                                      cWidth,F);
        }
        if(o& 64){
            if(F->Function(&node2->children[6],node1) && node2->children[6].children )
                __ProcessTerminatingNodeAdjacentNodes(dx1,dy2,dz2,node1,radius1,
                                                      &node2->children[6],radius,
                                                      cWidth,F);
        }
        if(o&128){
            if(F->Function(&node2->children[7],node1) && node2->children[7].children )
                __ProcessTerminatingNodeAdjacentNodes(dx2,dy2,dz2,node1,radius1,
                                                      &node2->children[7],radius,
                                                      cWidth,F);
        }
    }
}

template<class PointAdjacencyFunction>
void OctNode::__ProcessPointAdjacentNodes(const float& dx,const float& dy,const float& dz,
                                          OctNode* node2,const float& radius2,const float& cWidth2,
                                          PointAdjacencyFunction* F)
{
    /**     $cWidth takes half of node2's radius because it's convenient to
      *      move the center position to generate new (dx, dy, dz).
      *      Another half of node2's radius is offset by $radius.    */
    float cWidth=cWidth2/2;
    float radius=radius2/2;

    int o=ChildOverlap(dx,dy,dz,radius,cWidth);
    if(o){
        float dx1=dx-cWidth;
        float dx2=dx+cWidth;
        float dy1=dy-cWidth;
        float dy2=dy+cWidth;
        float dz1=dz-cWidth;
        float dz2=dz+cWidth;
        if(o&  1){
            F->Function(&node2->children[0]);
            if(node2->children[0].children)
                __ProcessPointAdjacentNodes(dx1,dy1,dz1,
                                            &node2->children[0],radius,
                                            cWidth,F);
        }
        if(o&  2){
            F->Function(&node2->children[1]);
            if(node2->children[1].children)
                __ProcessPointAdjacentNodes(dx2,dy1,dz1,
                                            &node2->children[1],radius,
                                            cWidth,F);
        }
        if(o&  4){
            F->Function(&node2->children[2]);
            if(node2->children[2].children)
                __ProcessPointAdjacentNodes(dx1,dy2,dz1,
                                            &node2->children[2],radius,
                                            cWidth,F);
        }
        if(o&  8){
            F->Function(&node2->children[3]);
            if(node2->children[3].children)
                __ProcessPointAdjacentNodes(dx2,dy2,dz1,
                                            &node2->children[3],radius,
                                            cWidth,F);
        }
        if(o& 16){
            F->Function(&node2->children[4]);
            if(node2->children[4].children)
                __ProcessPointAdjacentNodes(dx1,dy1,dz2,
                                            &node2->children[4],radius,
                                            cWidth,F);
        }
        if(o& 32){
            F->Function(&node2->children[5]);
            if(node2->children[5].children)
                __ProcessPointAdjacentNodes(dx2,dy1,dz2,
                                            &node2->children[5],radius,
                                            cWidth,F);
        }
        if(o& 64){
            F->Function(&node2->children[6]);
            if(node2->children[6].children)
                __ProcessPointAdjacentNodes(dx1,dy2,dz2,
                                            &node2->children[6],radius,
                                            cWidth,F);
        }
        if(o&128){
            F->Function(&node2->children[7]);
            if(node2->children[7].children)
                __ProcessPointAdjacentNodes(dx2,dy2,dz2,
                                            &node2->children[7],radius,
                                            cWidth,F);
        }
    }
}

template<class PointAdjacencyFunction>
void OctNode::__ProcessPointAdjacentNodes(const float& dx,const float& dy,const float& dz,
                                          const float& radius1,
                                          OctNode* node2,const float& radius2,const float& cWidth2,
                                          PointAdjacencyFunction* F)
{
    /**     $cWidth takes half of node2's radius because it's convenient to
     *      move the center position to generate new (dx, dy, dz).
     *      Another half of node2's radius is offset by $radius.    */
    float cWidth=cWidth2/2;
    float radius=radius2/2;

    int o=ChildOverlap(dx,dy,dz,radius1+radius,cWidth);
    if(o){
        float dx1=dx-cWidth;
        float dx2=dx+cWidth;
        float dy1=dy-cWidth;
        float dy2=dy+cWidth;
        float dz1=dz-cWidth;
        float dz2=dz+cWidth;
        if(o&  1){
            F->Function(&node2->children[0]);
            if(node2->children[0].children)
                __ProcessPointAdjacentNodes(dx1,dy1,dz1,
                                            radius1,
                                            &node2->children[0],radius,
                                            cWidth,F);
        }
        if(o&  2){
            F->Function(&node2->children[1]);
            if(node2->children[1].children)
                __ProcessPointAdjacentNodes(dx2,dy1,dz1,
                                            radius1,
                                            &node2->children[1],radius,
                                            cWidth,F);
        }
        if(o&  4){
            F->Function(&node2->children[2]);
            if(node2->children[2].children)
                __ProcessPointAdjacentNodes(dx1,dy2,dz1,
                                            radius1,
                                            &node2->children[2],radius,
                                            cWidth,F);
        }
        if(o&  8){
            F->Function(&node2->children[3]);
            if(node2->children[3].children)
                __ProcessPointAdjacentNodes(dx2,dy2,dz1,
                                            radius1,
                                            &node2->children[3],radius,
                                            cWidth,F);
        }
        if(o& 16){
            F->Function(&node2->children[4]);
            if(node2->children[4].children)
                __ProcessPointAdjacentNodes(dx1,dy1,dz2,
                                            radius1,
                                            &node2->children[4],radius,
                                            cWidth,F);
        }
        if(o& 32){
            F->Function(&node2->children[5]);
            if(node2->children[5].children)
                __ProcessPointAdjacentNodes(dx2,dy1,dz2,
                                            radius1,
                                            &node2->children[5],radius,
                                            cWidth,F);
        }
        if(o& 64){
            F->Function(&node2->children[6]);
            if(node2->children[6].children)
                __ProcessPointAdjacentNodes(dx1,dy2,dz2,
                                            radius1,
                                            &node2->children[6],radius,
                                            cWidth,F);
        }
        if(o&128){
            F->Function(&node2->children[7]);
            if(node2->children[7].children)
                __ProcessPointAdjacentNodes(dx2,dy2,dz2,
                                            radius1,
                                            &node2->children[7],radius,
                                            cWidth,F);
        }
    }
}

const OctNode* OctNode::__faceNeighbor(const int& dir,const int& off) const{
    if(!parent) return NULL;    // there is no neighbor outside this node
    int pIndex=int(this-parent->children);  //get its children index
    pIndex^=(1<<dir);    // reverse the bit in direction bit to get neighbor's right index in its sibling
    if( (pIndex & (1<<dir) ) == (off<<dir) )
        return &parent->children[pIndex];
    const OctNode* temp=parent->__faceNeighbor(dir,off);
    if(!temp || !temp->children)
        return temp;
    return &temp->children[pIndex];
}

OctNode* OctNode::__faceNeighbor(const int&dir, const int& off,const int& forceChildren) {
    if(!parent) return NULL;
    int pIndex=int(this-parent->children);
    pIndex^=(1<<dir);
    if( (pIndex & (1<<dir) ) == (off<<dir) )
        return &parent->children[pIndex];
    OctNode* temp=parent->__faceNeighbor(dir,off,forceChildren);
    if(!temp) return NULL;
    if(!temp->children){
        if(forceChildren) temp->initChildren();
        else return temp;
    }
    return &temp->children[pIndex];
}


const OctNode* OctNode::__edgeNeighbor(const int& o,const int i[2],const int idx[2]) const {
    if(!parent) return NULL;
    int pIndex=int(this-parent->children);
    int aIndex,x[DIMENSION];

    Cube::FactorCornerIndex(pIndex,x[0],x[1],x[2]);
    /**     the first bit of aIndex denotes x coords ^ */
    aIndex=( ~( (i[0] ^ x[idx[0]]) | ( (i[1] ^ x[idx[1]] ) << 1 ) )  ) & 3;
    pIndex^=(7 ^ (1<<o));

    /**     get the edge neighbor from the parent's face adjacent neighbor  */
    if(aIndex==1){
        const OctNode* temp=parent->__faceNeighbor(idx[0],i[0]);
        if(!temp || !temp->children) return NULL;
        else return &temp->children[pIndex];
    }
    if(aIndex==2){
        const OctNode* temp=parent->__faceNeighbor(idx[1],i[1]);
        if(!temp || !temp->children) return NULL;
        else return &temp->children[pIndex];
    }
    /**     get the edge neighbor from the parent   */
    if(aIndex==0){
        return &parent->children[pIndex];
    }
    /**     get the edge neighbor from the parent's edge adjacent neighbor  */
    if(aIndex==3){
        const OctNode* temp= parent->__edgeNeighbor(o,i,idx);
        if(!temp || !temp->children) return temp;
        return &temp->children[pIndex];
    }
    return NULL;
}
OctNode* OctNode::__edgeNeighbor(const int& o,const int i[2],const int idx[2],const int& forceChildren) {
    if(!parent) return NULL;
    int pIndex=int(this-parent->children);
    int aIndex,x[DIMENSION];

    Cube::FactorCornerIndex(pIndex,x[0],x[1],x[2]);
    aIndex=( ~( (i[0] ^ x[idx[0]]) | ( (i[1] ^ x[idx[1]] ) << 1 ) )  ) & 3;
    pIndex^=(7 ^ (1<<o));

    /**     get the edge neighbor from the parent's face adjacent neighbor  */
    if(aIndex==1){
        OctNode* temp=parent->__faceNeighbor(idx[0],i[0],0);
        if(!temp || !temp->children) return NULL;
        else return &temp->children[pIndex];
    }
    if(aIndex==2){
        OctNode* temp=parent->__faceNeighbor(idx[1],i[1],0);
        if(!temp || !temp->children) return NULL;
        else return &temp->children[pIndex];
    }
    /**     get the edge neighbor from the parent   */
    if(aIndex==0){
        return &parent->children[pIndex];
    }
    /**     get the edge neighbor from the parent's edge adjacent neighbor  */
    if(aIndex==3){
        OctNode* temp= parent->__edgeNeighbor(o,i,idx,forceChildren);
        if(!temp) return NULL;
        if(!temp->children) {
            if(forceChildren) temp->initChildren();
            else return temp;
        }
        return &temp->children[pIndex];
    }
    return NULL;

}

inline bool OctNode::Overlap(const Point3D<float>& c1, const Point3D<float>& c2,const float& dWidth) {
    if(fabs(c1.coords[0]-c2.coords[0])>=dWidth || fabs(c1.coords[1]-c2.coords[1])>=dWidth || fabs(c1.coords[2]-c2.coords[2])>=dWidth)
        return false;
    return true;
}
inline bool OctNode::Overlap(const float& c1, const float& c2, const float& c3, const float& dWidth) {
    if(c1>=dWidth || c1<=-dWidth || c2>=dWidth || c2<=-dWidth || c3>=dWidth || c3<=-dWidth)
        return false;
    return true;
}


inline int OctNode::ChildOverlap(const float& dx, const float& dy, const float& dz, const float& d, const float& cRadius2){
    float w1=d-cRadius2;
    float w2=d+cRadius2;
    int overlap=0;

    int test=0, test1=0;
    /**     -d <(dx-cRadius2)<d     */
    if(dx < w2 && dx > -w1) test =1;
    /**     -d <(dx+cRadius2)<d     */
    if(dx < w1 && dx > -w2) test|=2;

    if(!test) return 0;
    /**     -d <(dz-cRadius2)<d     */
    if(dz < w2 && dz > -w1) test1 =test;
    /**     -d <(dz+cRadius2)<d     */
    if(dz < w1 && dz > -w2) test1|=test<<4;

    if(!test1) return 0;
    /**     -d <(dy-cRadius2)<d     */
    if(dy < w2 && dy > -w1) overlap =test1;
    /**     -d <(dy+cRadius2)<d     */
    if(dy < w1 && dy > -w2) overlap|=test1<<2;
    return overlap;
}

// OctNode public member
const int OctNode::DepthShift=5;
const int OctNode::OffsetShift=19;
const int OctNode::OffsetShift1=OffsetShift;
const int OctNode::OffsetShift2=OffsetShift+OffsetShift1;
const int OctNode::OffsetShift3=OffsetShift+OffsetShift2;
const int OctNode::DepthMask=(1<<DepthShift)-1;
const int OctNode::OffsetMask=(1<<OffsetShift)-1;

inline void OctNode::DepthAndOffset(const long long& index, int& depth,int offset[3]){
    depth=int(index & DepthMask);
    offset[0]=(int((index>>OffsetShift1)&OffsetMask)+1)&(~(1<<depth));
    offset[1]=(int((index>>OffsetShift2)&OffsetMask)+1)&(~(1<<depth));
    offset[2]=(int((index>>OffsetShift3)&OffsetMask)+1)&(~(1<<depth));
}

inline void OctNode::CenterAndWidth(const long long &index, Point3D<float> &center, float &width) {
    int depth,offset[3];
    DepthAndOffset(index, depth, offset);
    width=float(1.0/(1<<depth));
    for(int dim=0;dim<DIMENSION;++dim){
        center.coords[dim]=float(0.5+offset[dim])*width;
    }
}

inline int OctNode::Depth(const long long& index){
    return int(index & DepthMask);
}

OctNode::OctNode() {
    parent=children=NULL;
    d=off[0]=off[1]=off[2]=0;
}

OctNode::~OctNode() {
    if(children) delete [] children;
    parent=children=NULL;
}

bool OctNode::initChildren() {
    int i,j,k;
    if(children) delete [] children;
    children=NULL;
    children=new OctNode[Cube::CORNERS];

    if(!children){
        std::cout<<"Failed to init children in OctNode::initChildren\n";
        exit(0);
    }

    int d1,off1[3];
    depthAndOffset(d1,off1);
    for(int i=0;i<2;++i){
        for(int j=0;j<2;++j){
            for(int k=0;k<2;++k){
                int idx=Cube::CornerIndex(i,j,k);
                children[idx].parent=this;
                children[idx].children=NULL;
                int off2[3];
                off2[0]=(off1[0]<<1)+i;
                off2[1]=(off1[1]<<1)+j;
                off2[2]=(off1[2]<<1)+k;
                Index(d1+1,off2,children[idx].d,children[idx].off);
            }
        }
    }
    return true;
}

inline int OctNode::depth(void) const {
    return d;
}

inline Point3D<float> OctNode::center(void) const {
    Point3D<float> center;
    int depth,offset[3];
    depth=d;
    for(int i=0;i<DIMENSION;++i){
        offset[i]=(off[i]+1)&(~(1<<depth));
    }
    float width=float(1.0/(1<<depth));
    for(int dim=0;dim<DIMENSION;++dim){
        center.coords[dim]=float(0.5+offset[dim])*width;
    }
    return center;
}

inline float OctNode::width(void) const {
    return float(1.0/(1<<d));
}

inline void OctNode::depthAndOffset(int &depth, int offset[3]) const {
    depth=int(d);
    for(int i=0;i<DIMENSION;++i){
        offset[i]=(int(off[i])+1)&(~(1<<depth));
    }
}

inline void OctNode::centerAndWidth(Point3D<float>& center, float& width) const{
    int depth,offset[3];
    depth=d;
    for(int i=0;i<DIMENSION;++i){
        offset[i]=(off[i]+1)&(~(1<<depth));
    }
    width=float(1.0/(1<<depth));
    for(int dim=0;dim<DIMENSION;++dim){
        center.coords[dim]=float(0.5+offset[dim])*width;
    }
}

inline void OctNode::Index(const int& depth, const int offset[3], int& d, int off[3]){
    d=depth;
    for(int i=0;i<DIMENSION;++i){
        off[i]=(1<<depth)+offset[i]-1;
    }
}

int OctNode::leaves() const {
    if(!children) return 1;
    int c=0;
    for(int i=0;i<Cube::CORNERS;++i){
        c+=children[i].leaves();
    }
    return c;
}

int OctNode::maxDepthLeaves(const int& maxDepth) const {
    if(depth()>maxDepth) return 0;
    if(!children) return 1;
    int c=0;
    for(int i=0;i<Cube::CORNERS;++i){
        c+=children[i].maxDepthLeaves(maxDepth);
    }
    return c;
}

int OctNode::nodes(void) const {
    if(!children) return 1;
    int c=0;
    for(int i=0;i<Cube::CORNERS;++i){
        c+=children[i].nodes();
    }
    return c+1;
}

int OctNode::maxDepth(void) const {
    if(!children) return 0;
    int c=0,d;
    for(int i=0;i<Cube::CORNERS;++i) {
        d=children[i].maxDepth();
        if(!i || d>c) c=d;
    }
    return c+1;
}

const OctNode* OctNode::nextLeaf(const OctNode* current) const {
    if(!current){
        const OctNode* temp=this;
        while(temp->children) temp=temp->children;
        return temp;
    }
    if(current->children) return current->nextLeaf();
    const OctNode* temp=nextBranch(current);
    if(!temp) return NULL;
    return temp->nextLeaf();
}

OctNode* OctNode::nextLeaf(OctNode* current) {
    if(!current){
        OctNode* temp=this;
        while(temp->children) temp=temp->children;
        return temp;
    }
    if(current->children) return current->nextLeaf();
    OctNode* temp=nextBranch(current);
    if(!temp) return NULL;
    return temp->nextLeaf();
}

const OctNode* OctNode::nextNode(const OctNode* currentNode) const{
    if(!currentNode) return this;
    else if(currentNode->children) return currentNode->children;
    return nextBranch(currentNode);
}

OctNode* OctNode::nextNode(OctNode* currentNode){
    if(!currentNode) return this;
    else if(currentNode->children) return currentNode->children;
    return nextBranch(currentNode);
}
const OctNode* OctNode::nextBranch(const OctNode* current) const {
    if(!current->parent || current == this) return NULL;
    if(current-current->parent->children==Cube::CORNERS-1) return nextBranch(current->parent);
    return current+1;
}

OctNode* OctNode::nextBranch(OctNode* current) {
    if(!current->parent || current == this) return NULL;
    if(current-current->parent->children==Cube::CORNERS-1) return nextBranch(current->parent);
    return current+1;
}

void OctNode::setFullDepth(const int& maxDepth) {
    if(maxDepth){
        if(!children) this->initChildren();
        for(int i=0;i<8;++i){
            this->children[i].setFullDepth(maxDepth-1);
        }
    }
}

template<class NodeAdjacencyFunction>
void OctNode::processNodeNodes(OctNode* node,NodeAdjacencyFunction* F,const int& processCurrent){
    if(processCurrent){F->Function(this,node);}
    if(!children){return;}
    __processNodeNodes(node,F);
}

const OctNode* OctNode::root(void) const {
    const OctNode* temp=this;
    while(temp->parent) {temp=temp->parent;}
    return temp;
}

bool OctNode::write(const char* fileName) const {
    FILE* fp=fopen(fileName,"wb");
    if(!fp) return false;
    bool ret=write(fp);
    fclose(fp);
    return ret;
}

bool OctNode::write(FILE* fp) const {
    fwrite(this,sizeof(OctNode),1,fp);
    if(children){
        for(int i=0;i<Cube::CORNERS;++i){
            children[i].write(fp);
        }
    }
    return true;
}

bool OctNode::read(const char* fileName) {
    FILE* fp=fopen(fileName,"rb");
    if(!fp) return false;
    bool ret=read(fp);
    fclose(fp);
    return ret;
}

bool OctNode::read(FILE* fp) {
    fread(this,sizeof(OctNode),1,fp);
    parent=NULL;
    if(children){
        children=NULL;
        initChildren();
        for(int i=0;i<Cube::CORNERS;++i){
            children[i].read(fp);
            children[i].parent=this;
        }
    }
    return true;
}

OctNode& OctNode::operator = (const OctNode& node){
    int i;
    if(children) delete [] children;
    children=NULL;

    d=node.depth();
    for(int i=0;i<DIMENSION;++i) this->off[i]=node.off[i];
    if(node.children){
        initChildren();
        for(int i=0;i<Cube::CORNERS;++i){
            this->children[i]=node.children[i];
        }
    }
    return *this;
}


template<class NodeAdjacencyFunction>
void OctNode::ProcessNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                     OctNode* node1, const float& radius1,
                                     OctNode* node2, const float& radius2, const float& width2,
                                     NodeAdjacencyFunction* F,const int& processCurrent)
{
    /**     SetLaplacianWeights(): DivergenceFunction
     *      (dx, dy, dz) is the vector(node2->node1)
     *      radius1:    1.5 * width of node1
     *      radius2:    1.5 * width of node2
     *      width2:     width of node2              */

    /**     SolveFixedDepthMatrix(): LaplacianProjectionFunction
     *      radius1:    1.4449  * width
     *      radius2:    1.0001  * width
     *      width2:     width                        */
    if(!Overlap(dx,dy,dz,radius1+radius2)) return;
    if(processCurrent) F->Function(node2,node1);
    if(!node2->children) return;
    /**     SetLaplacianWeights(): DivergenceFunction
     *      (dx, dy, dz) is the vector(node1->node2)
     *      radius1:    1.5 * width of node1
     *      radius2:    1.5 * width of node2
     *      cWidth2:    radius of node2 = width of node2 / 2    */

    /**     SolveFixedDepthMatrix(): LaplacianProjectionFunction
     *      radius1:    1.4449  * width
     *      radius2:    1.0001  * width
     *      width2:     radius = width / 2                      */
    __ProcessNodeAdjacentNodes(-dx,-dy,-dz,node1,radius1,node2,radius2,width2/2,F);
}

template<class NodeAdjacencyFunction>
void OctNode::ProcessNodeAdjacentNodes(OctNode* node1, const float& radius1,
                                       OctNode* node2, const float& radius2,
                                       NodeAdjacencyFunction* F,const int& processCurrent)
{
    /**     SetLaplacianWeights(): DivergenceFunction
     *      radius1:    1.5
     *      radius2:    1.5                         */
    Point3D<float> c1,c2;
    float w1,w2;

    node1->centerAndWidth(c1,w1);

    /**     c2 is the &tree's center,
     *      w2 is the width of [0, 1]^3 node = 1    */
    node2->centerAndWidth(c2,w2);

    /**     SetLaplacianWeights(): DivergenceFunction
     *      (dx, dy, dz) is the vector(node2->node1)
     *      radius1:    1.5 * width of node1
     *      radius2:    1.5 * width of node2
     *      width2:     width of node2              */
    ProcessNodeAdjacentNodes(c1.coords[0]-c2.coords[0],
                             c1.coords[1]-c2.coords[1],
                             c1.coords[2]-c2.coords[2],
                             node1,radius1*w1,
                             node2,radius2*w2,w2,
                             F,processCurrent);
}

template<class NodeAdjacencyFunction>
void OctNode::ProcessFixedDepthNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                 OctNode* node1, const float& radius1,
                                                 OctNode* node2, const float& radius2, const float& width2,
                                                 const int& depth,
                                                 NodeAdjacencyFunction* F,
                                                 const int& processCurrent)
{
    int d=node2->depth();
    if(d>depth) return;
    if(!Overlap(dx,dy,dz,radius1+radius2)) return;
    if(d==depth) {
        if(processCurrent)
            F->Function(node2,node1);
    }else{
        if(!node2->children) return;
        __ProcessFixedDepthNodeAdjacentNodes(-dx,-dy,-dz,
                                             node1,radius1,
                                             node2,radius2,width2/2,
                                             depth-1,F);
    }
}

template<class NodeAdjacencyFunction>
void OctNode::ProcessFixedDepthNodeAdjacentNodes(OctNode* node1, const float& radius1,
                                                 OctNode* node2, const float& radius2,
                                                 const int& depth,
                                                 NodeAdjacencyFunction* F,
                                                 const int& processCurrent)
{
    Point3D<float> c1,c2;
    float w1,w2;
    node1->centerAndWidth(c1,w1);
    node2->centerAndWidth(c2,w2);
    ProcessFixedDepthNodeAdjacentNodes(c1.coords[0]-c2.coords[0],
                                       c1.coords[1]-c2.coords[1],
                                       c1.coords[2]-c2.coords[2],
                                       node1,radius1*w1,
                                       node2,radius2*w2,w2,
                                       depth,F,processCurrent);
}

template<class NodeAdjacencyFunction>
void OctNode::ProcessMaxDepthNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                               OctNode* node1, const float& radius1,
                                               OctNode* node2, const float& radius2, const float& width2,
                                               const int& depth,
                                               NodeAdjacencyFunction* F,
                                               const int& processCurrent)
{
    int d=node2->depth();
    if(d > depth) return;

    /**     radius1:    width of node1 * 3
     *      radius2:    width of tree node2 * 0.5 = radius of node2
     *      width2:     width of tree node2
     *      If c1 and c2 is not close enough, don't process.    */
    if(!Overlap(dx,dy,dz,radius1+radius2))
    /**      This return almost won't be triggered               */
        return;


    /**     if node2 is not deep enough,
     *      then call initChildren()        */
    if(processCurrent) F->Function(node2,node1);


    /**     (-dx, -dy, -dz) is the vector(center1->center2)
     *      radius1:    width of node1 * 3
     *      radius2:    width of tree node2 * 0.5 = radius of node2
     *      cWidth2:    radius of node2 = width of node2 / 2    */
    if(d < depth && node2->children) __ProcessMaxDepthNodeAdjacentNodes(-dx,-dy,-dz,node1,radius1,node2,radius2,width2/2,depth-1,F);
    /**     Because above function execute on node2's children, so depth need to minus 1    */
}

template<class NodeAdjacencyFunction>
void OctNode::ProcessMaxDepthNodeAdjacentNodes(OctNode* node1, const float& radius1,
                                               OctNode* node2, const float& radius2,
                                               const int& depth,
                                               NodeAdjacencyFunction* F,
                                               const int& processCurrent)
{
    /**     radius1:    3.0
     *      radius2:    0.5     */

    Point3D<float> c1,c2;
    float w1,w2;
    node1->centerAndWidth(c1,w1);

    /**     c2 is the &tree's center,
     *      w2 is the width of [0, 1]^3 node = 1     */
    node2->centerAndWidth(c2,w2);

    /**     (dx, dy, dz) is the vector(center2->center1)
     *      radius1:    width of node1 * 3
     *      radius2:    width of tree node2 * 0.5 = radius of node2 */
    ProcessMaxDepthNodeAdjacentNodes(c1.coords[0]-c2.coords[0],
                                     c1.coords[1]-c2.coords[1],
                                     c1.coords[2]-c2.coords[2],
                                     node1,radius1*w1,
                                     node2,radius2*w2,w2,
                                     depth,F,processCurrent);
}

template<class TerminatingNodeAdjacencyFunction>
void OctNode::ProcessTerminatingNodeAdjacentNodes(const float& dx, const float& dy, const float& dz,
                                                  OctNode* node1, const float& radius1,
                                                  OctNode* node2, const float& radius2, const float& width2,
                                                  TerminatingNodeAdjacencyFunction* F,
                                                  const int& processCurrent)
{
    /**     (dx, dy, dz) is the vector(center2->center1)
     *      radius1:    width of node1 * 1.500001
     *      radius2:    width of tree node2 * 0.5 = radius of node2
     *      width2:     width of tree node2                     */
    if(!Overlap(dx,dy,dz,radius1+radius2)) return;
    if(processCurrent) F->Function(node2,node1);
    if(!node2->children) return;

    /**     (-dx, -dy, -dz) is the vector(center1->center2)
     *      radius1:    width of node1 * 1.500001
     *      radius2:    width of tree node2 * 0.5 = radius of node2
     *      cWidth2:    radius of node2 = width of node2 / 2    */
    __ProcessTerminatingNodeAdjacentNodes(-dx,-dy,-dz,node1,radius1,node2,radius2,width2/2,F);
}

template<class TerminatingNodeAdjacencyFunction>
void OctNode::ProcessTerminatingNodeAdjacentNodes(OctNode* node1, const float& radius1,
                                                  OctNode* node2, const float& radius2,
                                                  TerminatingNodeAdjacencyFunction* F,
                                                  const int& processCurrent)
{
    /**     GetFixedDepthLaplacian():
     *      radius1:    1.500001
     *      radius2:    0.5         */
    Point3D<float> c1,c2;
    float w1,w2;
    node1->centerAndWidth(c1,w1);
    node2->centerAndWidth(c2,w2);

    /**     (dx, dy, dz) is the vector(center2->center1)
     *      radius1:    width of node1 * 1.500001
     *      radius2:    width of tree node2 * 0.5 = radius of node2 */
    ProcessTerminatingNodeAdjacentNodes(c1.coords[0]-c2.coords[0],
                                        c1.coords[1]-c2.coords[1],
                                        c1.coords[2]-c2.coords[2],
                                        node1,radius1*w1,
                                        node2,radius2*w2,w2,
                                        F,processCurrent);
}


template<class PointAdjacencyFunction>
void OctNode::ProcessPointAdjacentNodes(const float& dx,const float& dy,const float& dz,
                                        OctNode* node2,const float& radius2,const float& width2,
                                        PointAdjacencyFunction* F,const int& processCurrent)
{
    if(!Overlap(dx,dy,dz,radius2)) return;
    if(processCurrent) F->Function(node2);
    if(!node2->children) return;
    __ProcessPointAdjacentNodes(-dx,-dy,-dz,node2,radius2,width2/2,F);
}


template<class PointAdjacencyFunction>
void OctNode::ProcessPointAdjacentNodes(const Point3D<float>& center1,
                                        OctNode* node2, const float& radius2,
                                        PointAdjacencyFunction* F, const int& processCurrent)
{
    Point3D<float> c2;
    float w2;
    node2->centerAndWidth(c2,w2);
    ProcessPointAdjacentNodes(center1.coords[0]-c2.coords[0],
                              center1.coords[1]-c2.coords[1],
                              center1.coords[2]-c2.coords[2],
                              node2,radius2*w2,w2,
                              F,processCurrent);
}

template<class PointAdjacencyFunction>
void OctNode::ProcessPointAdjacentNodes(const float& dx,const float& dy,const float& dz,
                                        const float& radius1,
                                        OctNode* node2,const float& radius2,const float& width2,
                                        PointAdjacencyFunction* F,const int& processCurrent)
{
    if(!Overlap(dx,dy,dz,radius1+radius2)) return;
    if(processCurrent) F->Function(node2);
    if(!node2->children) return;
    __ProcessPointAdjacentNodes(-dx,-dy,-dz,radius1,node2,radius2,width2/2,F);
}

template<class PointAdjacencyFunction>
void OctNode::ProcessPointAdjacentNodes(const Point3D<float>& center1, const float& radius1,
                                        OctNode* node2, const float& radius2,
                                        PointAdjacencyFunction* F, const int& processCurrent)
{
    Point3D<float> c2;
    float w2;
    node2->centerAndWidth(c2,w2);
    ProcessPointAdjacentNodes(center1.coords[0]-c2.coords[0],
                              center1.coords[1]-c2.coords[1],
                              center1.coords[2]-c2.coords[2],
                              radius1,
                              node2,radius2*w2,w2,
                              F,processCurrent);

}

int OctNode::CornerIndex(const Point3D<float>& center,const Point3D<float>& p) {
    int cIndex=0;
    if(p.coords[0]>center.coords[0]) cIndex|=1;
    if(p.coords[1]>center.coords[1]) cIndex|=2;
    if(p.coords[2]>center.coords[2]) cIndex|=4;
    return cIndex;
}

const OctNode* OctNode::faceNeighbor(const int& faceIndex) const{
    return __faceNeighbor(faceIndex>>1,faceIndex&1);
}

OctNode* OctNode::faceNeighbor(const int& faceIndex, const int& forceChildren) {
    return __faceNeighbor(faceIndex>>1,faceIndex&1,forceChildren);
}

const OctNode* OctNode::edgeNeighbor(const int& edgeIndex) const{
    int idx[2],o,i[2];
    Cube::FactorEdgeIndex(edgeIndex,o,i[0],i[1]);
    switch (o) {
        case 0: idx[0]=1;   idx[1]=2;   break;
        case 1: idx[0]=0;   idx[1]=2;   break;
        case 2: idx[0]=0;   idx[1]=1;   break;
    }
    return __edgeNeighbor(o,i,idx);
}

OctNode* OctNode::edgeNeighbor(const int& edgeIndex,const int& forceChildren) {
    int idx[2],o,i[2];
    Cube::FactorEdgeIndex(edgeIndex,o,i[0],i[1]);
    switch (o) {
        case 0: idx[0]=1;   idx[1]=2;   break;
        case 1: idx[0]=0;   idx[1]=2;   break;
        case 2: idx[0]=0;   idx[1]=1;   break;
    }
    return __edgeNeighbor(o,i,idx,forceChildren);
}


const OctNode* OctNode::cornerNeighbor(const int& cornerIndex) const {
    int pIndex,aIndex=0;
    if(!parent){return NULL;}

    pIndex=int(this-parent->children);
    aIndex=(cornerIndex ^ pIndex);	// The disagreement bits
    pIndex=(~pIndex)&7;				// The antipodal point
    if(aIndex==7){					// Agree on no bits
        return &parent->children[pIndex];
    }
    else if(aIndex==0){				// Agree on all bits
        const OctNode* temp=((const OctNode*)parent)->cornerNeighbor(cornerIndex);
        if(!temp || !temp->children){return temp;}
        else{return &temp->children[pIndex];}
    }
    else if(aIndex==6){				// Agree on face 0
        const OctNode* temp=((const OctNode*)parent)->__faceNeighbor(0,cornerIndex & 1);
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==5){				// Agree on face 1
        const OctNode* temp=((const OctNode*)parent)->__faceNeighbor(1,(cornerIndex & 2)>>1);
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==3){				// Agree on face 2
        const OctNode* temp=((const OctNode*)parent)->__faceNeighbor(2,(cornerIndex & 4)>>2);
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==4){				// Agree on edge 2
        const OctNode* temp=((const OctNode*)parent)->edgeNeighbor(8 | (cornerIndex & 1) | (cornerIndex & 2) );
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==2){				// Agree on edge 1
        const OctNode* temp=((const OctNode*)parent)->edgeNeighbor(4 | (cornerIndex & 1) | ((cornerIndex & 4)>>1) );
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==1){				// Agree on edge 0
        const OctNode* temp=((const OctNode*)parent)->edgeNeighbor(((cornerIndex & 2) | (cornerIndex & 4))>>1 );
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    return NULL;
}

OctNode* OctNode::cornerNeighbor(const int& cornerIndex,const int& forceChildren) {
    int pIndex,aIndex=0;
    if(!parent){return NULL;}

    pIndex=int(this-parent->children);
    aIndex=(cornerIndex ^ pIndex);	// The disagreement bits
    pIndex=(~pIndex)&7;				// The antipodal point
    if(aIndex==7){					// Agree on no bits
        return &parent->children[pIndex];
    }
    else if(aIndex==0){				// Agree on all bits
        OctNode* temp=((OctNode*)parent)->cornerNeighbor(cornerIndex,forceChildren);
        if(!temp){return NULL;}
        if(!temp->children){
            if(forceChildren){temp->initChildren();}
            else{return temp;}
        }
        return &temp->children[pIndex];
    }
    else if(aIndex==6){				// Agree on face 0
        OctNode* temp=((OctNode*)parent)->__faceNeighbor(0,cornerIndex & 1,0);
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==5){				// Agree on face 1
        OctNode* temp=((OctNode*)parent)->__faceNeighbor(1,(cornerIndex & 2)>>1,0);
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==3){				// Agree on face 2
        OctNode* temp=((OctNode*)parent)->__faceNeighbor(2,(cornerIndex & 4)>>2,0);
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==4){				// Agree on edge 2
        OctNode* temp=((OctNode*)parent)->edgeNeighbor(8 | (cornerIndex & 1) | (cornerIndex & 2) );
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==2){				// Agree on edge 1
        OctNode* temp=((OctNode*)parent)->edgeNeighbor(4 | (cornerIndex & 1) | ((cornerIndex & 4)>>1) );
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    else if(aIndex==1){				// Agree on edge 0
        OctNode* temp=((OctNode*)parent)->edgeNeighbor(((cornerIndex & 2) | (cornerIndex & 4))>>1 );
        if(!temp || !temp->children){return NULL;}
        else{return & temp->children[pIndex];}
    }
    return NULL;
}


const OctNode* OctNode::getNearestLeaf(const Point3D<float>& p) const {
    int nearest;
    float temp,dist2;
    if(!children) return this;
    for(int i=0;i<Cube::CORNERS;++i){
        temp=SquareDistance(children[i].center(),p);
        if(!i || temp<dist2){
            dist2=temp;
            nearest=i;
        }
    }
    return children[nearest].getNearestLeaf(p);
}

OctNode* OctNode::getNearestLeaf(const Point3D<float>& p) {
    Point3D<float> center;
    float width;
    int cIndex;
    if(!children) return this;
    centerAndWidth(center,width);
    OctNode* temp=this;
    while(temp->children){
        cIndex=CornerIndex(center,p);
        temp=&temp->children[cIndex];
        width/=2;
        if(cIndex&1) center.coords[0]+=width/2;
        else         center.coords[0]-=width/2;
        if(cIndex&2) center.coords[1]+=width/2;
        else         center.coords[1]-=width/2;
        if(cIndex&4) center.coords[2]+=width/2;
        else         center.coords[2]-=width/2;
    }
    return temp;
}

bool OctNode::CommonEdge(const OctNode* node1, const int& eIndex1, const OctNode* node2, const int& eIndex2) {
    int o1,o2,i1,i2,j1,j2;

    Cube::FactorEdgeIndex(eIndex1,o1,i1,j1);
    Cube::FactorEdgeIndex(eIndex2,o2,i2,j2);
    if(o1!=o2) return 0;

    int dir[2];
    switch (o1) {
        case 0: dir[0]=1;   dir[1]=2;   break;
        case 1: dir[0]=0;   dir[1]=2;   break;
        case 2: dir[0]=0;   dir[1]=1;   break;
    };
    int d1,d2,off1[3],off2[3];
    node1->depthAndOffset(d1,off1);
    node2->depthAndOffset(d2,off2);

    /**     conclude from the calculating process of OctNode::off[3]    */
    int idx1[2],idx2[2];
    idx1[0]=off1[dir[0]]+(1<<d1)+i1;
    idx1[1]=off1[dir[1]]+(1<<d1)+j1;

    idx2[0]=off2[dir[0]]+(1<<d2)+i2;
    idx2[1]=off2[dir[1]]+(1<<d2)+j2;

    if(d1>d2){
        idx2[0]<<=(d1-d2);
        idx2[1]<<=(d1-d2);
    }else{
        idx1[0]<<=(d2-d1);
        idx1[1]<<=(d2-d1);
    }
    if(idx1[0]==idx2[0] && idx1[1]==idx2[1])
        return true;
    return false;
}

int OctNode::CompareForwardPointerDepths(const void* v1,const void* v2){
    const OctNode *n1,*n2;
    n1=(*(const OctNode**)v1);
    n2=(*(const OctNode**)v2);
    if(n1->d!=n2->d){return int(n1->d)-int(n2->d);}
    while(n1->parent != n2->parent){
        n1=n1->parent;
        n2=n2->parent;
    }
    if(n1->off[0]!=n2->off[0]){return int(n1->off[0])-int(n2->off[0]);}
    if(n1->off[1]!=n2->off[1]){return int(n1->off[1])-int(n2->off[1]);}
    return int(n1->off[2])-int(n2->off[2]);
}

int OctNode::CompareBackwardPointerDepths(const void* v1,const void* v2){
    const OctNode *n1,*n2;
    n1=(*(const OctNode**)v1);
    n2=(*(const OctNode**)v2);
    if(n1->d!=n2->d){return -(int(n1->d)-int(n2->d));}
    while(n1->parent != n2->parent){
        n1=n1->parent;
        n2=n2->parent;
    }
    if(n1->off[0]!=n2->off[0]){return -(int(n1->off[0])-int(n2->off[0]));}
    if(n1->off[1]!=n2->off[1]){return -(int(n1->off[1])-int(n2->off[1]));}
    return -(int(n1->off[2])-int(n2->off[2]));
//    return (*(const OctNode**)v2)->depth()-(*(const OctNode**)v1)->depth();
}

inline int OctNode::Overlap2(const int& depth1, const int offSet1[DIMENSION], const float& multiplier1,
                             const int& depth2, const int offSet2[DIMENSION], const float& multiplier2)
{
    /**     depth1:     depth-startingDepth
     *      depth2:     depth
     *      offSet1:    node1's offset in GetFixedDepthLaplacian()
     *      off1:       the node close to AdjacencySetFunction's node
     *      multiplier1:    0.5
     *      multiplier2:    2.001                                     */
    int d=depth2-depth1;    //  $startingDepth in GetFixedDepthLaplacian()
    float w=multiplier2+multiplier1*(1<<d);
    float w2=float( (1<<(d-1)) - 0.5 );
    if(fabs(float(offSet2[0]-(offSet1[0]<<d) ) - w2 ) >= w ||
       fabs(float(offSet2[1]-(offSet1[1]<<d) ) - w2 ) >= w ||
       fabs(float(offSet2[2]-(offSet1[2]<<d) ) - w2 ) >= w)
        return 0;
    return 1;
}

// OctNode::Neighbors
OctNode::Neighbors::Neighbors(void) {
    clear();
}

void OctNode::Neighbors::clear(void) {
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) for(int k=0;k<3;++k) neighbors[i][j][k]=NULL;
}


// OctNode::NeighborKey
OctNode::NeighborKey::NeighborKey(void) {
    neighbors=NULL;
}

OctNode::NeighborKey::~NeighborKey(void) {
    if(neighbors) delete [] neighbors;
    neighbors=NULL;
}

void OctNode::NeighborKey::set(const int& d){
    if(neighbors) delete [] neighbors;
    neighbors=NULL;
    if(d<0) return;
    neighbors = new Neighbors[d+1];
}

OctNode::Neighbors& OctNode::NeighborKey::setNeighbors(OctNode *node) {
    int d=node->depth();
    /**      if $node the center of neighbors[d],
     *       means that neighbors[d] has been set to be the neighbors of $node  */

    if(node!=neighbors[d].neighbors[1][1][1]){
        neighbors[d].clear();

        if(!node->parent) neighbors[d].neighbors[1][1][1]=node;
        else{
            int i,j,k,x1,y1,z1,x2,y2,z2;
            int idx=int(node-node->parent->children);
            Cube::FactorCornerIndex(idx, x1, y1, z1);
            Cube::FactorCornerIndex((~idx)&7, x2, y2, z2);
            for(i=0;i<2;++i){
                for(j=0;j<2;++j){
                    for(k=0;k<2;++k){
                        /**     if i=x1, j=y1, k=z1,
                          *     x2+i=1, y2+j=1, z2+k=1
                          *     set neighbors[1][1][1] to be itself */
                        neighbors[d].neighbors[x2+i][y2+j][z2+k]=&node->parent->children[Cube::CornerIndex(i,j,k)];
                    }
                }
            }
            /**     temp contains parent's neighbors    */
            Neighbors& temp=setNeighbors(node->parent);

            // set the neighbors from across the faces of its parent
            /**     if x1=0, then x2=1,
             *      x=1,2 will be set,
             *      but x=0 is still empty.
             *      if x1=1, then x2=0,
             *      x=0,1 will be set,
             *      but x=2 is still empty.
             */
            i=x1<<1;
            if(temp.neighbors[i][1][1]){
                if(!temp.neighbors[i][1][1]->children) temp.neighbors[i][1][1]->initChildren();
                for(j=0;j<2;++j){
                    for(k=0;k<2;++k){
                        neighbors[d].neighbors[i][y2+j][z2+k]=&temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];
                    }
                }
            }

            j=y1<<1;
            if(temp.neighbors[1][j][1]){
                if(!temp.neighbors[1][j][1]->children) temp.neighbors[1][j][1]->initChildren();
                for(i=0;i<2;++i){
                    for(k=0;k<2;++k){
                        neighbors[d].neighbors[x2+i][j][z2+k]=&temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];
                    }
                }
            }

            k=z1<<1;
            if(temp.neighbors[1][1][k]){
                if(!temp.neighbors[1][1][k]->children) temp.neighbors[1][1][k]->initChildren();
                for(i=0;i<2;++i){
                    for(j=0;j<2;++j){
                        neighbors[d].neighbors[x2+i][y2+j][k]=&temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];
                    }
                }
            }

            // Set the neighbors from across the edges
            i=x1<<1;    j=y1<<1;
            if(temp.neighbors[i][j][1]){
                if(!temp.neighbors[i][j][1]->children) temp.neighbors[i][j][1]->initChildren();
                for(k=0;k<2;++k)
                    neighbors[d].neighbors[i][j][z2+k]=&temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];
            }
            i=x1<<1;    k=z1<<1;
            if(temp.neighbors[i][1][k]){
                if(!temp.neighbors[i][1][k]->children) temp.neighbors[i][1][k]->initChildren();
                for(j=0;j<2;++j)
                    neighbors[d].neighbors[i][y2+j][k]=&temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];
            }
            j=y1<<1;    k=z1<<1;
            if(temp.neighbors[1][j][k]){
                if(!temp.neighbors[1][j][k]->children) temp.neighbors[1][j][k]->initChildren();
                for(i=0;i<2;++i)
                    neighbors[d].neighbors[x2+i][j][k]=&temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];
            }

            // Set the neighbor from across the corner
            i=x1<<1;    j=y1<<1;    k=z1<<1;
            if(temp.neighbors[i][j][k]){
                if(!temp.neighbors[i][j][k]->children) temp.neighbors[i][j][k]->initChildren();
                neighbors[d].neighbors[i][j][k]=&temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
            }
        }
    }
    return neighbors[d];
}

OctNode::Neighbors& OctNode::NeighborKey::getNeighbors(OctNode *node) {
    int d=node->depth();
    if(node!=neighbors[d].neighbors[1][1][1]){
        neighbors[d].clear();

        if(!node->parent){neighbors[d].neighbors[1][1][1]=node;}
        else{
            int i,j,k,x1,y1,z1,x2,y2,z2;
            int idx=int(node-node->parent->children);
            Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
            Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);
            for(i=0;i<2;i++){
                for(j=0;j<2;j++){
                    for(k=0;k<2;k++){
                        neighbors[d].neighbors[x2+i][y2+j][z2+k]=&node->parent->children[Cube::CornerIndex(i,j,k)];
                    }
                }
            }
            Neighbors& temp=getNeighbors(node->parent);

            // Set the neighbors from across the faces
            i=x1<<1;
            if(temp.neighbors[i][1][1] && temp.neighbors[i][1][1]->children){
                for(j=0;j<2;j++){for(k=0;k<2;k++){neighbors[d].neighbors[i][y2+j][z2+k]=&temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];}}
            }
            j=y1<<1;
            if(temp.neighbors[1][j][1] && temp.neighbors[1][j][1]->children){
                for(i=0;i<2;i++){for(k=0;k<2;k++){neighbors[d].neighbors[x2+i][j][z2+k]=&temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];}}
            }
            k=z1<<1;
            if(temp.neighbors[1][1][k] && temp.neighbors[1][1][k]->children){
                for(i=0;i<2;i++){for(j=0;j<2;j++){neighbors[d].neighbors[x2+i][y2+j][k]=&temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];}}
            }

            // Set the neighbors from across the edges
            i=x1<<1;	j=y1<<1;
            if(temp.neighbors[i][j][1] && temp.neighbors[i][j][1]->children){
                for(k=0;k<2;k++){neighbors[d].neighbors[i][j][z2+k]=&temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];}
            }
            i=x1<<1;	k=z1<<1;
            if(temp.neighbors[i][1][k] && temp.neighbors[i][1][k]->children){
                for(j=0;j<2;j++){neighbors[d].neighbors[i][y2+j][k]=&temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];}
            }
            j=y1<<1;	k=z1<<1;
            if(temp.neighbors[1][j][k] && temp.neighbors[1][j][k]->children){
                for(i=0;i<2;i++){neighbors[d].neighbors[x2+i][j][k]=&temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];}
            }

            // Set the neighbor from across the corner
            i=x1<<1;	j=y1<<1;	k=z1<<1;
            if(temp.neighbors[i][j][k] && temp.neighbors[i][j][k]->children){
                neighbors[d].neighbors[i][j][k]=&temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
            }
        }
    }
    return neighbors[node->depth()];
}