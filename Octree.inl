

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
    return (eSegmentCount >> (faceIndex*2) & 3);
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
    } else {
        isoNode=IsoNodeData();
    }
    value=0;
}

NodeData::~NodeData() {}

// OctNode private member

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

int OctNode::initChildren() {
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
                Index(d+1,off2,children[idx].d,children[idx].off);
            }
        }
    }
    return 1;
}

inline int OctNode::depth(void) const {
    return d;
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

const OctNode* OctNode::root(void) const {
    const OctNode* temp=this;
    while(temp->parent) {temp=temp->parent;}
    return temp;
}

int OctNode::write(const char* fileName) const {
    FILE* fp=fopen(fileName,"wb");
    if(!fp) return 0;
    int ret=write(fp);
    fclose(fp);
    return ret;
}

int OctNode::write(FILE* fp) const {
    fwrite(this,sizeof(OctNode),1,fp);
    if(children){
        for(int i=0;i<Cube::CORNERS;++i){
            children[i].write(fp);
        }
    }
    return 1;
}

int OctNode::read(const char* fileName) {
    FILE* fp=fopen(fileName,"rb");
    if(!fp) return 0;
    int ret=read(fp);
    fclose(fp);
    return ret;
}

int OctNode::read(FILE* fp) {
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
    return 1;
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
                          *     set neighbors[1][1][1] to be itself
                          */
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