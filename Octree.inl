

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

// OctNode

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