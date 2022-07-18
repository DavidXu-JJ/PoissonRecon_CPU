
#include "Octree.h"

#define FORCE_UNIT_NORMALS 1
#define USE_DOT_RATIOS 1

const float EPSILON=float(1e-6);
const float ROUND_EPS=float(1e-5);

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

// Octree
template<int Degree>
float Octree<Degree>::GetLaplacian(const int idx[DIMENSION]) const
{
#if USE_DOT_RATIOS
    return Real(fData.dotTable[idx[0]]*fData.dotTable[idx[1]]*fData.dotTable[idx[2]]*(fData.d2DotTable[idx[0]]+fData.d2DotTable[idx[1]]+fData.d2DotTable[idx[2]]));
#else // !USE_DOT_RATIOS
    double dot[3];
    dot[0]=fData.dotTable[idx[0]];
    dot[1]=fData.dotTable[idx[1]];
    dot[2]=fData.dotTable[idx[2]];
    return Real(
            fData.d2DotTable[idx[0]]*dot[1]*dot[2]+
            fData.d2DotTable[idx[1]]*dot[0]*dot[2]+
            fData.d2DotTable[idx[2]]*dot[0]*dot[1]
    );
#endif // USE_DOT_RATIOS
}

template <int Degree>
float Octree<Degree>::GetDivergence(const int idx[DIMENSION],const Point3D<float>& normal) const
{
#if USE_DOT_RATIOS
    double dot=fData.dotTable[idx[0]]*fData.dotTable[idx[1]]*fData.dotTable[idx[2]];
    return Real(dot*(fData.dDotTable[idx[0]]*normal.coords[0]+fData.dDotTable[idx[1]]*normal.coords[1]+fData.dDotTable[idx[2]]*normal.coords[2]));
#else // !USE_DOT_RATIOS
    double dot[DIMENSION];
	dot[0]=fData.dotTable[idx[0]];
	dot[1]=fData.dotTable[idx[1]];
	dot[2]=fData.dotTable[idx[2]];
	return Real(
		fData.dDotTable[idx[0]]*normal.coords[0]*dot[1]*dot[2]+
		fData.dDotTable[idx[1]]*normal.coords[1]*dot[0]*dot[2]+
		fData.dDotTable[idx[2]]*normal.coords[2]*dot[0]*dot[1]
	);
#endif // !USE_DOT_RATIOS
}

template<int Degree>
void Octree<Degree>::DivergenceFunction::Function(OctNode *node1, OctNode *node2) {
    for(int i=0;i<3;++i)
        scratch[i]=index[i]+int(node1->off[i]);
    node1->nodeData.value+=ot->GetDivergence(scratch,normal);
}

template<int Degree>
void Octree<Degree>::LaplacianProjectionFunction::Function(OctNode *node1, OctNode *node2) {
    for(int i=0;i<3;++i)
        scratch[i]=index[i]+int(node1->off[i]);
    node1->nodeData.value-=float(ot->GetLaplacian(scratch)*value);
}

/*
template<int Degree>
int Octree<Degree>::LaplacianMatrixFunction::Function(OctNode *node1, OctNode *node2) {
    float temp;
    int d1=int(node1->d);
    int x1,y1,z1;
    x1=int(node1->off[0]);
    y1=int(node1->off[1]);
    z1=int(node1->off[2]);

    int dDepth=d2-d1;
    int dcoords;
    dcoords=(x2>>dDepth)-x1;

}
*/

template<int Degree>
int Octree<Degree>::NonLinearUpdateWeightContribution(OctNode* node, const Point3D<float>& position) {
    int i,j,k;
    /**     Get the neighbor of $node,
     *      initChildren is called and all elements is valid.   */
    OctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);

    /**     dx[0] denote dim x,
     *      dx[1] denote dim y,
     *      dx[2] denote dim z.                 */
    double x,dx[DIMENSION][3];
    double width;
    /**     Get the center and width of $node   */
    Point3D<float> center;
    float w;
    node->centerAndWidth(center,w);
    width=w;

    for(i=0;i<DIMENSION;++i){
        /**     x denote the distance to index 0 in this dimension,
         *      x bigger, $position is more close to neighbor with index 0 in this dimension */
        x=(center.coords[i]-position.coords[i]-width)/width;
        dx[i][0]=1.125+1.500*x+0.500*x*x;

        /**     x denote the distance to index 1 in this dimension,
         *      abs(x) bigger, less contribution to neighbor with index 1 in this dimension */
        x=(center.coords[i]-position.coords[i])/width;
        dx[i][1]=0.750        -      x*x;

        dx[i][2]=1.0-dx[i][0]-dx[i][1];
    }

    for(i=0;i<3;++i){
        for(j=0;j<3;++j){
            for(k=0;k<3;++k){
                if(neighbors.neighbors[i][j][k])
                    neighbors.neighbors[i][j][k]->nodeData.centerWeightContribution+=float(dx[0][i]*dx[1][j]*dx[2][k]);
            }
        }
    }
    return 0;
}


template<int Degree>
float Octree<Degree>::NonLinearGetSampleWeight(OctNode* node, const Point3D<float>& position) {
    float weight=0;
    double x,dx[DIMENSION][3];
    int i,j,k;
    OctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);
    double width;
    Point3D<float> center;
    float w;
    node->centerAndWidth(center,w);
    width=w;

    for(i=0;i<DIMENSION;++i){
        /**     x denote the distance to index 0 in this dimension,
         *      x bigger, $position is more close to neighbor with index 0 in this dimension */
        x=(center.coords[i]-position.coords[i]-width)/width;
        dx[i][0]=1.125+1.500*x+0.500*x*x;

        /**     x denote the distance to index 1 in this dimension,
         *      abs(x) bigger, less contribution to neighbor with index 1 in this dimension */
        x=(center.coords[i]-position.coords[i])/width;
        dx[i][1]=0.750        -      x*x;

        dx[i][2]=1.0-dx[i][0]-dx[i][1];
    }

    for(i=0;i<DIMENSION;++i){
        for(j=0;j<DIMENSION;++j){
            for(k=0;k<DIMENSION;++k){
                if(neighbors.neighbors[i][j][k])
                    weight+=float(dx[0][i]*dx[1][j]*dx[2][k]*
                                  neighbors.neighbors[i][j][k]->nodeData.centerWeightContribution);
            }
        }
    }
    return float(1.0/weight);
}

template<int Degree>
void Octree<Degree>::NonLinearGetSampleDepthAndWeight(OctNode* node, const Point3D<float>& position,
                                                      const float& samplesPerNode,
                                                      float& depth,
                                                      float& weight)
{
    OctNode* temp=node;
    weight=float(1.0)/NonLinearGetSampleWeight(temp,position);
    if(weight>=samplesPerNode+1)
        depth=float(temp->depth()+log( weight / (samplesPerNode+1) ) / log(double(1<<(DIMENSION-1) ) ) );
    else{
        float oldAlpha,newAlpha;
        oldAlpha=newAlpha=weight;
        while(newAlpha<(samplesPerNode+1) && temp->parent){
            temp=temp->parent;
            oldAlpha=newAlpha;
            newAlpha=float(1.0)/NonLinearGetSampleWeight(temp,position);
        }
        depth=float(temp->depth()+log( newAlpha / (samplesPerNode+1) ) / log( newAlpha / oldAlpha ) );
    }
    weight=float(pow(double(1<<(DIMENSION-1) ), -double(depth) ) );
}

template<int Degree>
int Octree<Degree>::NonLinearSplatOrientedPoint(OctNode* node, const Point3D<float>& position, const Point3D<float>& normal) {
    double x,dxdydz,dx[DIMENSION][3];
    int i,j,k;
    OctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);
    double width;
    Point3D<float> center;
    float w;
    node->centerAndWidth(center,w);
    width=w;

    for(i=0;i<DIMENSION;++i){
        /**     x denote the distance to index 0 in this dimension,
         *      x bigger, $position is more close to neighbor with index 0 in this dimension */
        x=(center.coords[i]-position.coords[i]-width)/width;
        dx[i][0]=1.125+1.500*x+0.500*x*x;

        /**     x denote the distance to index 1 in this dimension,
         *      abs(x) bigger, less contribution to neighbor with index 1 in this dimension */
        x=(center.coords[i]-position.coords[i])/width;
        dx[i][1]=0.750        -      x*x;

        dx[i][2]=1.0-dx[i][0]-dx[i][1];
    }

    for(i=0;i<DIMENSION;++i){
        for(j=0;j<DIMENSION;++j){
            for(k=0;k<DIMENSION;++k){
                if(neighbors.neighbors[i][j][k]){
                    dxdydz=dx[0][i]*dx[1][j]*dx[2][k];
                    int idx=neighbors.neighbors[i][j][k]->nodeData.nodeIndex;
                    if(idx<0){
                        Point3D<float> n;
                        n.coords[0]=n.coords[1]=n.coords[2]=0;
                        idx=neighbors.neighbors[i][j][k]->nodeData.nodeIndex=int(normals->size());
                        normals->push_back(n);
                    }
                    (*normals)[idx].coords[0]+=float(normal.coords[0]*dxdydz);
                    (*normals)[idx].coords[1]+=float(normal.coords[1]*dxdydz);
                    (*normals)[idx].coords[2]+=float(normal.coords[2]*dxdydz);
                }
            }
        }
    }
    return 0;
}

template<int Degree>
void Octree<Degree>::NonLinearSplatOrientedPoint(const Point3D<float>& position,
                                                 const Point3D<float>& normal,
                                                 const int& splatDepth,
                                                 const float& samplesPerNode,
                                                 const int& minDepth,
                                                 const int& maxDepth)
{
    double dx;
    Point3D<float> n;
    OctNode* temp;
    int i;
    double width;
    Point3D<float> myCenter;
    float myWidth;
    myCenter.coords[0]=myCenter.coords[1]=myCenter.coords[2]=float(0.5);
    myWidth=float(1.0);

    temp=&tree;
    /**     reach the nearest node to $position point at splatDepth */
    while(temp->depth()<splatDepth){
        if(!temp->children){
            printf("Can't reach splatDepth.\nOctree<Degree>::NonLinearSplatOrientedPoint error.\n");
            return;
        }
        int cIndex=OctNode::CornerIndex(myCenter,position);
        temp=&temp->children[cIndex];
        myWidth/=2;
        for(i=0;i<DIMENSION;++i){
            if(cIndex & (1<<DIMENSION)) myCenter.coords[i]+=myWidth/2;
            else                        myCenter.coords[i]-=myWidth/2;
        }
    }

    float alpha,newDepth;
    NonLinearGetSampleDepthAndWeight(temp,position,samplesPerNode,newDepth,alpha);

    if(newDepth<minDepth)
        newDepth=float(minDepth);
    if(newDepth>maxDepth)
        newDepth=float(maxDepth);
    int topDepth=int(ceil(newDepth));

    dx=1.0-(topDepth-newDepth);
    if(topDepth<=minDepth){
        topDepth=minDepth;
        dx=1;
    }
    else if(topDepth>maxDepth){
        topDepth=maxDepth;
        dx=1;
    }

    /**     make $temp be the nearest node at $topDepth to $position point  */
    while(temp->depth()>topDepth) temp=temp->parent;
    while(temp->depth()<topDepth) {
        if(!temp->children) temp->initChildren();
        int cIndex=OctNode::CornerIndex(myCenter,position);
        temp=&temp->children[cIndex];
        myWidth/=2;
        for(i=0;i<DIMENSION;++i){
            if(cIndex & (1<<DIMENSION)) myCenter.coords[i]+=myWidth/2;
            else                        myCenter.coords[i]-=myWidth/2;
        }
    }

    /**     scale the original normal   */
    width=1.0/(1<<temp->depth());
    for(i=0;i<DIMENSION;++i) {
        n.coords[i] = normal.coords[i] * alpha
                      / float(pow(width, 3)
                      * float(dx));
    }

    /**     update the contribution to normals[temp] at $topDepth   */
    NonLinearSplatOrientedPoint(temp,position,n);
    /**     In case that newDepth is between [topDepth-1, topDepth],
     *      update the contribution to normals[temp->parent] at $topDepth-1 */
    if(fabs(1.0-dx)>EPSILON){
        dx=float(1.0/dx);
        temp=temp->parent;
        width=1.0/(1<<temp->depth());

        for(i=0;i<DIMENSION;++i) {
            n.coords[i] = normal.coords[i] * alpha
                          / float(pow(width,3))
                          * float(dx);
        }
        NonLinearSplatOrientedPoint(temp,position,n);
    }
}

template<int Degree>
int Octree<Degree>::HasNormals(OctNode *node, const float &epsilon) {
    int hasNormals=0;
    int& idx=node->nodeData.nodeIndex;
    Point3D<float>& normal=(*normals)[idx];
    if(idx >= 0 && Length(normal) > epsilon)
        hasNormals=1;
    if(node->children)
        for(int i=0;i<Cube::CORNERS && !hasNormals; ++i){
            hasNormals|=HasNormals(&node->children[i],epsilon);
        }
    return hasNormals;
}

// Octree public function
template<int Degree>
Octree<Degree>::Octree(void) {
    radius=0;
    postNormalSmooth=0;
}

template<int Degree>
void Octree<Degree>::setFunctionData(const PPolynomial<Degree>& ReconstructionFunction,
                     const int& maxDepth,
                     const int& normalize,
                     const float& normalSmooth)
{
    radius=float(fabs(ReconstructionFunction.polys[0].start));
    if(normalSmooth>0) postNormalSmooth=normalSmooth;
    fData.set(maxDepth,ReconstructionFunction,normalize,USE_DOT_RATIOS);
}

template<int Degree>
int Octree<Degree>::setTree(char* fileName,
                            const int& maxDepth,
                            const int& binary,
                            const int& kernelDepth,
                            const float& samplesPerNode,
                            const float& scaleFactor,
                            Point3D<float>& center,
                            float& scale,
                            const int& resetSamples)
{
    FILE* fp;
    Point3D<float> min,max,position,normal,myCenter;
    float myWidth;
    int i,cnt=0;
    float c[2*DIMENSION];
    OctNode* temp;
    int splatDepth=0;

    NodeData::UseIndex=1;
    neighborKey.set(maxDepth);
    splatDepth=std::max(splatDepth,kernelDepth);
    if(binary) fp=fopen(fileName,"rb");
    else fp=fopen(fileName,"r");
    if(!fp) return 0;


    /**     Get the min and max in each dimension   */
    while(1){
        if(binary){
            if(fread(c,sizeof(float),2*DIMENSION,fp)!=2*DIMENSION)
                break;
        } else{
            if(fscanf(fp," %f %f %f %f %f %f ",&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])!=2*DIMENSION)
                break;
        }
        for(i=0;i<DIMENSION;++i){
            if(!cnt || c[i]<min.coords[i]) min.coords[i]=c[i];
            if(!cnt || c[i]>max.coords[i]) max.coords[i]=c[i];
        }
        ++cnt;
    }

    /**     center and scale used for normalization is ready     */
    for(i=0;i<DIMENSION;++i){
        if(!i || scale<max.coords[i]-min.coords[i]) scale=float(max.coords[i]-min.coords[i]);
        center.coords[i]=float(max.coords[i]+min.coords[i])/2;
    }
    scale*=scaleFactor;
    for(i=0;i<DIMENSION;++i)
        center.coords[i]-=scale/2;


    /**     Initialize weight contribution of the node in this tree     */
    if(splatDepth>0){
        /**     set position indicator  */
        fseek(fp,0,SEEK_SET);
        cnt=0;
        while(1){
            if(binary){
                if(fread(c,sizeof(float),2*DIMENSION,fp)!=2*DIMENSION)
                    break;
            }else{
                if(fscanf(fp," %f %f %f %f %f %f ",&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])!=2*DIMENSION)
                    break;
            }
            /**     project point to [0,1]^3    */
            for(i=0;i<DIMENSION;++i)
                position.coords[i]=(c[i]-center.coords[i])/scale;

            /**     center at (0.5, 0.5, 0.5)   */
            myCenter.coords[0]=myCenter.coords[1]=myCenter.coords[2]=float(0.5);
            myWidth=float(1.0);

            /**     kick out the point which coords is smaller than 0 or bigger than 1  */
            for(i=0;i<DIMENSION;++i)
                if(position.coords[i]<myCenter.coords[i]-myWidth/2 || position.coords[i]>myCenter.coords[i]+myWidth/2)
                    break;
            if(i!=DIMENSION) continue;

            /**     Now, the point is fine to use.
             *      Update the $position contribution to every neighbor node with depth < splatDepth.
             *      myCenter and myWidth use to track the trace to its subnode. */
            temp=&tree;
            int d=0;
            while(d<splatDepth){
                NonLinearUpdateWeightContribution(temp,position);
                if(!temp->children) temp->initChildren();
                /**     cIndex = cornerIndex close to position  */
                int cIndex=OctNode::CornerIndex(myCenter,position);
                /**     step to the nearest subnode,
                 *      adjust the myCenter and myWidth     */
                temp=&temp->children[cIndex];
                myWidth/=2;
                for(i=0;i<DIMENSION;++i) {
                    if (cIndex & (1<<i)) myCenter.coords[i] += myWidth / 2;
                    else                 myCenter.coords[i] -= myWidth / 2;
                }
                ++d;
            }
            NonLinearUpdateWeightContribution(temp,position);
            ++cnt;
        }
    }

    normals=new std::vector<Point3D<float> >();
    fseek(fp,0,SEEK_SET);
    cnt=0;

    /**     Update normals in each related node */
    while(1){
        if(binary){
            if(fread(c,sizeof(float),2*DIMENSION,fp)!=2*DIMENSION)
                break;
        }else{
            if(fscanf(fp," %f %f %f %f %f %f ",&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])!=2*DIMENSION)
                break;
        }
        for(i=0;i<DIMENSION;++i){
            position.coords[i]=(c[i]-center.coords[i])/scale;
            normal.coords[i]=c[DIMENSION+i];
        }
        myCenter.coords[0]=myCenter.coords[1]=myCenter.coords[2]=float(0.5);
        myWidth=float(1.0);

        for(i=0;i<DIMENSION;++i)
            if(position.coords[i]<myCenter.coords[i]-myWidth/2 || position.coords[i]>myCenter.coords[i]+myWidth/2)
                break;
        if(i!=DIMENSION) continue;

#if FORCE_UNIT_NORMALS
        float len=float(Length(normal));
        if(len>EPSILON)
            len=1.0f/len;
        len*=(2<<maxDepth);
        for(i=0;i<DIMENSION;++i){
            normal.coords[i]*=len;
        }
#endif  // FORCE_UNIT_NORMALS

        /**     Now, position and normal has been normalized    */

        if(resetSamples && samplesPerNode>0 && splatDepth){
            /**     Update the normals of the node at depth defined by $samplesPerNode
             *      and nodes' weight contribution initialized  */
            NonLinearUpdateWeightContribution(position,normal,splatDepth,samplesPerNode,1,maxDepth);
        }else{

            /**     calculate the weight $alpha */
            float alpha=1;
            temp=&tree;
            if(splatDepth){
                int d=0;
                /**     reach the splatDepth node   */
                while(d<splatDepth){
                    int cIndex=OctNode::CornerIndex(myCenter,position);
                    temp=&temp->children[cIndex];
                    myWidth/=2;
                    for(i=0;i<DIMENSION;++i){
                        if(cIndex & (1<<DIMENSION)) myCenter.coords[i]+=myWidth/2;
                        else                        myCenter.coords[i]-=myWidth/2;
                    }
                    d++;
                }
                alpha=NonLinearGetSampleWeight(temp,position);
            }

            /**     scale normal    */
            for(i=0;i<DIMENSION;++i)
                normal.coords[i]*=alpha;

            /**     Update the normals of the node only at maxDepth  */
            int d=0;
            while(d<maxDepth){
                if(!temp->children) temp->initChildren();
                int cIndex=OctNode::CornerIndex(myCenter,position);
                temp=&temp->children[cIndex];
                myWidth/=2;
                for(i=0;i<DIMENSION;++i){
                    if(cIndex & (1<<DIMENSION)) myCenter.coords[i]+=myWidth/2;
                    else                        myCenter.coords[i]-=myWidth/2;
                }
                d++;
            }
            NonLinearSplatOrientedPoint(temp,position,normal);
        }
    }
    fclose(fp);
    return cnt;
}

template<int Degree>
void Octree<Degree>::ClipTree(void) {
    OctNode* temp;
    temp=tree.nextNode();
    while(temp){
        if(temp->children){
            int hasNormals=0;
            for(int i=0;i<Cube::CORNERS && !hasNormals;++i){
                hasNormals=HasNormals(&temp->children[i],EPSILON);
            }
            if(!hasNormals)
                temp->children=NULL;
        }
        temp=tree.nextNode(temp);
    }
}

