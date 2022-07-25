#include "PPolynomial.h"
#include "Geometry.h"
#include "Octree.h"
#include <cmath>
#include "PlyFile.h"
#include "MultiGridOctreeData.h"



int main(int argc,char* argv[])
{
    char fileName[]="/Users/davidxu/horse.npts";
    char outName[]="/Users/davidxu/horse.ply";

    Octree<2> tree;
    PPolynomial<2> ReconstructionFunction=PPolynomial<2>::GaussianApproximation();
    int depth=10,kernelDepth=8;
    Point3D<float> center;
    float scale=1.0f;

    for(int i=0;i<3;++i)
        center.coords[i]=0;

    tree.setFunctionData(ReconstructionFunction,10,0,float(1.0)/(1<<depth));

    tree.setTree(fileName,depth,kernelDepth,1.0,1.25,center,scale,1);

    SortedTreeNodes s;

    tree.ClipTree();

    tree.finalize1(3);

    tree.SetLaplacianWeights();

    tree.finalize2(3);

    tree.LaplacianMatrixIteration(0);

    float isoValue=tree.GetIsoValue();


    CoredVectorMeshData mesh;

    tree.GetMCIsoTriangles(isoValue,&mesh,0);

    PlyWriteTriangles(outName,&mesh,PLY_BINARY_NATIVE,center,scale,NULL,0);


}

