#ifndef AVRO_LIB_GRAPHICS_BSP_H_
#define AVRO_LIB_GRAPHICS_BSP_H_

namespace avro
{

template <typename type> class Topology;
class Simplex;

namespace graphics
{

class BSPTriangle {
public:
  BSPTriangle( const Topology<Simplex>& topology , index_t k );



};

class BSPTree {

public:
  BSPTree( const Topology<Simplex>& triangles );

  void build();

  void add_triangle();

private:
  const Topology<Simplex>& triangles_;

  std::shared_ptr< BSPTree > front_;
  std::shared_ptr< BSPTree > back_;
};

/*

static void gl2psBuildBspTree(GL2PSbsptree *tree, GL2PSlist *primitives)
{
  GL2PSprimitive *prim, *frontprim = NULL, *backprim = NULL;
  GL2PSlist *frontlist, *backlist;
  GLint i, idx;

  tree->front = NULL;
  tree->back = NULL;
  tree->primitives = gl2psListCreate(1, 2, sizeof(GL2PSprimitive*));
  idx = gl2psFindRoot(primitives, &prim);
  gl2psGetPlane(prim, tree->plane);
  gl2psAddPrimitiveInList(prim, tree->primitives);

  frontlist = gl2psListCreate(1, 2, sizeof(GL2PSprimitive*));
  backlist = gl2psListCreate(1, 2, sizeof(GL2PSprimitive*));

  for(i = 0; i < gl2psListNbr(primitives); i++){
    if(i != idx){
      prim = *(GL2PSprimitive**)gl2psListPointer(primitives,i);
      switch(gl2psSplitPrimitive(prim, tree->plane, &frontprim, &backprim)){
      case GL2PS_COINCIDENT:
        gl2psAddPrimitiveInList(prim, tree->primitives);
        break;
      case GL2PS_IN_BACK_OF:
        gl2psAddPrimitiveInList(prim, backlist);
        break;
      case GL2PS_IN_FRONT_OF:
        gl2psAddPrimitiveInList(prim, frontlist);
        break;
      case GL2PS_SPANNING:
        gl2psAddPrimitiveInList(backprim, backlist);
        gl2psAddPrimitiveInList(frontprim, frontlist);
        gl2psFreePrimitive(&prim);
        break;
      }
    }
  }

  if(gl2psListNbr(tree->primitives)){
    gl2psListSort(tree->primitives, gl2psTrianglesFirst);
  }

  if(gl2psListNbr(frontlist)){
    gl2psListSort(frontlist, gl2psTrianglesFirst);
    tree->front = (GL2PSbsptree*)gl2psMalloc(sizeof(GL2PSbsptree));
    gl2psBuildBspTree(tree->front, frontlist);
  }
  else{
    gl2psListDelete(frontlist);
  }

  if(gl2psListNbr(backlist)){
    gl2psListSort(backlist, gl2psTrianglesFirst);
    tree->back = (GL2PSbsptree*)gl2psMalloc(sizeof(GL2PSbsptree));
    gl2psBuildBspTree(tree->back, backlist);
  }
  else{
    gl2psListDelete(backlist);
  }

  gl2psListDelete(primitives);
}

static void gl2psTraverseBspTree(GL2PSbsptree *tree, GL2PSxyz eye, GLfloat epsilon,
                                 GLboolean (*compare)(GLfloat f1, GLfloat f2),
                                 void (*action)(void *data), int inverse)
{
  GLfloat result;

  if(!tree) return;

  result = gl2psComparePointPlane(eye, tree->plane);

  if(GL_TRUE == compare(result, epsilon)){
    gl2psTraverseBspTree(tree->back, eye, epsilon, compare, action, inverse);
    if(inverse){
      gl2psListActionInverse(tree->primitives, action);
    }
    else{
      gl2psListAction(tree->primitives, action);
    }
    gl2psTraverseBspTree(tree->front, eye, epsilon, compare, action, inverse);
  }
  else if(GL_TRUE == compare(-epsilon, result)){
    gl2psTraverseBspTree(tree->front, eye, epsilon, compare, action, inverse);
    if(inverse){
      gl2psListActionInverse(tree->primitives, action);
    }
    else{
      gl2psListAction(tree->primitives, action);
    }
    gl2psTraverseBspTree(tree->back, eye, epsilon, compare, action, inverse);
  }
  else{
    gl2psTraverseBspTree(tree->front, eye, epsilon, compare, action, inverse);
    gl2psTraverseBspTree(tree->back, eye, epsilon, compare, action, inverse);
  }
}


static GLint gl2psSplitPrimitive(GL2PSprimitive *prim, GL2PSplane plane,
                                 GL2PSprimitive **front, GL2PSprimitive **back)
{
  GLshort i, j, in = 0, out = 0, in0[5], in1[5], out0[5], out1[5];
  GLint type;
  GLfloat d[5];

  type = GL2PS_COINCIDENT;

  for(i = 0; i < prim->numverts; i++){
    d[i] = gl2psComparePointPlane(prim->verts[i].xyz, plane);
  }

  switch(prim->type){
  case GL2PS_POINT :
    if(d[0] > GL2PS_EPSILON)       type = GL2PS_IN_BACK_OF;
    else if(d[0] < -GL2PS_EPSILON) type = GL2PS_IN_FRONT_OF;
    else                           type = GL2PS_COINCIDENT;
    break;
  default :
    for(i = 0; i < prim->numverts; i++){
      j = gl2psGetIndex(i, prim->numverts);
      if(d[j] > GL2PS_EPSILON){
        if(type == GL2PS_COINCIDENT)      type = GL2PS_IN_BACK_OF;
        else if(type != GL2PS_IN_BACK_OF) type = GL2PS_SPANNING;
        if(d[i] < -GL2PS_EPSILON){
          gl2psAddIndex(in0, in1, &in, i, j);
          gl2psAddIndex(out0, out1, &out, i, j);
          type = GL2PS_SPANNING;
        }
        gl2psAddIndex(out0, out1, &out, j, -1);
      }
      else if(d[j] < -GL2PS_EPSILON){
        if(type == GL2PS_COINCIDENT)       type = GL2PS_IN_FRONT_OF;
        else if(type != GL2PS_IN_FRONT_OF) type = GL2PS_SPANNING;
        if(d[i] > GL2PS_EPSILON){
          gl2psAddIndex(in0, in1, &in, i, j);
          gl2psAddIndex(out0, out1, &out, i, j);
          type = GL2PS_SPANNING;
        }
        gl2psAddIndex(in0, in1, &in, j, -1);
      }
      else{
        gl2psAddIndex(in0, in1, &in, j, -1);
        gl2psAddIndex(out0, out1, &out, j, -1);
      }
    }
    break;
  }

  if(type == GL2PS_SPANNING){
    *back = (GL2PSprimitive*)gl2psMalloc(sizeof(GL2PSprimitive));
    *front = (GL2PSprimitive*)gl2psMalloc(sizeof(GL2PSprimitive));
    gl2psCreateSplitPrimitive(prim, plane, *back, out, out0, out1);
    gl2psCreateSplitPrimitive(prim, plane, *front, in, in0, in1);
  }

  return type;
}



*/

} // graphics

} // avro

#endif
