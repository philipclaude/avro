#include "adaptation/geometry.h"

namespace avro
{

// retrieves the geometric parameter coordinates given a facet on the entity
void
geometry_params( Entity* e0 , const Points& points , const index_t* v , const index_t nv , real_t* params )
{
  avro_assert( nv>0 );

  coord_t udim = points.udim();

  EGADS::Object* e = (EGADS::Object*) e0;

  // get all the entities these points are on
  index_t count = 0;
  std::vector<EGADS::Object*> entities( nv );
  std::vector<bool> ufound( nv, false );
  for (index_t k=0;k<nv;k++)
  {
    entities[k] = (EGADS::Object*) points.entity( v[k] );

    // if the vertex entity is that of the requested one, set the uv coords
    if (entities[k]==e)
    {
      for (int i=0;i<udim;i++)
         params[k*udim+i] = points.u( v[k] )[i];
      ufound[k] = true;
      count++;
    }

    if (entities[k]==nil) points.print(true);

    avro_assert(entities[k] != nil); // all points must have entities
  }

  if (count==nv) return; // all parameters have been retrieved

  //---------------------------------------//
  if (e->object_class()==FACE)
  {
    avro_assert(udim == 2);

    std::vector<std::array<double,2>> uvp(nv), uvm(nv);

    EGADS::Object* face = e;
    for (index_t k=0;k<nv;k++)
    {
      // skip any parametric coordinates already set
      if (ufound[k]) continue;

      // try to find an edge and t-value to get uv
      double t = 0;
      EGADS::Object* edge = nil;

      if (entities[k]->object_class()==NODE)
      {
        // find an edge that is a parent of the node and a child of the face
        EGADS::Object* node = entities[k];
        for (index_t i=0;i<node->nb_parents();i++)
        {
          edge = (EGADS::Object*) node->parents(i);

          if (edge->object_class()!=EDGE) continue;

          // skip edges that define periodicity in parameter space if possible
          if (edge->sense_required() && i<node->nb_parents()-1 ) continue;

          // if the face is a parent of the this edge
          if (edge->has_parent(face))
          {
            int periodic;
            double trange[2];
            EGADS_ENSURE_SUCCESS( EG_getRange(*edge->object(), trange, &periodic) );

            #if 0
            if (periodic==1)
            {
              t = trange[0]; // TODO: Is this right?!?!?

              // @marshall: changed to the same logic as non-periodic
              if (edge->egchild(0) == *node->object()) t = trange[0];
              if (edge->egchild(1) == *node->object()) t = trange[1];
            }
            else
            {
              if (edge->egchild(0) == *node->object()) t = trange[0];
              if (edge->egchild(1) == *node->object()) t = trange[1];
            }
            #else
            if (edge->egchild(0) == *node->object()) t = trange[0];
            if (edge->egchild(1) == *node->object()) t = trange[1];
            #endif

            break;
          }
          else
            edge = nil;
        }
      }
      else if (entities[k]->object_class()==EDGE)
      {
        edge = entities[k];
        t = points.u( v[k] )[0];
      }
      else
      {
        printf("could not find parameter coordinates along entity for vertex (%lu) %lu\n",k,v[k]);
        points.print(v[k],true);
        printf("%s\n",EGADS::utilities::object_class_name(entities[k]->object_class()).c_str());
        entities[k]->print_header();
        avro_assert_not_reached;
      }

      avro_assert( edge!=nil );
      if (edge->sense_required())
      {
        // need to sort out which uv should be used later
        EGADS_ENSURE_SUCCESS( EG_getEdgeUV( *face->object() , *edge->object() , -1 , t , uvm[k].data() ) );
        EGADS_ENSURE_SUCCESS( EG_getEdgeUV( *face->object() , *edge->object() ,  1 , t , uvp[k].data() ) );
      }
      else
      {
        // get the uv value from the edge t value. No need to worry about periodicity in uv
        double uv[2] = {0,0};
        int status = EG_getEdgeUV( *face->object() , *edge->object() , 0 , t , uv );

        // special treatment of interior edges that might not be connected to the face
        if (status == EGADS_NOTFOUND && edge->interior())
        {
          bool foundGuess = false;
          for (index_t i=0;i<nv;i++)
          {
            if (!ufound[i]) continue;
            uv[0] = params[i*udim  ];
            uv[1] = params[i*udim+1];
            avro_implement;
            foundGuess = true;
            break;
          }

          double x_inv[3];
          double* x = const_cast<double*>(points[ v[k] ]);
          if (foundGuess)
          {
            EGADS_ENSURE_SUCCESS( EG_invEvaluateGuess( *face->object() , x , uv , x_inv ) );
          }
          else
          {
            EGADS_ENSURE_SUCCESS( EG_invEvaluate( *face->object() , x , uv , x_inv ) );
          }

          real_t d = numerics::distance2(x,x_inv,3);
          real_t tol_edge=0, tol_face=0;
          EGADS_ENSURE_SUCCESS( EG_tolerance( *edge->object(), &tol_edge ) );
          EGADS_ENSURE_SUCCESS( EG_getTolerance( *face->object(), &tol_face ) );
          avro_assert_msg( d <= std::min(tol_edge,tol_face) ,
                           "d = %1.16e, tol_edge = %1.16e, tol_face = %1.16e" , d , tol_edge , tol_face );
        }
        else
          EGADS_ENSURE_SUCCESS( status );

        params[k*udim  ] = uv[0];
        params[k*udim+1] = uv[1];

        ufound[k] = true;
        count++;
      }
    }
    avro_assert(count > 0); // at least one uv value must be set...

    for (index_t k=0;k<nv;k++)
    {
      if (ufound[k]) continue;

      // get a uv-value that has already been set
      double uv[2] = {0,0};
      for (index_t j=0;j<nv;j++)
      {
        if (j==k || !ufound[k]) continue;
        uv[0] = params[j*udim  ];
        uv[1] = params[j*udim+1];
        break;
      }

      double mindist = std::numeric_limits<double>::max();
      index_t jmin = nv;
      int sense = 0;
      for (index_t j=0;j<nv;j++)
      {
        if (j == k) continue;

        double distm = sqrt(pow(uvm[j][0] - uv[0], 2) + pow(uvm[j][1] - uv[1], 2));
        double distp = sqrt(pow(uvp[j][0] - uv[0], 2) + pow(uvp[j][1] - uv[1], 2));
        if (distm < mindist)
        {
          mindist = distm;
          jmin = j;
          sense = -1;
        }
        if (distp < mindist)
        {
          mindist = distp;
          jmin = j;
          sense = 1;
        }
      }

      if (sense==1)
      {
        params[k*udim  ] = uvp[jmin][0];
        params[k*udim+1] = uvp[jmin][1];
      }
      else
      {
        params[k*udim  ] = uvm[jmin][0];
        params[k*udim+1] = uvm[jmin][1];
      }
      ufound[k] = true;
      count++;

      if (count == nv) break; // all parameters found
    }
  }

  //---------------------------------------//
  else if (e->object_class()==EDGE)
  {
    avro_assert(udim >= 1);

    EGADS::Object* edge = e;

    for (index_t k=0;k<nv;k++)
    {
      // skip any parametric coordinates already set
      if (ufound[k]) continue;

      // try to find an edge and t-value to get uv
      double t = 0;

      if (entities[k]->object_class()==NODE)
      {
        // find an edge that is a parent of the node and a child of the face
        EGADS::Object* node = entities[k];

        int periodic;
        double trange[2];
        EGADS_ENSURE_SUCCESS( EG_getRange(*edge->object(), trange, &periodic) );
        #if 0
        if (periodic==1)
        {
          t = trange[0]; // TODO: Is this right?!?!?

          // @marshall: changed to the same logic as non-periodic
          if (edge->egchild(0) == *node->object()) t = trange[0];
          if (edge->egchild(1) == *node->object()) t = trange[1];
        }
        else
        {
          if (edge->egchild(0) == *node->object()) t = trange[0];
          if (edge->egchild(1) == *node->object()) t = trange[1];
        }
        #else
        if (edge->egchild(0) == *node->object()) t = trange[0];
        if (edge->egchild(1) == *node->object()) t = trange[1];
        #endif

        // set the parametric value
        params[k*udim] = t;

        ufound[k] = true;
        count++;
      }
      else
        avro_assert(false); // this should not happen...

      if (count == nv) break; // all parameters found
    }
  }

  avro_assert(count == nv);

#if 0
  coord_t dim  = points.dim();
  for (index_t k=0;k<nv;k++)
  {
    // evaluate the coordinates for the parameters we found
    double x_eval[18];
    EGADS_ENSURE_SUCCESS( EG_evaluate( *e->object() , &params[k*udim] , x_eval ) );

    // check the evaluated coordinates are close to the true ones
    const real_t* x = points[ v[k] ];
    real_t d = numerics::distance2(x,x_eval,dim);
    real_t tol = 0;
    EGADS_ENSURE_SUCCESS( EG_getTolerance( *e->object(), &tol ) );

    if (d>tol)
    {
      printf("Failure on facet vertex: %lu\n",v[k]);
      if (nv == 2)
        printf("facet = (%lu,%lu)\n",v[0],v[1]);
      else if (nv == 3)
        printf("facet = (%lu,%lu,%lu)\n",v[0],v[1],v[2]);
      printf("x_true = (%g,%g,%g), x_eval = (%g,%g,%g) with param coordinates (%g,%g)\n",
                x[0],x[1],x[2],x_eval[0],x_eval[1],x_eval[2],params[k*udim],params[k*udim+1]);
      printf("Getting parameters for entity:\n");
      e->print(false);
      printf("Entity storing vertex with param = (%g,%g):\n", points.u( v[k] )[0], points.u( v[k] )[1]);
      entities[k]->print(false);
    }

    avro_assert_msg( d <= tol , "d = %1.16e, tol = %1.16e" , d , tol );
  }
#endif
}

}
