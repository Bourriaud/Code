#include <stdlib.h>
#include <stdio.h>
#include <p4est.h>
#include <p4est_connectivity.h>
#include <p4est_geometry.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#include <p4est_nodes.h>
#include <p4est_bits.h>

static sc_MPI_Comm mpicomm;

typedef struct var
{
  double u;
  double du[P4EST_DIM];
}
var_t;

static void initial_condition(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
  var_t *data = (var_t *) q->p.user_data;
  data->u = 1.;
}

static void get_value(p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_t *p4est = info->p4est;
  p4est_quadrant_t *q = info->quad;
  p4est_topidx_t      which_tree = info->treeid;
  var_t *data = (var_t *) q->p.user_data;
  sc_array_t *u_interp = (sc_array_t *) user_data;
  double this_u;
  double *this_u_ptr;
  p4est_locidx_t local_id = info->quadid;
  p4est_tree_t *tree;
  p4est_locidx_t arrayoffset;

  tree = p4est_tree_array_index (p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
  arrayoffset = P4EST_CHILDREN * local_id;
  if(local_id>100){this_u=1.;}
  else{this_u=0.;}
  for (int i=0;i<P4EST_CHILDREN;i++)
  {
    this_u_ptr = (double *) sc_array_index (u_interp,arrayoffset+i);
    this_u_ptr[0] = this_u;
  }
}

void p4_test()
{
  sc_MPI_Comm mpicomm;
  p4est_connectivity_t *conn;
  p4est_geometry_t *geom;
  p4est_t *p4est;
  p4est_vtk_context_t *context;
  p4est_locidx_t      numquads;
  const char* filename="/home/imb/abourriaud/Documents/Code/results/salam";
  const char* scalar_name="variable1";
  sc_array_t* values;

  mpicomm = sc_MPI_COMM_WORLD;
  conn = p4est_connectivity_new_periodic();
  geom = p4est_geometry_new_connectivity(conn);
  p4est = p4est_new_ext(mpicomm,conn,0,4,1,sizeof(var_t),initial_condition,NULL);

  numquads = p4est->local_num_quadrants;
  values = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN);
  p4est_iterate (p4est, NULL, (void *) values, get_value, NULL, NULL);
  context = p4est_vtk_context_new (p4est, filename);
  p4est_vtk_context_set_scale (context, 0.99);
  context = p4est_vtk_write_header (context);
  context = p4est_vtk_write_point_dataf (context,1,0,"solution",values,context);
  p4est_vtk_write_footer (context);

  p4est_destroy(p4est);
  p4est_connectivity_destroy(conn);
}

void p4_build_mesh(int level, p4est_t** p4est_out, p4est_topidx_t* tt_out, p4est_mesh_t** mesh_out, sc_array_t** quadrants_out, p4est_nodes_t** nodes_out, int** edges_out, int* np, int* nc, int* ne)
{
  sc_MPI_Comm mpicomm;
  p4est_connectivity_t *conn;
  p4est_t *p4est;
  p4est_topidx_t  tt;
  sc_array_t *trees;
  p4est_tree_t *tree;
  p4est_nodes_t *nodes;
  sc_array_t *quadrants;
  p4est_mesh_t* mesh;
  p4est_ghost_t* ghost;
  int i;
  int k;
  int neigh;
  p4est_quadrant_t *quad;
  p4est_quadrant_t *quad_neigh;
  int* edges;

  mpicomm = sc_MPI_COMM_WORLD;
  conn = p4est_connectivity_new_periodic();
  p4est = p4est_new_ext(mpicomm,conn,0,level,1,sizeof(var_t),initial_condition,NULL);
  *p4est_out=p4est;
  ghost = p4est_ghost_new (p4est,P4EST_CONNECT_FULL);
  mesh = p4est_mesh_new_ext(p4est,ghost,1,1,P4EST_CONNECT_FULL);
  *mesh_out=mesh;

  tt = p4est->first_local_tree;
  *tt_out = tt;
  trees = p4est->trees;
  tree = p4est_tree_array_index (trees, tt);
  nodes = p4est_nodes_new (p4est,NULL);
  *nodes_out = nodes;
  quadrants = &tree->quadrants;
  *quadrants_out = quadrants;
  *np = nodes->indep_nodes.elem_count;
  *nc = quadrants->elem_count;

  *ne=0;
  edges = (int*) malloc(sizeof(int)*4*(*nc));
  for (i=0;i<mesh->local_num_quadrants;i++)
  {
    quad=sc_array_index(quadrants,i);
    for (k=0;k<4;k++)
    {
      quad_neigh=sc_array_index(quadrants,mesh->quad_to_quad[4*i+k]);
      neigh=mesh->quad_to_quad[4*i+k];
      edges[4*i+k]=-1;

      if (i<neigh)
      {
        edges[4*i+k]=*ne;
        *ne=*ne+1;
      }
      else
      {
        switch(k)
        {
          case 0:
          {
            if (quad_neigh->x>quad->x)
            {
              edges[4*i+k]=*ne;
              *ne=*ne+1;
            }
          }
          break;
          case 1:
          {
            if (quad_neigh->x<quad->x)
            {
              edges[4*i+k]=*ne;
              *ne=*ne+1;
            }
          }
          break;
          case 2:
          {
            if (quad_neigh->y>quad->y)
            {
              edges[4*i+k]=*ne;
              *ne=*ne+1;
            }
          }
          break;
          case 3:
          {
            if (quad_neigh->y<quad->y)
            {
              edges[4*i+k]=*ne;
              *ne=*ne+1;
            }
          }
          break;
        }
      }
    }
  }
  *edges_out=edges;
}

void p4_get_node(p4est_t* p4est, p4est_topidx_t tt, p4est_nodes_t* nodes, int i, double* X_out, double* Y_out)
{
  double vxyz[3];
  p4est_indep_t *indep;

  double* X = (double*) malloc(sizeof(double)*nodes->indep_nodes.elem_count);
  double* Y = (double*) malloc(sizeof(double)*nodes->indep_nodes.elem_count);

  indep = (p4est_indep_t *) sc_array_index (&nodes->indep_nodes, (size_t) i);
  p4est_qcoord_to_vertex (p4est->connectivity, tt, indep->x, indep->y, vxyz);

  *X_out=vxyz[0];
  *Y_out=vxyz[1];
}

void p4_get_cell(p4est_t* p4est, p4est_mesh_t* mesh, p4est_locidx_t tt, sc_array_t* quadrants, p4est_nodes_t* nodes, int* edges, int i, double* Xc_out, double* Yc_out, double* dX_out, double* dY_out, int** corners_out, int** neighbors_out, int* lev)
{
  p4est_quadrant_t *quad;
  size_t num_quads;
  p4est_qcoord_t len;
  double Xc[3],dX[3];
  int* corners;
  int* neighbors;
  int* neigh;
  int k;
  int index;
  p4est_quadrant_t* quad_neigh;

  num_quads = quadrants->elem_count;

  quad = p4est_quadrant_array_index (quadrants,(size_t) i);
  *lev = quad->level;
  len = P4EST_QUADRANT_LEN (*lev);
  p4est_qcoord_to_vertex (p4est->connectivity, tt, len, len, dX);
  p4est_qcoord_to_vertex (p4est->connectivity, tt, quad->x, quad->y, Xc);
  Xc[0]=Xc[0]+dX[0]/2.;
  Xc[1]=Xc[1]+dX[1]/2.;

  *Xc_out=Xc[0];
  *Yc_out=Xc[1];
  *dX_out=dX[0];
  *dY_out=dX[1];

  corners = (int*) malloc(sizeof(int)*4);
  corners[0]=(int) nodes->local_nodes[4*i]+1;
  corners[1]=(int) nodes->local_nodes[4*i+1]+1;
  corners[2]=(int) nodes->local_nodes[4*i+2]+1;
  corners[3]=(int) nodes->local_nodes[4*i+3]+1;
  *corners_out=corners;

  neighbors = (int*) malloc(sizeof(int)*8);
  for (k=0;k<4;k++)
  {
    quad_neigh=sc_array_index(quadrants,mesh->quad_to_quad[4*i+k]);
    neighbors[k]=mesh->quad_to_quad[4*i+k]+1;
    switch(k)
    {
      case 0:
      {
        if (quad_neigh->x>quad->x){neighbors[k]=-neighbors[k];}
      }
      break;
      case 1:
      {
        if (quad_neigh->x<quad->x){neighbors[k]=-neighbors[k];}
      }
      break;
      case 2:
      {
        if (quad_neigh->y>quad->y){neighbors[k]=-neighbors[k];}
      }
      break;
      case 3:
      {
        if (quad_neigh->y<quad->y){neighbors[k]=-neighbors[k];}
      }
      break;
    }
  }
  for (k=4;k<8;k++)
  {
    index = mesh->quad_to_corner[4*i+k-4];
    if (index+1<=mesh->local_num_quadrants){neighbors[k] = index+1;}
    else
    {
      index = index-mesh->local_num_quadrants-mesh->ghost_num_quadrants;
      neigh = sc_array_index(mesh->corner_quad,index);
      neighbors[k]=*neigh+1;
    }
    quad_neigh = sc_array_index(quadrants,neighbors[k]-1);
    switch(k-4)
    {
      case 0:
      {
        if (quad_neigh->x>quad->x||quad_neigh->y>quad->y){neighbors[k]=-neighbors[k];}
      }
      break;
      case 1:
      {
        if (quad_neigh->x<quad->x||quad_neigh->y>quad->y){neighbors[k]=-neighbors[k];}
      }
      break;
      case 2:
      {
        if (quad_neigh->x>quad->x||quad_neigh->y<quad->y){neighbors[k]=-neighbors[k];}
      }
      break;
      case 3:
      {
        if (quad_neigh->x<quad->x||quad_neigh->y<quad->y){neighbors[k]=-neighbors[k];}
      }
      break;
    }
  }
  *neighbors_out=neighbors;
}

void p4_get_edge(p4est_t *p4est, p4est_mesh_t* mesh, sc_array_t* quadrants, int* edges, int ne, int k, int i, int* iedge, int* cell1, int* cell2, int** sub_out)
{
  p4est_quadrant_t *quad;
  p4est_quadrant_t *quad_neigh;
  int neigh;
  int half;
  int* sub;

  *iedge=edges[4*k+i]+1;
  neigh=mesh->quad_to_quad[4*k+i];
  *cell1=k+1;
  *cell2=neigh+1;

  quad=sc_array_index(quadrants,k);
  quad_neigh=sc_array_index(quadrants,mesh->quad_to_quad[4*k+i]);
  half=mesh->quad_to_face[4*k+i];
  sub = (int*) malloc(sizeof(int)*2);

  switch(i)
  {
    case 0:
    {
      *cell1=neigh+1;
      *cell2=k+1;
      if (quad_neigh->x>quad->x)
      {
        *cell1=-*cell1;
      }
      if (half>7) {sub[0]=0;sub[1]=1;}
      else if (half<0) {sub[0]=1;sub[1]=0;}
      else {sub[0]=0;sub[1]=0;}
    }
    break;
    case 1:
    {
      *cell1=k+1;
      *cell2=neigh+1;
      if (quad_neigh->x<quad->x)
      {
        *cell2=-*cell2;
      }
      if (half>7) {sub[0]=1;sub[1]=0;}
      else if (half<0) {sub[0]=0;sub[1]=1;}
      else {sub[0]=0;sub[1]=0;}
    }
    break;
    case 2:
    {
      *cell1=neigh+1;
      *cell2=k+1;
      if (quad_neigh->y>quad->y)
      {
        *cell1=-*cell1;
      }
      if (half>7) {sub[0]=0;sub[1]=1;}
      else if (half<0) {sub[0]=1;sub[1]=0;}
      else {sub[0]=0;sub[1]=0;}
    }
    break;
    case 3:
    {
      *cell1=k+1;
      *cell2=neigh+1;
      if (quad_neigh->y<quad->y)
      {
        *cell2=-*cell2;
      }
      if (half>7) {sub[0]=1;sub[1]=0;}
      else if (half<0) {sub[0]=0;sub[1]=1;}
      else {sub[0]=0;sub[1]=0;}
    }
    break;
  }
  *sub_out=sub;
}
