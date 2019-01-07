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
  double* u;
}
var_t;

static void initial_condition(p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* q)
{
  //var_t* data = (var_t *) q->p.user_data;
  //data->u[0] = 1.;
}

void p4_new (int level, p4est_t** p4est_out)
{
  sc_MPI_Comm mpicomm;
  p4est_connectivity_t* conn;
  p4est_t* p4est;

  mpicomm = sc_MPI_COMM_WORLD;
  conn = p4est_connectivity_new_periodic();
  p4est = p4est_new_ext(mpicomm,conn,0,level,1,sizeof(var_t),initial_condition,NULL);
  *p4est_out=p4est;
}

void p4_build_mesh(p4est_t* p4est, p4est_topidx_t* tt_out, p4est_mesh_t** mesh_out, sc_array_t** quadrants_out, p4est_nodes_t** nodes_out, int** edges_out, int* np, int* nc, int* ne)
{
  p4est_topidx_t  tt;
  sc_array_t* trees;
  p4est_tree_t* tree;
  p4est_nodes_t* nodes;
  sc_array_t* quadrants;
  p4est_mesh_t* mesh;
  p4est_ghost_t* ghost;
  int i;
  int k;
  int neigh;
  p4est_quadrant_t* quad;
  p4est_quadrant_t* quad_neigh;
  int* edges;
  int* half_neigh;

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

  *ne=1;
  edges = (int*) malloc(sizeof(int)*8*(*nc));

  for (i=0;i<mesh->local_num_quadrants;i++)
  {
    quad=sc_array_index(quadrants,i);
    for (k=0;k<4;k++)
    {
      if (mesh->quad_to_face[4*i+k]>-1)
      {
        quad_neigh=sc_array_index(quadrants,mesh->quad_to_quad[4*i+k]);
        neigh=mesh->quad_to_quad[4*i+k];
        if (i<neigh)
        {
          edges[4*i+k]=*ne;
          switch(k)
          {
            case 0:
            {edges[4*neigh+1]=-*ne;}
            break;
            case 1:
            {edges[4*neigh+0]=-*ne;}
            break;
            case 2:
            {edges[4*neigh+3]=-*ne;}
            break;
            case 3:
            {edges[4*neigh+2]=-*ne;}
            break;
          }
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
      else
      {
        half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i+k]);
        quad_neigh=sc_array_index(quadrants,half_neigh[0]);
        if (i<half_neigh[0])
        {
          edges[4*i+k]=*ne;
          switch(k)
          {
            case 0:
            {
              edges[4*half_neigh[0]+1]=-*ne;
              edges[4*half_neigh[1]+1]=-*ne-1;
            }
            break;
            case 1:
            {
              edges[4*half_neigh[0]]=-*ne;
              edges[4*half_neigh[1]]=-*ne-1;
            }
            break;
            case 2:
            {
              edges[4*half_neigh[0]+3]=-*ne;
              edges[4*half_neigh[1]+3]=-*ne-1;
            }
            break;
            case 3:
            {
              edges[4*half_neigh[0]+2]=-*ne;
              edges[4*half_neigh[1]+2]=-*ne-1;
            }
            break;
          }
          *ne=*ne+2;
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
                *ne=*ne+2;
              }
            }
            break;
            case 1:
            {
              if (quad_neigh->x<quad->x)
              {
                edges[4*i+k]=*ne;
                *ne=*ne+2;
              }
            }
            break;
            case 2:
            {
              if (quad_neigh->y>quad->y)
              {
                edges[4*i+k]=*ne;
                *ne=*ne+2;
              }
            }
            break;
            case 3:
            {
              if (quad_neigh->y<quad->y)
              {
                edges[4*i+k]=*ne;
                *ne=*ne+2;
              }
            }
            break;
          }
        }
      }
    }
  }
  *edges_out=edges;
  *ne=*ne-1;
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

void p4_get_cell(p4est_t* p4est, p4est_mesh_t* mesh, p4est_locidx_t tt, sc_array_t* quadrants, p4est_nodes_t* nodes, int* edges, int i, double* Xc_out, double* Yc_out, double* dX_out, double* dY_out, int** corners_out, int* Nneigh, int** neighbors_out, int* Nnodes, int** nodes_out, int* N_edge, int* lev)
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
  int* half_neigh;
  int* cell_nodes;
  int inode;

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
  cell_nodes = (int*) malloc(sizeof(int)*8);
  neighbors = (int*) malloc(sizeof(int)*12);
  *Nneigh = 0;
  inode = 0;

  corners[0]=(int) nodes->local_nodes[4*i]+1;
  corners[1]=(int) nodes->local_nodes[4*i+1]+1;
  corners[2]=(int) nodes->local_nodes[4*i+2]+1;
  corners[3]=(int) nodes->local_nodes[4*i+3]+1;
  if (mesh->quad_to_face[4*i+2]<0)
  {
    half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i+2]);
    quad_neigh=sc_array_index(quadrants,half_neigh[0]);
    if (quad_neigh->y<quad->y)
    {
      cell_nodes[inode]=nodes->local_nodes[4*half_neigh[0]+2]+1;
      cell_nodes[inode+1]=nodes->local_nodes[4*half_neigh[1]+2]+1;
      inode=inode+2;
    }
    else
    {
      cell_nodes[inode]=nodes->local_nodes[4*i]+1;
      inode=inode+1;
    }
  }
  else
  {
    cell_nodes[inode]=nodes->local_nodes[4*i]+1;
    inode=inode+1;
  }
  if (mesh->quad_to_face[4*i+1]<0)
  {
    half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i+1]);
    quad_neigh=sc_array_index(quadrants,half_neigh[0]);
    if (quad_neigh->x>quad->x)
    {
      cell_nodes[inode]=nodes->local_nodes[4*half_neigh[0]]+1;
      cell_nodes[inode+1]=nodes->local_nodes[4*half_neigh[1]]+1;
      inode=inode+2;
    }
    else
    {
      cell_nodes[inode]=nodes->local_nodes[4*i+1]+1;
      inode=inode+1;
    }
  }
  else
  {
    cell_nodes[inode]=nodes->local_nodes[4*i+1]+1;
    inode=inode+1;
  }
  if (mesh->quad_to_face[4*i+3]<0)
  {
    half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i+3]);
    quad_neigh=sc_array_index(quadrants,half_neigh[0]);
    if (quad_neigh->y>quad->y)
    {
      cell_nodes[inode]=nodes->local_nodes[4*half_neigh[1]+1]+1;
      cell_nodes[inode+1]=nodes->local_nodes[4*half_neigh[0]+1]+1;
      inode=inode+2;
    }
    else
    {
      cell_nodes[inode]=nodes->local_nodes[4*i+3]+1;
      inode=inode+1;
    }
  }
  else
  {
    cell_nodes[inode]=nodes->local_nodes[4*i+3]+1;
    inode=inode+1;
  }
  if (mesh->quad_to_face[4*i]<0)
  {
    half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i]);
    quad_neigh=sc_array_index(quadrants,half_neigh[0]);
    if (quad_neigh->x<quad->x)
    {
      cell_nodes[inode]=nodes->local_nodes[4*half_neigh[1]+3]+1;
      cell_nodes[inode+1]=nodes->local_nodes[4*half_neigh[0]+3]+1;
      inode=inode+2;
    }
    else
    {
      cell_nodes[inode]=nodes->local_nodes[4*i+2]+1;
      inode=inode+1;
    }
  }
  else
  {
    cell_nodes[inode]=nodes->local_nodes[4*i+2]+1;
    inode=inode+1;
  }
  *Nnodes=inode;

  for (k=0;k<4;k++)
  {
    if (mesh->quad_to_face[4*i+k]<0)
    {
      half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i+k]);
      neighbors[*Nneigh]=half_neigh[0]+1;
      neighbors[*Nneigh+1]=half_neigh[1]+1;
      *Nneigh=*Nneigh+2;
    }
    else
    {
      quad_neigh=sc_array_index(quadrants,mesh->quad_to_quad[4*i+k]);
      neighbors[*Nneigh]=mesh->quad_to_quad[4*i+k]+1;
      switch(k)
      {
        case 0:
        {
          if (quad_neigh->x>quad->x){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
        case 1:
        {
          if (quad_neigh->x<quad->x){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
        case 2:
        {
          if (quad_neigh->y>quad->y){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
        case 3:
        {
          if (quad_neigh->y<quad->y){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
      }
      *Nneigh=*Nneigh+1;
    }
  }
  *N_edge=*Nneigh;
  for (k=4;k<8;k++)
  {
    index = mesh->quad_to_corner[4*i+k-4];
    if (index>-1)
    {
      if (index+1<=mesh->local_num_quadrants)
      {
        neighbors[*Nneigh] = index+1;
      }
      else
      {
        index = index-mesh->local_num_quadrants-mesh->ghost_num_quadrants;
        neigh = sc_array_index(mesh->corner_quad,index);
        neighbors[*Nneigh]=*neigh+1;
      }
      quad_neigh = sc_array_index(quadrants,neighbors[*Nneigh]-1);
      switch(k-4)
      {
        case 0:
        {
          if (quad_neigh->x>quad->x||quad_neigh->y>quad->y){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
        case 1:
        {
          if (quad_neigh->x<quad->x||quad_neigh->y>quad->y){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
        case 2:
        {
          if (quad_neigh->x>quad->x||quad_neigh->y<quad->y){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
        case 3:
        {
          if (quad_neigh->x<quad->x||quad_neigh->y<quad->y){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
      }
      *Nneigh=*Nneigh+1;
    }
  }
  *neighbors_out=neighbors;
  *corners_out=corners;
  *nodes_out=cell_nodes;
}

void p4_get_edge(p4est_t *p4est, p4est_mesh_t* mesh, sc_array_t* quadrants, int* edges, int k, int i, int** iedge_out, int* nedge, int** cell1_out, int** cell2_out, int** sub_out, int** period_out)
{
  p4est_quadrant_t *quad;
  p4est_quadrant_t *quad_neigh;
  int neigh;
  int half;
  int* sub;
  int* cell1;
  int* cell2;
  int* half_neigh;
  int* iedge;
  int* period;

  half=mesh->quad_to_face[4*k+i];
  if (half>-1)
  {
    iedge = (int*) malloc(sizeof(int));
    iedge[0]=edges[4*k+i];
    *nedge=1;
    neigh=mesh->quad_to_quad[4*k+i];

    quad=sc_array_index(quadrants,k);
    quad_neigh=sc_array_index(quadrants,mesh->quad_to_quad[4*k+i]);
    cell1 = (int*) malloc(sizeof(int));
    cell2 = (int*) malloc(sizeof(int));
    sub = (int*) malloc(sizeof(int)*2);
    period = (int*) malloc(sizeof(int));

    switch(i)
    {
      case 0:
      {
        cell1[0]=neigh+1;
        cell2[0]=k+1;
        if (quad_neigh->x>quad->x)
        {
          cell1[0]=-cell1[0];
          period[0]=edges[4*neigh+1];
        }
        else {period[0]=edges[4*k];}
        if (half>7)
        {
          sub[0]=1;
          sub[1]=0;
          half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*neigh+1]);
          if (k==half_neigh[1]) {period[0]=period[0]+1;}
        }
        else {sub[0]=0;sub[1]=0;}
      }
      break;
      case 1:
      {
        cell1[0]=k+1;
        cell2[0]=neigh+1;
        if (quad_neigh->x<quad->x)
        {
          cell2[0]=-cell2[0];
          period[0]=edges[4*neigh];
        }
        else {period[0]=edges[4*k+1];}
        if (half>7)
        {
          sub[0]=0;
          sub[1]=1;
          half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*neigh]);
          if (k==half_neigh[1]) {period[0]=period[0]+1;}
        }
        else {sub[0]=0;sub[1]=0;}
      }
      break;
      case 2:
      {
        cell1[0]=neigh+1;
        cell2[0]=k+1;
        if (quad_neigh->y>quad->y)
        {
          cell1[0]=-cell1[0];
          period[0]=edges[4*neigh+3];
        }
        else {period[0]=edges[4*k+2];}
        if (half>7)
        {
          sub[0]=1;
          sub[1]=0;
          half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*neigh+3]);
          if (k==half_neigh[1]) {period[0]=period[0]+1;}
        }
        else {sub[0]=0;sub[1]=0;}
      }
      break;
      case 3:
      {
        cell1[0]=k+1;
        cell2[0]=neigh+1;
        if (quad_neigh->y<quad->y)
        {
          cell2[0]=-cell2[0];
          period[0]=edges[4*neigh+2];
        }
        else {period[0]=edges[4*k+3];}
        if (half>7)
        {
          sub[0]=0;
          sub[1]=1;
          half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*neigh+2]);
          if (k==half_neigh[1]) {period[0]=period[0]+1;}
        }
        else {sub[0]=0;sub[1]=0;}
      }
      break;
    }
  }
  else
  {
    half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*k+i]);
    iedge = (int*) malloc(sizeof(int)*2);
    iedge[0]=edges[4*k+i];
    iedge[1]=edges[4*k+i]+1;
    quad=sc_array_index(quadrants,k);
    quad_neigh=sc_array_index(quadrants,half_neigh[0]);
    sub = (int*) malloc(sizeof(int)*2);
    cell1 = (int*) malloc(sizeof(int)*2);
    cell2 = (int*) malloc(sizeof(int)*2);
    period = (int*) malloc(sizeof(int)*2);

    switch(i)
    {
      case 0:
      {
        *nedge=2;
        cell1[0]=half_neigh[0]+1;
        cell1[1]=half_neigh[1]+1;
        cell2[0]=k+1;
        cell2[1]=k+1;
        if (quad_neigh->x>quad->x)
        {
          cell1[0]=-cell1[0];
          cell1[1]=-cell1[1];
          period[0]=edges[4*half_neigh[0]+1];
          period[1]=edges[4*half_neigh[1]+1];
        }
        else
        {
          period[0]=edges[4*k];
          period[1]=edges[4*k]+1;
        }
        sub[0]=0;sub[1]=1;
      }
      break;
      case 1:
      {
        *nedge=2;
        cell1[0]=k+1;
        cell1[1]=k+1;
        cell2[0]=half_neigh[0]+1;
        cell2[1]=half_neigh[1]+1;
        if (quad_neigh->x<quad->x)
        {
          cell2[0]=-cell2[0];
          cell2[1]=-cell2[1];
          period[0]=edges[4*half_neigh[0]];
          period[1]=edges[4*half_neigh[1]];
        }
        else
        {
          period[0]=edges[4*k+1];
          period[1]=edges[4*k+1]+1;
        }
        sub[0]=1;sub[1]=0;
      }
      break;
      case 2:
      {
        *nedge=2;
        cell1[0]=half_neigh[0]+1;
        cell1[1]=half_neigh[1]+1;
        cell2[0]=k+1;
        cell2[1]=k+1;
        if (quad_neigh->y>quad->y)
        {
          cell1[0]=-cell1[0];
          cell1[1]=-cell1[1];
          period[0]=edges[4*half_neigh[0]+3];
          period[1]=edges[4*half_neigh[1]+3];
        }
        else
        {
          period[0]=edges[4*k+2];
          period[1]=edges[4*k+2]+1;
        }
        sub[0]=0;sub[1]=1;
      }
      break;
      case 3:
      {
        *nedge=2;
        cell1[0]=k+1;
        cell1[1]=k+1;
        cell2[0]=half_neigh[0]+1;
        cell2[1]=half_neigh[1]+1;
        if (quad_neigh->y<quad->y)
        {
          cell2[0]=-cell2[0];
          cell2[1]=-cell2[1];
          period[0]=edges[4*half_neigh[0]+2];
          period[1]=edges[4*half_neigh[1]+2];
        }
        else
        {
          period[0]=edges[4*k+3];
          period[1]=edges[4*k+3]+1;
        }
        sub[0]=1;sub[1]=0;
      }
      break;
    }
  }
  *sub_out=sub;
  *cell1_out=cell1;
  *cell2_out=cell2;
  *iedge_out=iedge;
  *period_out=period;
}

static int refine_test (p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* q)
{
  double X[3];
  p4est_topidx_t tt;
  var_t* data;

  tt = p4est->first_local_tree;
  data = (var_t *) q->p.user_data;
  if (data->u[0]>0.5) {return 1;}
  //p4est_qcoord_to_vertex (p4est->connectivity, tt, q->x, q->y, X);
  //if (X[0]<0.2&&X[1]<0.2){return 1;}
  else {return 0;}
}

void p4_refine (p4est_t* p4est, sc_array_t* quadrants, double* sol, int nsol, int refine_recursive, char* refine_fn, char* init_fn)
{
  int k;
  int isol;
  p4est_quadrant_t* q;
  var_t* data;

  printf("Refine function : |%s| \ninit_fn : |%s| \n",refine_fn,init_fn);
  for (k=0;k<quadrants->elem_count;k++)
  {
    q=sc_array_index(quadrants,k);
    data = (var_t *) q->p.user_data;
    data->u = (double*) malloc(sizeof(double)*nsol);
    for (isol=0;isol<nsol;isol++)
    {
      data->u[isol]=sol[isol*quadrants->elem_count+k];
    }
  }
  p4est_refine(p4est,refine_recursive,refine_test,NULL);
}
