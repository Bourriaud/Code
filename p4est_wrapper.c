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
  double *u,*u1,*u2,*u3,*u4;
  double *pol;
  int coarsen;
  int refine;
  int nsol;
  int size_pol;
  int k;
}
var_t;

static void initial_condition(p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* q)
{
  //var_t* data = (var_t *) q->p.user_data;
  //data->u[0] = 1.;
}

void p4_new (int level, p4est_connectivity_t** conn_out, p4est_t** p4est_out)
{
  sc_MPI_Comm mpicomm;
  p4est_connectivity_t* conn;
  p4est_t* p4est;

  mpicomm = sc_MPI_COMM_WORLD;
  conn = p4est_connectivity_new_periodic();
  p4est = p4est_new_ext(mpicomm,conn,0,level,1,sizeof(var_t),initial_condition,NULL);
  *conn_out=conn;
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
  var_t* data;

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
    data=quad->p.user_data;
    data->k=i;
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
  p4est_ghost_destroy(ghost);
}

void p4_get_node(p4est_t* p4est, p4est_topidx_t tt, p4est_nodes_t* nodes, int i, double* X_out, double* Y_out)
{
  double vxyz[3];
  p4est_indep_t *indep;

  indep = sc_array_index (&nodes->indep_nodes, (size_t) i);
  p4est_qcoord_to_vertex (p4est->connectivity, tt, indep->x, indep->y, vxyz);

  *X_out=vxyz[0];
  *Y_out=vxyz[1];
}

void p4_get_cell(p4est_t* p4est, p4est_mesh_t* mesh, p4est_locidx_t tt, sc_array_t* quadrants, p4est_nodes_t* nodes, int* edges, int i, double* Xc_out, double* Yc_out, double* dX_out, double* dY_out, int** corners_out, int** corners_cell_out, int* Nneigh, int** neighbors_out, int* Nnodes, int** nodes_out, int* N_edge, int* lev)
{
  p4est_quadrant_t *quad;
  size_t num_quads;
  p4est_qcoord_t len;
  double Xc[3],dX[3];
  int* corners;
  int* corners_cell;
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
  corners_cell = (int*) malloc(sizeof(int)*4);
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
          if (quad_neigh->y+dX[1]*2.>quad->y+dX[1])
          {
            corners_cell[2] = mesh->quad_to_quad[4*i]+1;
            if (quad_neigh->x>quad->x){corners_cell[2]=-corners_cell[2];}
          }
          if (quad_neigh->y+dX[1]*2.<quad->y+dX[1])
          {
            corners_cell[0] = mesh->quad_to_quad[4*i]+1;
            if (quad_neigh->x>quad->x){corners_cell[0]=-corners_cell[0];}
          }
          if (quad_neigh->x>quad->x){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
        case 1:
        {
          if (quad_neigh->y+dX[1]*2.>quad->y+dX[1])
          {
            corners_cell[3] = mesh->quad_to_quad[4*i+1]+1;
            if (quad_neigh->x<quad->x){corners_cell[3]=-corners_cell[3];}
          }
          if (quad_neigh->y+dX[1]*2.<quad->y+dX[1])
          {
            corners_cell[1] = mesh->quad_to_quad[4*i+1]+1;
            if (quad_neigh->x<quad->x){corners_cell[1]=-corners_cell[1];}
          }
          if (quad_neigh->x<quad->x){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
        case 2:
        {
          if (quad_neigh->x+dX[0]*2.>quad->x+dX[0])
          {
            corners_cell[1] = mesh->quad_to_quad[4*i+2]+1;
            if (quad_neigh->y>quad->y){corners_cell[1]=-corners_cell[1];}
          }
          if (quad_neigh->x+dX[0]*2.<quad->x+dX[0])
          {
            corners_cell[0] = mesh->quad_to_quad[4*i+2]+1;
            if (quad_neigh->y>quad->y){corners_cell[0]=-corners_cell[0];}
          }
          if (quad_neigh->y>quad->y){neighbors[*Nneigh]=-neighbors[*Nneigh];}
        }
        break;
        case 3:
        {
          if (quad_neigh->x+dX[0]*2.>quad->x+dX[0])
          {
            corners_cell[3] = mesh->quad_to_quad[4*i+3]+1;
            if (quad_neigh->y<quad->y){corners_cell[3]=-corners_cell[3];}
          }
          if (quad_neigh->x+dX[0]*2.<quad->x+dX[0])
          {
            corners_cell[2] = mesh->quad_to_quad[4*i+3]+1;
            if (quad_neigh->y<quad->y){corners_cell[2]=-corners_cell[2];}
          }
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
        corners_cell[k-4] = index+1;
      }
      else
      {
        index = index-mesh->local_num_quadrants-mesh->ghost_num_quadrants;
        neigh = sc_array_index(mesh->corner_quad,index);
        neighbors[*Nneigh] = *neigh+1;
        corners_cell[k-4] = *neigh+1;
      }
      quad_neigh = sc_array_index(quadrants,neighbors[*Nneigh]-1);
      switch(k-4)
      {
        case 0:
        {
          if (quad_neigh->x>quad->x||quad_neigh->y>quad->y)
          {
            neighbors[*Nneigh]=-neighbors[*Nneigh];
            corners_cell[k-4] = -corners_cell[k-4];
          }
        }
        break;
        case 1:
        {
          if (quad_neigh->x<quad->x||quad_neigh->y>quad->y)
          {
            neighbors[*Nneigh]=-neighbors[*Nneigh];
            corners_cell[k-4] = -corners_cell[k-4];
          }
        }
        break;
        case 2:
        {
          if (quad_neigh->x>quad->x||quad_neigh->y<quad->y)
          {
            neighbors[*Nneigh]=-neighbors[*Nneigh];
            corners_cell[k-4] = -corners_cell[k-4];
          }
        }
        break;
        case 3:
        {
          if (quad_neigh->x<quad->x||quad_neigh->y<quad->y)
          {
            neighbors[*Nneigh]=-neighbors[*Nneigh];
            corners_cell[k-4] = -corners_cell[k-4];
          }
        }
        break;
      }
      *Nneigh=*Nneigh+1;
    }
  }
  *neighbors_out=neighbors;
  *corners_out=corners;
  *corners_cell_out=corners_cell;
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

static int coarsen_value (p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t** q)
{
  double X[3];
  p4est_topidx_t tt;
  var_t* data1;
  var_t* data2;
  var_t* data3;
  var_t* data4;
  int m;

  tt = p4est->first_local_tree;
  data1 = (var_t *) q[0]->p.user_data;
  data2 = (var_t *) q[1]->p.user_data;
  data3 = (var_t *) q[2]->p.user_data;
  data4 = (var_t *) q[3]->p.user_data;
  m=data1->coarsen+data2->coarsen+data3->coarsen+data4->coarsen;
  if (m>1) {return 1;}
  else {return 0;}
}

static int refine_value (p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* q)
{
  p4est_topidx_t tt;
  var_t* data;

  tt = p4est->first_local_tree;
  data = (var_t *) q->p.user_data;
  if (data->refine==1) {return 1;}
  else {return 0;}
}

static void replace_fn (p4est_t* p4est, p4est_topidx_t which_tree, int num_outgoing, p4est_quadrant_t* outgoing[], int num_incoming, p4est_quadrant_t* incoming[])
{
  var_t *parent_data,*child_data;
  int i,j,isol;
  double h;

  if (num_outgoing > 1)
  {
    /* this is coarsening */
    parent_data = (var_t *) incoming[0]->p.user_data;
    child_data = (var_t *) outgoing[0]->p.user_data;
    parent_data->u = (double*) malloc(sizeof(double)*child_data->nsol);
    parent_data->pol = (double*) malloc(sizeof(double)*child_data->nsol*child_data->size_pol);
    parent_data->nsol=child_data->nsol;
    parent_data->size_pol=child_data->size_pol;
    parent_data->refine=0;
    for (isol=0;isol<child_data->nsol;isol++)
    {
      parent_data->u[isol] = 0.;
      for (j=0;j<parent_data->size_pol;j++)
      {
        parent_data->pol[isol*parent_data->size_pol+j]=0.;
      }
    }
    for (i=0;i<P4EST_CHILDREN;i++)
    {
      child_data = (var_t *) outgoing[i]->p.user_data;
      for (isol=0;isol<child_data->nsol;isol++)
      {
        parent_data->u[isol] += child_data->u[isol]/P4EST_CHILDREN;
        for (j=0;j<parent_data->size_pol;j++)
        {
          parent_data->pol[isol*parent_data->size_pol+j] += child_data->pol[isol*parent_data->size_pol+j]/P4EST_CHILDREN;
        }
      }
      switch (i)
      {
        case 0:
          parent_data->u1 = child_data->u;
          break;
        case 1:
          parent_data->u2 = child_data->u;
          break;
        case 2:
          parent_data->u3 = child_data->u;
          break;
        case 3:
          parent_data->u4 = child_data->u;
          break;
      }
      free(child_data->u1);free(child_data->u2);free(child_data->u3);free(child_data->u4);
    }
  }
  else
  {
    /* this is refinement */
    parent_data = (var_t *) outgoing[0]->p.user_data;
    h=(double) P4EST_QUADRANT_LEN (outgoing[0]->level)/(double) P4EST_ROOT_LEN;
    for (i=0;i<P4EST_CHILDREN;i++)
    {
      child_data = (var_t *) incoming[i]->p.user_data;
      child_data->u = (double*) malloc(sizeof(double)*parent_data->nsol);
      child_data->u1 = (double*) malloc(sizeof(double)*parent_data->nsol);
      child_data->u2 = (double*) malloc(sizeof(double)*parent_data->nsol);
      child_data->u3 = (double*) malloc(sizeof(double)*parent_data->nsol);
      child_data->u4 = (double*) malloc(sizeof(double)*parent_data->nsol);
      child_data->pol = (double*) malloc(sizeof(double)*parent_data->nsol*parent_data->size_pol);
      child_data->nsol=parent_data->nsol;
      child_data->size_pol=parent_data->size_pol;
      child_data->coarsen=0;
      for (isol=0;isol<parent_data->nsol;isol++)
      {
        for (j=0;j<child_data->size_pol;j++)
        {
          child_data->pol[isol*child_data->size_pol+j]=parent_data->pol[isol*child_data->size_pol+j];
        }
      }
      switch (i)
      {
        case 0:
          for (isol=0;isol<parent_data->nsol;isol++)
          {
            child_data->u[isol] = parent_data->u1[isol];
            child_data->u1[isol] = parent_data->u1[isol];
            child_data->u2[isol] = parent_data->u1[isol];
            child_data->u3[isol] = parent_data->u1[isol];
            child_data->u4[isol] = parent_data->u1[isol];
          }
          break;
        case 1:
          for (isol=0;isol<parent_data->nsol;isol++)
          {
            child_data->u[isol] = parent_data->u2[isol];
            child_data->u1[isol] = parent_data->u2[isol];
            child_data->u2[isol] = parent_data->u2[isol];
            child_data->u3[isol] = parent_data->u2[isol];
            child_data->u4[isol] = parent_data->u2[isol];
          }
          break;
        case 2:
          for (isol=0;isol<parent_data->nsol;isol++)
          {
            child_data->u[isol] = parent_data->u3[isol];
            child_data->u1[isol] = parent_data->u3[isol];
            child_data->u2[isol] = parent_data->u3[isol];
            child_data->u3[isol] = parent_data->u3[isol];
            child_data->u4[isol] = parent_data->u3[isol];
          }
          break;
        case 3:
          for (isol=0;isol<parent_data->nsol;isol++)
          {
            child_data->u[isol] = parent_data->u4[isol];
            child_data->u1[isol] = parent_data->u4[isol];
            child_data->u2[isol] = parent_data->u4[isol];
            child_data->u3[isol] = parent_data->u4[isol];
            child_data->u4[isol] = parent_data->u4[isol];
          }
          break;
      }
    }
    free(parent_data->u);free(parent_data->u1);free(parent_data->u2);free(parent_data->u3);free(parent_data->u4);
    free(parent_data->pol);
  }
}

void p4_adapt (p4est_t* p4est, sc_array_t* quadrants, double* sol, double* sol_interp, double* pol_interp, int size_pol, int nsol, int* sol_coarsen, int* sol_refine, int maxlevel, int coarsen_recursive, int refine_recursive)
{
  int k;
  int i;
  int isol;
  p4est_quadrant_t* q;
  var_t* data;
  int callback_orphans=0;
  int nc=quadrants->elem_count;

  for (k=0;k<nc;k++)
  {
    q=sc_array_index(quadrants,k);
    data = q->p.user_data;
    data->u = (double*) malloc(sizeof(double)*nsol);
    data->u1 = (double*) malloc(sizeof(double)*nsol);
    data->u2 = (double*) malloc(sizeof(double)*nsol);
    data->u3 = (double*) malloc(sizeof(double)*nsol);
    data->u4 = (double*) malloc(sizeof(double)*nsol);
    data->pol = (double*) malloc(sizeof(double)*size_pol*nsol);
    data->nsol=nsol;
    data->size_pol=size_pol;
    data->coarsen=sol_coarsen[k];
    data->refine=sol_refine[k];

    for (isol=0;isol<nsol;isol++)
    {
      data->u[isol]=sol[isol*nc+k];
      data->u1[isol]=sol_interp[isol*4*nc+k];
      data->u2[isol]=sol_interp[isol*4*nc+k+nc];
      data->u3[isol]=sol_interp[isol*4*nc+k+2*nc];
      data->u4[isol]=sol_interp[isol*4*nc+k+3*nc];
      for (i=0;i<size_pol;i++)
      {
        data->pol[isol*size_pol+i]=pol_interp[isol*size_pol*nc+i*nc+k];
      }
    }
  }
  p4est_refine_ext(p4est,refine_recursive,maxlevel,refine_value,NULL,replace_fn);
  p4est_coarsen_ext(p4est,coarsen_recursive,callback_orphans,coarsen_value,NULL,replace_fn);
  p4est_balance_ext (p4est,P4EST_CONNECT_FULL,NULL,replace_fn);
}

void p4_new_sol (sc_array_t* quadrants, double** sol_out, double** pol_out)
{
  int nc;
  double* sol;
  double* pol;
  int k;
  int i;
  int isol;
  p4est_quadrant_t* q;
  var_t* data;

  nc=quadrants->elem_count;
  q=sc_array_index(quadrants,0);
  data=q->p.user_data;
  sol = (double*) malloc(sizeof(double)*data->nsol*nc);
  pol = (double*) malloc(sizeof(double)*data->nsol*nc*data->size_pol);

  for (k=0;k<nc;k++)
  {
    q=sc_array_index(quadrants,k);
    data=q->p.user_data;
    for (isol=0;isol<data->nsol;isol++)
    {
      sol[isol*nc+k]=data->u[isol];
      for (i=0;i<data->size_pol;i++)
      {
        pol[isol*data->size_pol*nc+i*nc+k]=data->pol[isol*data->size_pol+i];
      }
    }
    free(data->u);free(data->u1);free(data->u2);free(data->u3);free(data->u4);
    free(data->pol);
  }
  *sol_out=sol;
  *pol_out=pol;
}

void p4_free(void* ptr)
{
  free(ptr);
}

void p4_destroy(p4est_connectivity_t* conn, p4est_t* p4est)
{
  p4est_destroy (p4est);
  p4est_connectivity_destroy(conn);
}

void p4_destroy_mesh(p4est_mesh_t* mesh, p4est_nodes_t* nodes)
{
  p4est_mesh_destroy(mesh);
  p4est_nodes_destroy(nodes);
}
