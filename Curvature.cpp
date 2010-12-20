/* CGAL Delaunay Triangulation and GEOMVIEW */


/* #######################################################################################
   Code is for analysing a membrane surface, where the membrane comprises lipids which
   in turn are constructed from beads. The basic data structures should be simple enough
   to follow and bead positions are stored in arrays. Consider that a lipid comprises L
   beads and there are N lipids in the system. These data structures are mine and not to
   be found anyplace else, but as I said the principles can be easily understood.

   This is analysis code - written off the cuff and not designed to be robust or fast.
   Feel free to hack about and play with this code as I have done. It was originally
   intended to look at curvature of surfaces - in particular membranes, but the code should
   work for a variety of structures. Also this is a selection of the routines used, so do 
   not expect it to work, but it should give you some idea about how you might interact with
   CGAL.

   #######################################################################################

*/



#ifndef CGAL_USE_GEOMVIEW
int main()
{
  std::cout << "Geomview doesn't work on Windows, so..." << std::endl;
  return 0;
}
#else

#else

#include <fstream>
#include <unistd.h> 
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/intersections.h>

typedef CGAL::Cartesian<double>  K;
typedef K::Point_2 Point2;
typedef CGAL::Triangulation_euclidean_traits_xy_3<K> Gt3;
typedef Gt3::Point Point3;
typedef K::Segment_3 Segment3;
typedef K::Line_3 Line3;
typedef CGAL::Delaunay_triangulation_2<K>   Delaunay;
typedef CGAL::Delaunay_triangulation_2<Gt3> Terrain;
typedef CGAL::Delaunay_triangulation_3<K>   Delaunay3d;
typedef Terrain::Vertex_handle              Vertex_handle;
typedef Terrain::Face_handle                Face_handle;
typedef Terrain::Vertex_iterator            Vertex_iterator;
typedef Terrain::Face_iterator              Face_iterator;
typedef Terrain::Vertex                     Vertex;
typedef Terrain::Face                       Face;

void GrabSurface(int ***ptrLipidArray, double ***ptrBeadArray, int **ptrBeadType, int **ptrBeadSpecies, int MoleculeTotal, double Vertices[][3], int &NoVertices, int UpperLower, double Xbox, double Ybox, double Zbox, int Species[])
{
  int molecule, vcount, HeadBeadIndex, TailBeadIndex, loopspecies;
  double diff, zposHead, zposTail;
  vcount = 0;
  for (loopspecies = 0; loopspecies < NoSpecies; loopspecies++)
  {
    for (molecule=0;molecule<MoleculeTotal;molecule++)
    { 
      HeadBeadIndex = (*ptrLipidArray)[molecule][0];
      if ((*ptrBeadSpecies)[HeadBeadIndex] == Species[loopspecies])
      {
        
        TailBeadIndex = (*ptrLipidArray)[molecule][3];
        zposHead = (*ptrBeadArray)[HeadBeadIndex][2];
        zposTail = (*ptrBeadArray)[TailBeadIndex][2];
        diff = (zposHead - zposTail);
        
        if (((diff > 0) && (UpperLower == 0)) || ((diff < 0) && (UpperLower == 1)))
        {
          // Unfold loop
          Vertices[vcount][0] = (*ptrBeadArray)[HeadBeadIndex][0];
          Vertices[vcount][1] = (*ptrBeadArray)[HeadBeadIndex][1];
          Vertices[vcount][2] = (*ptrBeadArray)[HeadBeadIndex][2];
          vcount+=1;
        }
      }
    }
  }
  NoVertices = vcount;
}


void CreatePlane(int ***ptrLipidArray, double ***ptrBeadArray, int **ptrBeadType, int **ptrBeadSpecies, int MoleculeTotal, double Vertices[][3], int &NoVertices, int UpperLower, double Xbox, double Ybox, double Zbox, double Radius, double Length)
{
  int points_on_circum = 3;
  int points_on_line = 3;
  double sepCircum = Radius/points_on_circum;
  double sepLength = Length/points_on_line;
  int points, count, vertexcount;
  double x,y,z;
  vertexcount = 0;
  for (points=0;points<points_on_line;points++)
  {
    for (count=0;count<points_on_circum;count++)
    {
      z = 0.0;
      x = points*sepLength;
      y = count*sepCircum;
      Vertices[vertexcount][0] = x;
      Vertices[vertexcount][1] = y;
      Vertices[vertexcount][2] = z;
      vertexcount+=1;
    }
  } 
  NoVertices = vertexcount;
  Vertices[0][0] = 0.0;
  Vertices[0][1] = 0.0;
  Vertices[0][2] = 0.0;
  Vertices[1][0] = 10.0;
  Vertices[1][1] = 0.0;
  Vertices[1][2] = 0.0;
  Vertices[2][0] = 0.0;
  Vertices[2][1] = 10.0;
  Vertices[2][2] = 0.0;
  Vertices[3][0] = 10.0;
  Vertices[3][1] = 10.0;
  Vertices[3][2] = 0.0;
  NoVertices = 4;
}


void MulVector(double In[], double Out[], double scalar)
{
  Out[0] = In[0] * scalar;
  Out[1] = In[1] * scalar;
  Out[2] = In[2] * scalar;
}

double MagnitudeVector(double A[])
{
  return(sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]));
}


double DotProdvector(double A[], double B[])
{
  double Mag = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
  if (Mag > 1.0)
  {
    Mag = 1.0;
  }
  return(Mag);
}

void CrossProduct(double A[], double B[], double Out[])
{
  
  Out[0] = (A[1]*B[2] - A[2]*B[1]);
  Out[1] = -(A[0]*B[2] - A[2]*B[0]);
  Out[2] = (A[0]*B[1] - A[1]*B[0]);
}

void CalculateNormalToPlane(double VectorA[], double VectorB[], double Normal[], int UpperLower)
{
  // UpperLower gives the direction in the z - dimension. 0 = +ve, 1, = -ve
  double CrossProd[3];
  double MagCrossProd;
  CrossProduct(VectorA, VectorB, CrossProd);
  MagCrossProd = MagnitudeVector(CrossProd);
  MulVector(CrossProd, Normal, 1.0/MagCrossProd);
  if ((Normal[2] > 0) && (UpperLower == 1))
  {
    // flip the axes around
    Normal[0] = -Normal[0];
    Normal[1] = -Normal[1];
    Normal[2] = -Normal[2];
  }
  if ((Normal[2] <= 0) && (UpperLower == 0))
  {
    Normal[0] = -Normal[0];
    Normal[1] = -Normal[1];
    Normal[2] = -Normal[2];
  }
}

double AreaOfFace(double A[3], double B[3])
{
  double C[3];
  CrossProduct(A,B,C);
  return(fabs(magnitudevector(C)*0.5));
}


double SharedLength(double FaceA[][3], double FaceB[][3])
{
  double SharedVertex[2][3], Mag;
  int verticesobtained = 0, v0, v1;
  for (v0=0; v0<3; v0++)
  {
    for (v1=0; v1<3; v1++)
    {
      
      if ((FaceA[v0][0] == FaceB[v1][0]) && (FaceA[v0][1] == FaceB[v1][1]) && (FaceA[v0][2] == FaceB[v1][2]))
      {
        SharedVertex[verticesobtained][0] = FaceA[v0][0];
        SharedVertex[verticesobtained][1] = FaceA[v0][1];
        SharedVertex[verticesobtained][2] = FaceA[v0][2];
        verticesobtained+=1;
      }
    }
  }
  if ((SharedVertex[0][0] == SharedVertex[1][0]) && (SharedVertex[0][1] == SharedVertex[1][1]) && (SharedVertex[0][2] == SharedVertex[1][2]))
  {
    printf("Halt\n");
    getchar();
  }
  Mag = sqrt(SQR(SharedVertex[0][0]-SharedVertex[1][0]) + SQR(SharedVertex[0][1]-SharedVertex[1][1]) + SQR(SharedVertex[0][2]-SharedVertex[1][2]));
  return (Mag);
}

bool InOriginalMembrane(double FaceA[][3], double Xbox, double Ybox)
{
  // Lets do a check to see if Vertices lie within the original membrane. We dont want to overcount anything.
  bool OK=true;
  int vertex;
  for (vertex=0;vertex<3;vertex++)
  {
    if ((FaceA[vertex][0] > (Xbox/2.0)) || (FaceA[vertex][0] < (-Xbox/2.0)))
    {
      OK=false;
      //printf("failed...\n");
    }
    if ((FaceA[vertex][1] > (Ybox/2.0)) || (FaceA[vertex][1] < (-Ybox/2.0)))
    {
      OK=false;
    }

  }
  return(OK);
}

bool InOriginalMembrane(double FaceA[][3], double Xbox, double Ybox)
{
  // Lets do a check to see if Vertices lie within the original membrane. We dont want to overcount anything.
  bool OK=true;
  int vertex;
  for (vertex=0;vertex<3;vertex++)
  {
    if ((FaceA[vertex][0] > (Xbox/2.0)) || (FaceA[vertex][0] < (-Xbox/2.0)))
    {
      OK=false;
      //printf("failed...\n");
    }
    if ((FaceA[vertex][1] > (Ybox/2.0)) || (FaceA[vertex][1] < (-Ybox/2.0)))
    {
      OK=false;
    }

  }
  return(OK);
}

int main()
{
  Terrain T;
  JustGrabSurface(&LipidArray, &BeadArray, &BeadType, &BeadSpecies, TotalMolecules, Vertices, NoVertices, UpperLower, Xbox, Ybox,  Zbox, SpeciesSelection);
  for (vcount=0;vcount<NoVertices;vcount++)
  {
    T.insert( Point3(Vertices[vcount][0], Vertices[vcount][1], Vertices[vcount][2]) );
  }


  // This can be totally reworked below, rather hacky.
  if (1 == 1)
  {
    x_shift = Xbox;
    y_shift = 0.0;
    for (vcount=0;vcount<NoVertices;vcount++)
    {
      T.insert( Point3(Vertices[vcount][0]+x_shift, Vertices[vcount][1]+y_shift, Vertices[vcount][2]) );
    }
    x_shift = Xbox;
    y_shift = Ybox;
    for (vcount=0;vcount<NoVertices;vcount++)
    {
      T.insert( Point3(Vertices[vcount][0]+x_shift, Vertices[vcount][1]+y_shift, Vertices[vcount][2]) );
    }
    x_shift = 0.0;
    y_shift = Ybox;
    for (vcount=0;vcount<NoVertices;vcount++)
    {
     T.insert( Point3(Vertices[vcount][0]+x_shift, Vertices[vcount][1]+y_shift, Vertices[vcount][2]) );
    }
    x_shift = -Xbox;
    y_shift = 0.0;
    for (vcount=0;vcount<NoVertices;vcount++)
    {
      T.insert( Point3(Vertices[vcount][0]+x_shift, Vertices[vcount][1]+y_shift, Vertices[vcount][2]) );
    }
    x_shift = -Xbox;
    y_shift = -Ybox;
    for (vcount=0;vcount<NoVertices;vcount++)
    {
      T.insert( Point3(Vertices[vcount][0]+x_shift, Vertices[vcount][1]+y_shift, Vertices[vcount][2]) );
    }
    x_shift = 0.0;
    y_shift = -Ybox;
    for (vcount=0;vcount<NoVertices;vcount++)
    {
      T.insert( Point3(Vertices[vcount][0]+x_shift, Vertices[vcount][1]+y_shift, Vertices[vcount][2]) );
    }
    x_shift = Xbox;
    y_shift = -Ybox;
    for (vcount=0;vcount<NoVertices;vcount++)
    {
     T.insert( Point3(Vertices[vcount][0]+x_shift, Vertices[vcount][1]+y_shift, Vertices[vcount][2]) );
    }
    x_shift = -Xbox;
    y_shift = Ybox;
    for (vcount=0;vcount<NoVertices;vcount++)
    {
      T.insert( Point3(Vertices[vcount][0]+x_shift, Vertices[vcount][1]+y_shift, Vertices[vcount][2]) );
    }
  }
  if (1 == 2)
  {
    CGAL::Geomview_stream gv(CGAL::Bbox_3(-Xbox/2.0, -Ybox/2.0, -Zbox/2.0, Xbox/2.0, Ybox/2.0, Zbox/1.0));
    gv.set_bg_color(CGAL::Color(0, 200, 200));
    gv.set_line_width(10);
    gv.set_wired(true);     
    gv << T;
    sleep(50);
  }

  Vertex v0;
  Vertex v1;
  Vertex v2;
  Vertex v0_nn;
  Vertex v1_nn;
  Vertex v2_nn;
  Face face, neighbour_face;
  Face_handle neighbour;
  Face_iterator fit = T.faces_begin(), beyond = T.faces_end();
  Vertex_handle v_handle_0;
  Vertex_handle v_handle_1;
  Vertex_handle v_handle_2;
  Vertex_handle v_handle_0_nn;
  Vertex_handle v_handle_1_nn;
  Vertex_handle v_handle_2_nn;
  while (fit != beyond)
  {  
    face = *fit;
    ++fit;
    v_handle_0 = face.vertex(0);
    v_handle_1 = face.vertex(1);
    v_handle_2 = face.vertex(2);
    v0 = *v_handle_0;
    v1 = *v_handle_1;
    v2 = *v_handle_2;
    FaceA[0][0] = v0.point().x();
    FaceA[0][1] = v0.point().y();
    FaceA[0][2] = v0.point().z();
    FaceA[1][0] = v1.point().x();
    FaceA[1][1] = v1.point().y();
    FaceA[1][2] = v1.point().z();
    FaceA[2][0] = v2.point().x();
    FaceA[2][1] = v2.point().y();
    FaceA[2][2] = v2.point().z();
    VertexA[0] = FaceA[0][0];
    VertexA[1] = FaceA[0][1];
    VertexA[2] = FaceA[0][2];
    VertexB[0] = FaceA[1][0];
    VertexB[1] = FaceA[1][1];
    VertexB[2] = FaceA[1][2];
    VertexC[0] = FaceA[2][0];
    VertexC[1] = FaceA[2][1];
    VertexC[2] = FaceA[2][2];
    Centroid_Face(VertexA, VertexB, VertexC, CentroidA);
    if (InOriginalMembrane(FaceA, Xbox, Ybox))
    {
      LocalCurvaturePerFace = 0.0;
      A[0] = -v0.point().x() + v1.point().x();
      A[1] = -v0.point().y() + v1.point().y();
      A[2] = -v0.point().z() + v1.point().z();
      B[0] = -v0.point().x() + v2.point().x();
      B[1] = -v0.point().y() + v2.point().y();
      B[2] = -v0.point().z() + v2.point().z();
       
        areaofface = AreaOfFace(B,A);
        TotalArea+=areaofface;
       
        Mag = magnitudevector(A);
        Normalise(A, A, Mag);
        Mag = magnitudevector(B);
        Normalise(B, B, Mag);
        Mag = magnitudevector(A);
        CalculateNormalToPlane(A, B, Normal_base, UpperLower);
   
        count = 0;
        for (i=0;i<3;i++)
        {
          neighbour = face.neighbor(i);
          if (!T.is_infinite(neighbour))
          {
            neighbour_face = *neighbour;
            v_handle_0_nn = neighbour_face.vertex(0);
            v_handle_1_nn = neighbour_face.vertex(1);
            v_handle_2_nn = neighbour_face.vertex(2);
            v0_nn = *v_handle_0_nn;
            v1_nn = *v_handle_1_nn;
            v2_nn = *v_handle_2_nn;
            FaceB[0][0] = v0_nn.point().x();
            FaceB[0][1] = v0_nn.point().y();
            FaceB[0][2] = v0_nn.point().z();
            FaceB[1][0] = v1_nn.point().x();
            FaceB[1][1] = v1_nn.point().y();
            FaceB[1][2] = v1_nn.point().z();
            FaceB[2][0] = v2_nn.point().x();
            FaceB[2][1] = v2_nn.point().y();
            FaceB[2][2] = v2_nn.point().z();
            Points[i+1][0] = v0_nn.point().x();
            Points[i+1][1] = v0_nn.point().y();
            Points[i+1][2] = v0_nn.point().z();
            A_nn[0] = -v0_nn.point().x() + v1_nn.point().x();
            A_nn[1] = -v0_nn.point().y() + v1_nn.point().y();
            A_nn[2] = -v0_nn.point().z() + v1_nn.point().z();
            Mag = MgnitudeVector(A_nn);
            Normalise(A_nn, A_nn, Mag);   
            B_nn[0] = -v0_nn.point().x() + v2_nn.point().x();
            B_nn[1] = -v0_nn.point().y() + v2_nn.point().y();
            B_nn[2] = -v0_nn.point().z() + v2_nn.point().z();
            Mag = MagnitudeVector(B_nn);
            Normalise(B_nn, B_nn, Mag);
            CalculateNormalToPlane(A_nn, B_nn, Normal_neighbour, UpperLower);
            ScaleVectors(Normal_neighbour, scale_vec, Normal_neighbour_scaled);
            Points_off_plane[i+1][0] = Points[i+1][0] + Normal_neighbour_scaled[0];
            Points_off_plane[i+1][1] = Points[i+1][1] + Normal_neighbour_scaled[1];
            Points_off_plane[i+1][2] = Points[i+1][2] + Normal_neighbour_scaled[2];
            InverseCosine = DotProdVector(Normal_base, Normal_neighbour); 
            PartSum=PartSum + acos(inversecosine)*(SharedLength(FaceA, FaceB)/4.0);
            LocalCurvaturePerFace += (acos(inversecosine)*(SharedLength(FaceA, FaceB)/4.0));
            count+=1;
            paircount+=1;
            
          }
        }
        if (count == 3)
        {
          facecount+=1;
          Points[0][0] = CentroidA[0];
          Points[0][1] = CentroidA[1];
          Points[0][2] = CentroidA[2];
          ScaleVectors(Normal_base, scale_vec, Normal_base_scaled);
          Points_off_plane[0][0] = CentroidA[0] + Normal_base_scaled[0];
          Points_off_plane[0][1] = CentroidA[1] + Normal_base_scaled[1];
          Points_off_plane[0][2] = CentroidA[2] + Normal_base_scaled[2];
          mean_sep_on_plane = MeanSeparation(Points);
          
          mean_sep_off_plane = MeanSeparation(Points_off_plane);
       
          if ((mean_sep_off_plane > mean_sep_on_plane))
          {
            signage = +1.00;
            //printf("+ve\n");
          }else 
          {
            signage = -1.00;
            //printf("-ve\n");
          }
          if (UpperLower == 1)
          {
            signage = -signage;
          }

          // Histogram LocalCurvaturePerFace.
          // By local curvature...we mean "mean local curvature"
          LocalCurvaturePerFace = LocalCurvaturePerFace/areaofface;
          binwidth = (double) LocalCurvatureRange / (double) LocalCurvatureBins;
          bin = (int) ((signage * LocalCurvaturePerFace + offset) / binwidth);
          Histogram_Local_Curvature[bin]+=1;
          Hist_count+=1;
        
          Sum=Sum + 1.0*SQR(PartSum)*(1.0/areaofface);
          non_RMS_Sum += signage * LocalCurvaturePerFace;
 
          PartSum = 0.0;
          count = 0;
   
        }
      }
    }   
}
#endif






