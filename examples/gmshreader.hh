// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_GMSHREADER_HH
#define DUNE_GMSHREADER_HH

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<stdio.h>

#include<dune/common/geometrytype.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/grid/common/boundarysegment.hh>

namespace Dune {

  // Options for read operation
  struct GmshReaderOptions
  {
    enum GeometryOrder {
      /** @brief edges are straight lines. */
      firstOrder,
      /** @brief quadratic boundary approximation. */
      secondOrder 
    };
  };

  // arbitrary dimension, implementation is in specialization
  template<typename GridType, int dimension>
  class GmshReaderImp
  {
  };

  template<typename GridType>
  class GmshReaderImp<GridType,2>
  {
    // quadratic boundary segments in 1d
    class QuadraticBoundarySegment : public Dune::BoundarySegment<2>
    {
    public:
      QuadraticBoundarySegment (Dune::FieldVector<double,2> p0_, Dune::FieldVector<double,2> p1_, 
                                Dune::FieldVector<double,2> p2_)
        : p0(p0_), p1(p1_), p2(p2_)
      {}

      virtual Dune::FieldVector<double,2> operator() (const Dune::FieldVector<double,1>& local) const
      {
        Dune::FieldVector<double,2> y;
        y = 0.0;
        y.axpy(2.0*(local[0]-0.5)*(local[0]-1.0),p0);
        y.axpy(1.0-4.0*(local[0]-0.5)*(local[0]-0.5),p1);
        y.axpy(2.0*local[0]*(local[0]-0.5),p2);
        return y;
      }

    private:
      Dune::FieldVector<double,2> p0,p1,p2;
    };

    // quadratic boundary segments in 1d
    class LinearBoundarySegment : public Dune::BoundarySegment<2>
    {
    public:
      LinearBoundarySegment (Dune::FieldVector<double,2> p0_, Dune::FieldVector<double,2> p1_)
        : p0(p0_), p1(p1_)
      {}

      virtual Dune::FieldVector<double,2> operator() (const Dune::FieldVector<double,1>& local) const
      {
        Dune::FieldVector<double,2> y;
        y = 0.0;
        y.axpy(1.0-local[0],p0);
        y.axpy(local[0],p1);
        return y;
      }

    private:
      Dune::FieldVector<double,2> p0,p1;
    };

  public:

    static GridType* read(Dune::GridFactory<GridType>& factory, const std::string& fileName, 
                          GmshReaderOptions::GeometryOrder order, bool verbose)
    { 
      // check dimension
      const int dim = GridType::dimension;
      if (dim != 2)
        DUNE_THROW(Dune::NotImplemented, 
                   "Reading Gmsh format is not yet implemented for dimension " << dim);

      // open file name, we use C I/O
      FILE* file = fopen(fileName.c_str(),"r");
      if (file==0)
        DUNE_THROW(Dune::IOError, "Could not open " << fileName);

      // a read buffer
      char buf[512];

      //=========================================
      // Pass 1: Read vertices into vector
      //         Check vertices that are needed
      //         Renumber needed vertices
      //=========================================

      // process header
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      int version_number, file_type, data_size;
      fscanf(file,"%d %d %d\n",&version_number,&file_type,&data_size);
      if (version_number!=2)
         DUNE_THROW(Dune::IOError, "can only read version_number==2");
      if (verbose) std::cout << "version 2 Gmsh file detected" << std::endl; 
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");
      
      // node section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      int number_of_nodes;
      fscanf(file,"%d\n",&number_of_nodes);
      if (verbose) std::cout << "file contains " << number_of_nodes << " nodes" << std::endl;
      std::vector<Dune::FieldVector<double,dim> > nodes(number_of_nodes+1); // store positions
      for (int i=1; i<=number_of_nodes; i++)
        {
          int id;
          double x,y,z;
          fscanf(file,"%d %lg %lg %lg\n",&id,&x,&y,&z);
          //          if (verbose) std::cout << id << " " << x << " " << y << " " << z << std::endl;
          if (id!=i)
            DUNE_THROW(Dune::IOError, "expected id " << i);

          Dune::FieldVector<double,dim> position;
          position[0] = x; position[1] = y;
          nodes[i] = position; // just store node position
        }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndNodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndNodes");
  
      // element section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      int number_of_elements;
      fscanf(file,"%d\n",&number_of_elements);
      if (verbose) std::cout << "file contains " << number_of_elements << " elements" << std::endl; 
      unsigned int number_of_real_vertices=0;  // count number of vertices that are really needed
      std::map<int,unsigned int> renumber;
      for (int i=1; i<=number_of_elements; i++)
        {
          int id, elm_type, number_of_tags;
          fscanf(file,"%d %d %d",&id,&elm_type,&number_of_tags);
          int physical_entity, elementary_entity, mesh_partition;
          for (int k=1; k<=number_of_tags; k++)
            {
              int blub;
              fscanf(file,"%d",&blub);
              if (k==1) physical_entity = blub;
              if (k==2) elementary_entity = blub;
              if (k==3) mesh_partition = blub;
            }
          std::vector<int> simplexVertices(6);
          switch (elm_type)
            {
            case 1: // 2-node line
              simplexVertices.resize(2);
              fscanf(file,"%d %d\n",&(simplexVertices[0]),&(simplexVertices[1]));
              for (int i=0; i<2; i++)
                if (renumber.find(simplexVertices[i])==renumber.end()) 
                  {
                    //                   std::cout << "real vertex " << number_of_real_vertices << " at " << nodes[simplexVertices[i]] << " old number was " << simplexVertices[i] << std::endl;
                    renumber[simplexVertices[i]] = number_of_real_vertices++;
                    factory.insertVertex(nodes[simplexVertices[i]]);                  
                 }
              break;
            case 8: // 3-node line
              simplexVertices.resize(3);
              fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
              for (int i=0; i<2; i++)
                if (renumber.find(simplexVertices[i])==renumber.end()) 
                  {
                    renumber[simplexVertices[i]] = number_of_real_vertices++;
                    factory.insertVertex(nodes[simplexVertices[i]]);                  
                  }
              break;
            case 2: // 3-node triangle
              simplexVertices.resize(3);
              fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
              for (int i=0; i<3; i++)
                if (renumber.find(simplexVertices[i])==renumber.end()) 
                  {
                    renumber[simplexVertices[i]] = number_of_real_vertices++;
                    factory.insertVertex(nodes[simplexVertices[i]]);                  
                  }
              break;
            case 9: // 6-node triangle
              simplexVertices.resize(6);
              fscanf(file,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                                    &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]));
              for (int i=0; i<3; i++)
                if (renumber.find(simplexVertices[i])==renumber.end()) 
                  {
                    renumber[simplexVertices[i]] = number_of_real_vertices++;
                    factory.insertVertex(nodes[simplexVertices[i]]);                  
                  }
              break;
            default:            
              fgets(buf,512,file); // skip rest of line if no tetrahedron
            }
        }
      if (verbose) std::cout << "number of real vertices = " << number_of_real_vertices << std::endl;
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");

      //==============================================
      // Pass 2: Insert boundary segments and elements
      //==============================================

      // go to beginning of file
      rewind(file);

      // process header
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      fscanf(file,"%d %d %d\n",&version_number,&file_type,&data_size);
      if (version_number!=2)
         DUNE_THROW(Dune::IOError, "can only read version_number==2");
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");
      
      // node section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      fscanf(file,"%d\n",&number_of_nodes);
      for (int i=1; i<=number_of_nodes; i++)
        {
          int id;
          double x,y,z;
          fscanf(file,"%d %lg %lg %lg\n",&id,&x,&y,&z);
          if (id!=i)
            DUNE_THROW(Dune::IOError, "expected id " << i);
        }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndNodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndNodes");
  
      // element section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      fscanf(file,"%d\n",&number_of_elements);
      for (int i=1; i<=number_of_elements; i++)
        {
          int id, elm_type, number_of_tags;
          fscanf(file,"%d %d %d",&id,&elm_type,&number_of_tags);
          int physical_entity, elementary_entity, mesh_partition;
          for (int k=1; k<=number_of_tags; k++)
            {
              int blub;
              fscanf(file,"%d",&blub);
              if (k==1) physical_entity = blub;
              if (k==2) elementary_entity = blub;
              if (k==3) mesh_partition = blub;
            }
          std::vector<int> simplexVertices(6);
          std::vector<unsigned int> vertices(3);
          switch (elm_type)
            {
            case 1: // 2-node line
              simplexVertices.resize(2);
              fscanf(file,"%d %d\n",&(simplexVertices[0]),&(simplexVertices[1]));
              vertices.resize(2);
              for (int i=0; i<2; i++)
                vertices[i] = renumber[simplexVertices[i]]; // renumber vertices
              factory.insertBoundarySegment(vertices,
                                            new LinearBoundarySegment(nodes[simplexVertices[0]],
                                                                      nodes[simplexVertices[1]]));
              break;
            case 8: // 3-node line
              simplexVertices.resize(3);
              fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
              vertices.resize(2);
              for (int i=0; i<2; i++)
                vertices[i] = renumber[simplexVertices[i]]; // renumber vertices
              factory.insertBoundarySegment(vertices,
                                            new QuadraticBoundarySegment(nodes[simplexVertices[0]],
                                                                         nodes[simplexVertices[2]],
                                                                         nodes[simplexVertices[1]]));
              break;
            case 2: // 3-node triangle
              simplexVertices.resize(3);
              fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
              vertices.resize(3);
              for (int i=0; i<3; i++)
                vertices[i] = renumber[simplexVertices[i]]; // renumber vertices
              factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
              break;
            case 9: // 6-node triangle
              simplexVertices.resize(6);
              fscanf(file,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                                    &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]));
              vertices.resize(3);
              for (int i=0; i<3; i++)
                vertices[i] = renumber[simplexVertices[i]]; // renumber vertices
              factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
              break;
            default:            
              fgets(buf,512,file); // skip rest of line if no tetrahedron
            }
        }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");

      fclose(file);

      return factory.createGrid();
    }
  };

  template<typename GridType>
  class GmshReaderImp<GridType,3>
  {
  public:
    static GridType* read(Dune::GridFactory<GridType>& factory, const std::string& fileName, 
                          GmshReaderOptions::GeometryOrder order, bool verbose)
    { 
      // check dimension
      const int dim = GridType::dimension;
      if (dim != 3)
        DUNE_THROW(Dune::NotImplemented, 
                   "Reading Gmsh format is not yet implemented for dimension " << dim);

      // open file name, we use C I/O
      FILE* file = fopen(fileName.c_str(),"r");
      if (file==0)
        DUNE_THROW(Dune::IOError, "Could not open " << fileName);

      // a read buffer
      char buf[512];

      // process header
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      int version_number, file_type, data_size;
      fscanf(file,"%d %d %d\n",&version_number,&file_type,&data_size);
      if (version_number!=2)
         DUNE_THROW(Dune::IOError, "can only read version_number==2");
      if (verbose) std::cout << "version 2 Gmsh file detected" << std::endl; 
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");
      
      // node section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      int number_of_nodes;
      fscanf(file,"%d\n",&number_of_nodes);
      if (verbose) std::cout << "file contains " << number_of_nodes << " nodes" << std::endl; 
      for (int i=1; i<=number_of_nodes; i++)
        {
          int id;
          double x,y,z;
          fscanf(file,"%d %lg %lg %lg\n",&id,&x,&y,&z);
          //          if (verbose) std::cout << id << " " << x << " " << y << " " << z << std::endl;
          if (id!=i)
            DUNE_THROW(Dune::IOError, "expected id " << i);

          Dune::FieldVector<double,dim> position;
          position[0] = x; position[1] = y; position[2] = z;
          factory.insertVertex(position);
        }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndNodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndNodes");
  
      // element section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      int number_of_elements;
      fscanf(file,"%d\n",&number_of_elements);
      if (verbose) std::cout << "file contains " << number_of_elements << " elements" << std::endl; 
      int number_of_tetrahedra=0;
      for (int i=1; i<=number_of_elements; i++)
        {
          int id, elm_type, number_of_tags;
          fscanf(file,"%d %d %d",&id,&elm_type,&number_of_tags);
          if (elm_type!=4)
            fgets(buf,512,file); // skip rest of line if no tetrahedron
          else
            {
              int physical_entity, elementary_entity, mesh_partition;
              for (int k=1; k<=number_of_tags; k++)
                {
                  int blub;
                  fscanf(file,"%d",&blub);
                  if (k==1) physical_entity = blub;
                  if (k==2) elementary_entity = blub;
                  if (k==3) mesh_partition = blub;
                }
              std::vector<unsigned int> simplexVertices(4);
              fscanf(file,"%d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),
                                          &(simplexVertices[2]),&(simplexVertices[3]));
              for (int k=0; k<4; k++) simplexVertices[k] -= 1; 
              factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim), 
                                    simplexVertices);
            }
        }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");

      fclose(file);

      return factory.createGrid();
    }
  };

  template<typename GridType>
  class GmshReader
  {
  public:
    static GridType* read (const std::string& fileName, 
                           GmshReaderOptions::GeometryOrder order = GmshReaderOptions::firstOrder,
                           bool verbose = true)
    { 
      // make a grid factory
      Dune::GridFactory<GridType> factory;

      return GmshReaderImp<GridType,GridType::dimension>::read(factory,fileName,order,verbose);
    }

    static GridType* read (GridType& grid, const std::string& fileName, 
                           GmshReaderOptions::GeometryOrder order = GmshReaderOptions::firstOrder,
                           bool verbose = true)
    { 
      // make a grid factory
      Dune::GridFactory<GridType> factory(&grid);

      return GmshReaderImp<GridType,GridType::dimension>::read(factory,fileName,order,verbose);
    }
  };


} // namespace Dune

#endif
