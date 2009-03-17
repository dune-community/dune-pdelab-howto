// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_GMSHREADER_HH
#define DUNE_GMSHREADER_HH

#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>

#include<dune/common/geometrytype.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/referenceelements.hh>

namespace Dune {

  template<typename GridType>
  class GmshReader
  {
  public:
    static GridType* read (const std::string& fileName, bool verbose = true)
    { 
      // make a grid factory
      Dune::GridFactory<GridType> factory;

      return read(factory,fileName,verbose);
    }

    static GridType* read (GridType& grid, const std::string& fileName, bool verbose = true)
    { 
      // make a grid factory
      Dune::GridFactory<GridType> factory(&grid);

      return read(factory,fileName,verbose);
    }

  private:

    static GridType* read (Dune::GridFactory<GridType>& factory, const std::string& fileName, 
                           bool verbose = true)
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

      return factory.createGrid();
    }
  };

} // namespace Dune

#endif
