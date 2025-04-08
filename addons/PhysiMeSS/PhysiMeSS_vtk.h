#ifndef __PhysiMeSS_vtk__
#define __PhysiMeSS_vtk__

#include <iostream>
#include <fstream>
#include <string>

#include <vtkSmartPointer.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkNamedColors.h> 

#include <vtkLine.h>
#include <vtkVertex.h>


#include <vtkNew.h>          
#include <vtkLineSource.h>   
#include <vtkTubeFilter.h>   
#include <vtkSphereSource.h> 
#include <vtkAppendPolyData.h> 
#include <vtkXMLPolyDataWriter.h> 
#include <vtkPolyData.h>     
#include <vtkCellArray.h>    

#include "../../core/PhysiCell.h"
#include "PhysiMeSS_fibre.h"
#include "../../BioFVM/pugixml.hpp"  

namespace PhysiCell {
    extern bool enable_vtk_saves;  
};
void read_save_vtk_status(pugi::xml_node config_root);
void read_save_vtk_status(void);
std::vector<float> getColorForTube(float x_count);
bool vtp_save(std::string filename);

#endif