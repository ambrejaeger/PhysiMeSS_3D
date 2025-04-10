#include "PhysiMeSS_vtk.h"

namespace PhysiCell {
   bool enable_vtk_saves = false;  // Declare the variable (extern)
   std::vector<std::vector<float>> colors = {
      {0.447f, 0.0, 0.216f},  // Wine Red
      {1.000f, 0.647f, 0.000f},   // Orange
      {1.000f, 1.000f, 1.000f},  // White
      {1.000f, 0.753f, 0.796f},   // Pink
      {0.439f, 0.0f, 0.0f}   // Prune
   };
   std::map<int, std::vector<float>> colorDict;
};

using namespace PhysiCell;
//Define enable_vtk_saves to determine if .vtp files should be saved

void read_save_vtk_status(pugi::xml_node config_root) {
   pugi::xml_node node;
   node = xml_find_node(config_root, "save");
   node = xml_find_node(node, "VTK");
   if (node) {
      enable_vtk_saves = xml_get_bool_value(node, "enable");
      std::cout << "Enable vtk save: " << enable_vtk_saves << std::endl;
   }
   else {
      std::cout << "VTK save not specified" << std::endl;
   }
}

void read_save_vtk_status(void) {
   return read_save_vtk_status(physicell_config_root);
}

// Function to generate color based on X_crosslink_count
std::vector<float> getColorForTube(float x_count) {

   float r = 0.0f;
   float g = 0.0f;
   float b = 1.0f;

   if (x_count == 0 ) {
      r = 240;
      g = 255;
      b = 255;
   }

   else if (x_count == 1 ) {
      r = 135;
      g = 206;
      b = 250;
   }

   else if (x_count == 2 ) {
      r = 0;
      g = 0;
      b = 205;
   }

   else {
      r = 0;
      g = 0;
      b = 139;
   }

   return {r/255, g/255, b/255};
}



bool vtp_save(std::string filename) {

   vtkNew<vtkMultiBlockDataSet> mbDataset;
    
   int blockIdx = 0;
   // Process cells
   for (Cell* pcell : *all_cells) {

      vtkSmartPointer<vtkPolyDataAlgorithm> source;
      std::vector<float> color;
       
      if (isFibre(pcell)) {
         PhysiMeSS_Fibre* pfibre = static_cast<PhysiMeSS_Fibre*>(pcell);

         // Create line source
         vtkNew<vtkLineSource> lineSource;
         lineSource->SetPoint1(
            pcell->position[0] - pfibre->mLength * pcell->state.orientation[0],
            pcell->position[1] - pfibre->mLength * pcell->state.orientation[1],
            pcell->position[2] - pfibre->mLength * pcell->state.orientation[2]
         );
         lineSource->SetPoint2(
            pcell->position[0] + pfibre->mLength * pcell->state.orientation[0],
            pcell->position[1] + pfibre->mLength * pcell->state.orientation[1],
            pcell->position[2] + pfibre->mLength * pcell->state.orientation[2]
         );

         // Create tube filter over the line
         vtkNew<vtkTubeFilter> tubeFilter;
         tubeFilter->SetInputConnection(lineSource->GetOutputPort());
         tubeFilter->SetRadius(pfibre->mRadius);
         tubeFilter->SetNumberOfSides(30);
         tubeFilter->CappingOn();
         source = tubeFilter;

         // Get color for tube
         color = getColorForTube(std::min(3, pfibre->X_crosslink_count));
      } 
      else {
         // Create sphere source
         vtkNew<vtkSphereSource> sphereSource;
         sphereSource->SetCenter(pcell->position[0], pcell->position[1], pcell->position[2]);
         sphereSource->SetRadius(pcell->phenotype.geometry.radius);
         sphereSource->SetPhiResolution(50);
         sphereSource->SetThetaResolution(50);
         source = sphereSource;

         // Get color for sphere
         if (PhysiCell::colorDict.find(pcell->type) != PhysiCell::colorDict.end()) {
               color = PhysiCell::colorDict[pcell->type];
         } 
         else {
            color = PhysiCell::colors.back();
            PhysiCell::colorDict[pcell->type] = color;
            PhysiCell::colors.pop_back();
         }
      }
   

   source->Update();
   vtkSmartPointer<vtkPolyData> blockData = source->GetOutput();

   // Create color array for THIS BLOCK
   vtkNew<vtkFloatArray> blockColors;
   blockColors->SetName("Color");
   blockColors->SetNumberOfComponents(3);
   mbDataset->SetBlock(blockIdx++, source->GetOutput());

   // Fill colors for all cells in this block
   for (vtkIdType i = 0; i < blockData->GetNumberOfCells(); i++) {
      blockColors->InsertNextTuple3(color[0], color[1], color[2]);
   }

   // Attach colors to THIS BLOCK's data
   blockData->GetCellData()->AddArray(blockColors);
   }
   // Write the multiblock file
   vtkNew<vtkXMLMultiBlockDataWriter> writer;
   writer->SetFileName(filename.c_str());
   writer->SetInputData(mbDataset);
   writer->Write();
   
   return true;
}

