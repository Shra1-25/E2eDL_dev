#include "E2eDL/FrameProducers/interface/FrameCropping.h"

// Cropping frames with given width and height: frame_width, frame_height
void e2e::getFrame ( e2e::Frame2D& vframe, const e2e::seed& objSeed, const e2e::Frame2D* vdetector_image,
  int detImg_height, int detImg_width ) // can be removed if using Frame2D input
{

  //std::vector<std::vector<float>> vframe = std::vector<std::vector<float>> (frame_height,std::vector<float> (frame_width, 0.0));
  //std::cout << " >> vdetimg:" << vdetector_image->size() << "," << (*vdetector_image)[0].size() << std::endl;
  //std::cout << " >> vdetimg:" << vdetector_image->size() << std::endl;
  //std::cout << " >> vframe:" << vframe.size() << "," << vframe[0].size() << std::endl;
  //std::cout << " >> vframe:" << vframe.size() << "," << vframe.front().size() << std::endl;
  //std::cout << " >> vframe:" << vframe.size() << std::endl;
  int frame_height = vframe.size();
  int frame_width  = vframe.front().size();
  int ieta_seed = objSeed[0];
  int iphi_seed = objSeed[1];

  // Determining start and end indices for the frames to be cropped
  int start_x=0;
  int end_x=0;
  int start_y=0;
  //int end_y=0;
  int buff_x=0;
  int buff_y=0;
  int half_frame_width;
  int half_frame_height;
  if (frame_width%2==0){half_frame_width=frame_width/2-1;} else {half_frame_width=frame_width/2;}
  if(frame_height%2==0){half_frame_height=frame_height/2-1;} else {half_frame_height=frame_height/2;}
  if (iphi_seed<half_frame_width-1){
   start_y=0;
   buff_y=half_frame_width-iphi_seed;
   buff_y=detImg_width-buff_y;
  }
  else {
   start_y=iphi_seed-half_frame_width;
   buff_y=0;
  }
  if (iphi_seed>detImg_width-half_frame_width){
   //end_y=detImg_width-1;
  }
  else {
   //end_y=iphi_seed+frame_width/2;
  }
  if (ieta_seed<half_frame_height){
   start_x=0;
   buff_x=half_frame_height-ieta_seed;
  }
  else {
   start_x=ieta_seed-half_frame_height;
    buff_x=0;
  }
  if (ieta_seed>detImg_height-half_frame_height){
   end_x=detImg_height-1;
  }
  else {
   end_x=ieta_seed+half_frame_height;
  }

  // Cropping the flat input vector based on the above determined start and end indices
  for (int x_idx = start_x; x_idx<=end_x;x_idx++){
   for (int y_idx = 0/*start_y*/; y_idx<frame_width/*=end_y*/;y_idx++){
    vframe[x_idx-start_x+buff_x][y_idx/*y_idx-start_y+buff_y*/] = (*vdetector_image)[x_idx][(y_idx+buff_y+start_y)];
    //vEB_flat_frame[(x_idx-start_x+buff_x)*frame_width+y_idx/*-start_y+buff_y*/]=vdetector_image[x_idx*detImg_width+(y_idx+start_y+buff_y)%detImg_width];
    //std::cout<<"("<<x_idx-start_x+buff_x<<","<<y_idx<<"): "<<vframe[x_idx-start_x+buff_x][y_idx/*y_idx-start_y+buff_y*/]<<" "<<vdetector_image[x_idx*detImg_width+(y_idx+start_y+buff_y)%detImg_width];
   }
   //std::cout<<std::endl;
  }
  //std::cout<<" >> Size of frame is:"<<"("<<vframe.size()<<", "<<vframe[0].size()<<")"<<std::endl;
  //std::cout<<" >> Detector Image value at ("<<ieta_seed<<", "<<iphi_seed<<")is: "<<vframe[half_frame_height][half_frame_width]<<std::endl;
  return;
}
