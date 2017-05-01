#include<opencv2/highgui/highgui.hpp>
#include<opencv2/core/core.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
using namespace cv;

void augmentFrames(std::string filename, float label, std::ofstream &outfile);

int main()
{
    int omp_procs = omp_get_num_procs();
    printf("Processors: %d\n", omp_procs);
    omp_set_num_threads(omp_procs);
    //omp_set_num_threads(1);

    std::ifstream inLabels; 
    std::string data, filenames[1000], labelstr, filename;
    std::ofstream outfile ("Results/labels.txt",  std::ofstream::out);
    float labels[1000], label;
    inLabels.open("Set_125_offset_v/labels.txt");
    int count =0;
    
    do{
	inLabels >> data;
        filenames[count] = data.substr(0,data.find(":"));

	 // intermidiate parse step 
        labelstr = data.substr(data.find(":")+1, data.find(","));
        labelstr = labelstr.substr(0, labelstr.find(","));
        labels[count] = stof(labelstr);
        std::cout << filename << "     " << label << "\n";
	count +=1;
	    
    }while(!inLabels.eof());

    #pragma omp parallel for 
    for(int i=0; i<count;i++)
    {
	augmentFrames(filenames[i], labels[i], outfile);
    }

    outfile.close(); 
    return 0;
}


void augmentFrames(std::string filename, float label, std::ofstream &outfile)
{
    int dst_width=320, dst_height=180, warp_factor=20;
    float mirrored_label, lw1, rw1, lw2, rw2;
    std::string filedir; 
    
    if (outfile == NULL)
	    std::cout << "error\n";

    filedir = "Set_125_offset_v/Frame_" + filename + ".jpg";
    Mat img = imread(filedir,CV_LOAD_IMAGE_COLOR);
    
    Mat mirrored;
    Size mySize = Size(dst_width, dst_height);
    resize(img, img, mySize);

    flip(img, mirrored, 1);
    std::cout << "flipped label: " << label * -1 << "\n";

    imwrite("Results/Frame_" + filename  + "_2.jpg", mirrored);
  
    Rect crop(0, warp_factor, dst_width, dst_height-warp_factor*2);
    Point2f pts1[4], pts2[4];

    pts1[0] = Point2f(0, 0);
    pts1[1] = Point2f(320, 0);
    pts1[2] = Point2f(0, dst_height);
    pts1[3] = Point2f(320, dst_height);

    pts2[0] = Point2f(0, warp_factor);
    pts2[1] = Point2f(320, 0);
    pts2[2] = Point2f(0, dst_height-warp_factor);
    pts2[3] = Point2f(320, dst_height);

    Mat Transform = getPerspectiveTransform(pts1, pts2);
    Mat warped1, warped2;
    warpPerspective(img, warped1, Transform, mySize);
    warped1 = warped1(crop);
    resize(warped1, warped1, mySize);
    imwrite("Results/Frame_" + filename  + "_0.jpg", warped1);

    pts2[0] = Point2f(0, 0);
    pts2[1] = Point2f(320, warp_factor);
    pts2[2] = Point2f(0, dst_height);
    pts2[3] = Point2f(320, dst_height-warp_factor);

    Transform = getPerspectiveTransform(pts1, pts2);
    warpPerspective(img, warped2, Transform, mySize);
    warped2 = warped2(crop);
    resize(warped2, warped2, mySize);
    imwrite("Results/Frame_" + filename + "_1.jpg", warped2);

    lw1 = label + (warp_factor*.01);
    rw1 = label - (warp_factor*.01);

    /////////// Change warp factor and compute 2 new frames ///////////
    warp_factor = 40;
    Rect crop2(0, warp_factor, dst_width, dst_height-warp_factor*2);
    
    pts2[0] = Point2f(0, warp_factor);
    pts2[1] = Point2f(320, 0);
    pts2[2] = Point2f(0, dst_height-warp_factor);
    pts2[3] = Point2f(320, dst_height);

    Transform = getPerspectiveTransform(pts1, pts2);
    warpPerspective(img, warped1, Transform, mySize);
    warped1 = warped1(crop2);
    resize(warped1, warped1, mySize);
    imwrite("Results/Frame_" + filename  + "_3.jpg", warped1);

    pts2[0] = Point2f(0, 0);
    pts2[1] = Point2f(320, warp_factor);
    pts2[2] = Point2f(0, dst_height);
    pts2[3] = Point2f(320, dst_height-warp_factor);

    Transform = getPerspectiveTransform(pts1, pts2);
    warpPerspective(img, warped2, Transform, mySize);
    warped2 = warped2(crop2);
    resize(warped2, warped2, mySize);
    imwrite("Results/Frame_" + filename + "_4.jpg", warped2);

    lw2 = label + (warp_factor*.0185);
    rw2 = label - (warp_factor*.0185);

    ////////////// write to file ///////////////////////////////////	    

    std::cout << "label: " << label << "\n";
    std::cout << "left 1 warp label: " << lw1 << "\n";
    std::cout << "right 1 warp label: " << rw1 << "\n";
    std::cout << "left 2 warp label: " << lw2 << "\n";
    std::cout << "right 2 warp label: " << rw2 << "\n";

    // write label to file
    outfile << filename + "_0:" + std::to_string(lw1) + "\n";
    outfile << filename + "_1:" + std::to_string(rw1) + "\n";
    outfile << filename + "_2:" + std::to_string(label*-1) + "\n";
    outfile << filename + "_3:" + std::to_string(lw2) + "\n";
    outfile << filename + "_4:" + std::to_string(rw2) + "\n";
}


