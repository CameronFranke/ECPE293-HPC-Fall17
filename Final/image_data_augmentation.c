#include<opencv2/highgui/highgui.hpp>
#include<opencv2/core/core.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <fstream>
#include <string>
using namespace cv;

void augmentFrames(std::string filename, float label, std::ofstream &outfile);

int main()
{
    
    std::ifstream inLabels; 
    std::string data, filename, labelstr;
    std::ofstream outfile ("Results/labels.txt",  std::ofstream::out);
    float label;
    inLabels.open("Set_125_offset_v/labels.txt");
    
    do{
	//read in and parse labels 
     	inLabels >> data;
	filename = data.substr(0,data.find(":"));
	std::cout << filename << "     ";

	// intermidiate parse step 
	labelstr = data.substr(data.find(":")+1, data.find(","));
	labelstr = labelstr.substr(0, labelstr.find(","));
	label = stof(labelstr);

	std::cout << label << "\n";
	
	//function call to augment data
	augmentFrames(filename, label, outfile);

    }while(!inLabels.eof());

    outfile.close(); 
    return 0;
}


void augmentFrames(std::string filename, float label, std::ofstream &outfile)
{
    int dst_width=320, dst_height=180, warp_factor=20;
    float mirrored_label, lw1, rw1;
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

    std::cout << "label: " << label << "\n";
    std::cout << "left warp label: " << lw1 << "\n";
    std::cout << "right warp label: " << rw1 << "\n";

    // write label to file
    outfile << filename + "_0:" + std::to_string(lw1) + "\n";
    outfile << filename + "_1:" + std::to_string(rw1) + "\n";
    outfile << filename + "_2:" + std::to_string(label*-1) + "\n";
}


