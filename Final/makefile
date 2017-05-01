make: image_data_augmentation_parallel.c 
	g++ -ggdb `pkg-config --cflags opencv` -o `basename image_data_augmentation_parallel.c .c` image_data_augmentation_parallel.c `pkg-config --libs opencv` -std=c++11 -fopenmp

