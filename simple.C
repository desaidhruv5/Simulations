#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen/Eigen/Dense"
#include <cstdio>
#include <time.h>






int main()
{
  std::vector<int> v;
  v.push_back(5);

  int tot =1000000;
  for (int i = 0; i< tot; ++i)
  {
    std::vector<int> next_v;

    next_v.push_back(i);
    next_v.push_back(i+1);
    v = next_v;
    std::cout << "Percent complete: "<< float(i)/float(tot)*100 <<", size of next_v "<<float(next_v.size()) << ", and "<< next_v[0]<< "\r"<<std::flush;
    if (i==999999)
    {
      std::cout<<std::endl;
      std::cout<<"Done"<<std::endl;
    }
  }


  return 0;
}


