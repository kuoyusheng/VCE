//
//  vib_cluster.h
//  orbital_space_v4
//
//  Created by kuoyusheng on 2015/9/21.
//  Copyright (c) 2015年 kuoyusheng. All rights reserved.
//

#ifndef orbital_space_v4_vib_cluster_h
#define orbital_space_v4_vib_cluster_h
#include "clusters.h"
class vib_cluster{
public:
    vib_cluster();
    vib_cluster(Cluster proper, Cluster imp)
    ~Cluster();
    //for cluster expansion
    Cluster proper;
    //for lattice dynamics
    Cluster imp;
    float* fct;
private:
    int  no_of_ver=proper.no_of_ver+imp.no_of_ver;
    int  no_of_pro_ver=proper.no_of_ver;
    int  no_of_imp_ver=
};

#endif
