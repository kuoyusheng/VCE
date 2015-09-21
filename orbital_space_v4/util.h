//
//  util.h
//  orbital_space_v4
//
//  Created by kuoyusheng on 2015/9/21.
//  Copyright (c) 2015年 kuoyusheng. All rights reserved.
//

#ifndef orbital_space_v4_util_h
#define orbital_space_v4_util_h
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include </usr/include/Eigen/Eigen> // works ok
#include "null_space.h"
#include "clusters.h"
using namespace std;
using Eigen::MatrixXd;
int power(int val, int pow){
    int result=1;
    if(pow==0)
        return result;
    
    for (int i=0; i<(pow); i++) {
        result*=val;
    }
    return result;
}

void digit_combination(int dim,int i,int arr[],int digit)
{
    if (dim / 3 != 0) {
        digit_combination((dim / 3),i+1,arr,digit);
    }
    arr[digit-1-i]=(dim % 3);
}


void tensor_ind(int **tensor,int tensor_dim, int no_of_ver){
    for (int i=0; i<tensor_dim; i++) {
        digit_combination(i,0,tensor[i],no_of_ver);
    }
}
void digit_combination_2(int dim,int i,int arr[],int digit, int cut_off)
{
    if (dim / cut_off != 0) {
        digit_combination_2((dim / cut_off),i+1,arr,digit,cut_off);
    }
    arr[digit-1-i]=(dim % cut_off);
}

void find_np(int **tensor,int tensor_dim, int cut_off){
    for (int i=0; i<tensor_dim; i++) {
        digit_combination_2(i,0,tensor[i],3, cut_off);
    }
}

int factorial(int x, int result = 1) {
    if (x == 1||x==0) return result; else return factorial(x - 1, x * result);
}


int** combination_fixed_total(int tot, int groups){
    string s;
    int tot_permu=factorial(tot+groups-1)/(factorial(tot)*factorial(groups-1));
    
    for (int i=0; i<groups-1; i++) {
        s+="0";
    }
    for (int i=0; i<tot; i++) {
        s+="1";
    }
    int **arr=new int*[tot_permu];
    for (int i=0; i<tot_permu; i++) {
        arr[i]=new int[tot+groups-1];
    }
    int i=0;
    do
    {
        //std::cout << s << '\n';
        for (int j=0; j<s.length(); j++) {
            int a=s[j]-48;
            arr[i][j]=a;
        }
        i++;
    }
    while (std::next_permutation(s.begin(), s.end()));
    return arr;
}
void swap(int *fir, int *sec)
{
    int temp = *fir;
    *fir = *sec;
    *sec = temp;
}
static int i=0;
void permutation(int **arr2,int * arr, int curr, int size,  int tot)
{
    if(curr == size-1)
    {
        for(int a=0; a<size; a++)
            arr2[i][a]=arr[a];
        i++;
        if (i==tot) {
            i=0;
        }
    }
    
    else
    {
        for(int i=curr; i<size; i++)
        {
            swap(&arr[curr], &arr[i]);
            permutation(arr2,arr, curr+1, size,tot);
            swap(&arr[curr], &arr[i]);
        }
    }
}
int **permu_fun(int no_of_ver){
    int tot_permu=factorial(no_of_ver);
    int**arr2=new int*[tot_permu];
    for (int i=0; i<tot_permu; i++) {
        arr2[i]=new int[no_of_ver];
    }
    int *arr=new int[no_of_ver];
    for (int i=0; i<no_of_ver; i++) {
        arr[i]=i;
    }
    permutation(arr2, arr, 0, no_of_ver, tot_permu);
    delete [] arr;
    return arr2;
}
void symop_DtoC(float*** linear_tran, float*** linear_tran_c, int no_of_sym, Eigen::MatrixXd trans){
    for (int idx=0; idx<no_of_sym; idx++) {
        Eigen::MatrixXd temp(3,3);
        for (int i=0; i<3; i++) {
            for (int j=0;j<3; j++) {
                temp(i,j)=linear_tran[idx][i][j];
            }
        }
        temp=trans*temp*trans.inverse();
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                linear_tran_c[idx][i][j]=temp(i,j);
            }
        }
    }
}

#endif
