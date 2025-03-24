#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>

#include "huffman_basic.h"

#define Pi    3.14159265358979323846


double com_unit(int N,int index)
{
double aphla = 0.00;
if(index==0)
{
aphla=sqrt((double)1/(double)N);
}
else 
{
aphla=sqrt((double)2/(double)N);
}
return aphla;
}


//one dim
std::vector<double> dct_compute(std::vector<double> &data_)
{
int length_list=data_.size();
std::vector<double> results;
for(int index_=0;index_<length_list;index_++)
{
double temp_v =0.000;
double aphla_=com_unit(length_list,index_);
std::cout<<"aphla :"<<aphla_<<std::endl;
    for(int index=0;index<length_list;index++)
    {   
        //std::cout<<"data_ "<<data_[index]<<std::endl;
        temp_v+=data_[index]*cos(((2.0*index+1)*index_*Pi)/(2.0*length_list));
        //temp_v+=data_[index]*cos((index+0.5)*Pi*index_/length_list);
        //std::cout<<"temp_v :"<<temp_v<<std::endl;
    }
    results.push_back(temp_v*aphla_);
}
return results;
}

std::vector<double> inver_dct_compute(std::vector<double> &data_)
{
int length_list=data_.size();
std::vector<double> results;
for(int index_=0;index_<length_list;index_++)
{
double temp_v =0.000;

//std::cout<<"aphla :"<<aphla_<<std::endl;
    for(int index=0;index<length_list;index++)
    {   
        double aphla_=com_unit(length_list,index);
        //std::cout<<"data_ "<<data_[index]<<std::endl;
        temp_v+=aphla_*data_[index]*cos(((2.0*index_+1)*index*Pi)/(2.0*length_list));
        //temp_v+=data_[index]*cos((index+0.5)*Pi*index_/length_list);
        
    }
    results.push_back(temp_v);
    //std::cout<<"temp_v :"<<temp_v<<std::endl;
}
return results;
}

template <typename T>
void show(std::vector<T> values)
{
    std::cout<<"========================"<<std::endl;
    for(auto single_value:values)
    {
        std::cout<<single_value<<" ";
    }
    std::cout<<"\n";
    std::cout<<"========================"<<std::endl;
}



double fun_C(int N,int index_)
{
    double rs=0.0;
    if(index_==0)
    {
        rs=sqrt(double(1.0)/double(N));
    }
    else
    {
        rs=sqrt(double(2.0)/double(N));
    }

    return rs;
}

//
std::vector<std::vector<double>> dct2dimensions(std::vector<std::vector<double>> matrix_)
{

    int max_N =0;
    max_N=std::max(matrix_.size(),matrix_[0].size());
    std::vector<std::vector<double>> res_m;
    for(int row=0;row<matrix_.size();row++)
    {   
        double c_r=fun_C(max_N,row);
        std::vector<double> temp_row;
        for(int col=0;col<matrix_[row].size();col++)
        {   
            //
            
            double c_c=fun_C(max_N,col);
            double temp_sum = 0.00;
            for(int row_i=0;row_i<matrix_.size();row_i++)
            {
            for(int col_i=0;col_i<matrix_[row].size();col_i++)
            {
            temp_sum += matrix_[row_i][col_i]*cos(((2*row_i+1)*Pi*row)/(2*max_N))*cos(((2*col_i+1)*Pi*col)/(2*max_N));
            }
            }
            //x(k,l) 
            temp_row.push_back(temp_sum*c_c*c_r);
        }
        res_m.push_back(temp_row);
    }
    return res_m;
}




std::vector<std::vector<double>> in_dct2dimensions(std::vector<std::vector<double>> matrix_)
{

    int max_N =0;
    max_N=std::max(matrix_.size(),matrix_[0].size());
    std::vector<std::vector<double>> res_m;
    for(int row=0;row<matrix_.size();row++)
    {   
        
        std::vector<double> temp_row;
        for(int col=0;col<matrix_[row].size();col++)
        {   
            //
            
            
            double temp_sum = 0.00;
            for(int row_i=0;row_i<matrix_.size();row_i++)
            {
            double c_r=fun_C(max_N,row_i);
            for(int col_i=0;col_i<matrix_[row].size();col_i++)
            {
            double c_c=fun_C(max_N,col_i);
            temp_sum +=c_c*c_r* matrix_[row_i][col_i]*cos(((2*row+1)*Pi*row_i)/(2*max_N))*cos(((2*col+1)*Pi*col_i)/(2*max_N));
            }
            }
            //x(k,l) 
            temp_row.push_back(temp_sum);
        }
        res_m.push_back(temp_row);
    }
    return res_m;
}







std::vector<std::vector<double>> exchange_vector_to_square_maxrix(std::vector<std::vector<double>> matrix_)
{
    int rows=matrix_.size();
    int cols=matrix_[0].size();
    int real_n=std::max(rows,cols);
    if(rows>cols)
    {
        //add cols
        for(int index_row=0;index_row<matrix_.size();index_row++)
        {
        for(int index_f=0;index_f<(rows-cols);index_f++)
        {
            matrix_[index_row].push_back(0.0);
        }
        }
    }
    else if(rows<cols) 
    {
        //add rows
        std::vector<double> add_row(cols,0.0);
        for(int index_o=0;index_o<(cols-rows);index_o++)
        {
            matrix_.push_back(add_row);
        }
    }
    else 
    {
        //do nothing
    }
    return matrix_;
}

void show_2dimeon(std::vector<std::vector<double>> values_d)
{
    //
    std::cout<<"====================="<<std::endl;
    for(int row_=0;row_<values_d.size();row_++)
    {
        for(int col_=0;col_<values_d.size();col_++)
        {
            std::cout<<values_d[row_][col_]<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"====================="<<std::endl;
    //
}
















int main()
{
    std::vector<double> origin_ = {104.0,101.0,108.0,108.0,111.0};
    std::vector<double> after_dct=dct_compute(origin_);
    std::vector<double> inver_after=inver_dct_compute(after_dct);
    show(after_dct);
    show(inver_after);

    std::vector<std::vector<double>> temps_={{0.5961955,0.6660101,0.33207228},{0.06372824,0.99131946,0.74748041},{0.59761423,0.28287965,0.90013419}};
    show_2dimeon(temps_);

    std::vector<std::vector<double>> matrix_ = {{1,2,3},{4,5,6}};
    std::vector<std::vector<double>> matrix_reduce = {{1,2},{3,4},{5,6}};
    std::vector<std::vector<double>> matrix_after=exchange_vector_to_square_maxrix(matrix_);
    std::vector<std::vector<double>> matrix_after_=exchange_vector_to_square_maxrix(matrix_reduce);
    show_2dimeon(matrix_after);
    show_2dimeon(matrix_after_);

    std::vector<std::vector<double>> dct2_after=dct2dimensions(temps_);
    std::vector<std::vector<double>> origin_dct2=in_dct2dimensions(dct2_after);
    show_2dimeon(dct2_after);
    show_2dimeon(origin_dct2);


    

}