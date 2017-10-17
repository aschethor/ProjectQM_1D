#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
using namespace std;
// utility (random number between 0 and 1)
double drand(){
    return ((double) rand() / (RAND_MAX));
}

//even more utility
vector<double> add(vector<double>& a,vector<double>& b){
    vector<double> ret(a.size());
    for(int i=0;i<a.size();i++){
        ret[i]=a[i]+b[i];
    }
    return ret;
}

vector<double> multiply(vector<double> a,double b){
    vector<double> ret(a.size());
    for(int i=0;i<a.size();i++){
        ret[i]=a[i]*b;
    }
    return ret;
}

vector<double> multiply(vector<double> a,vector<double> b){
    vector<double> ret(a.size());
    for(int i=0;i<a.size();i++){
        ret[i]=a[i]*b[i];
    }
    return ret;
}

double sum(vector<double> a){
    double ret = 0;
    for(int i=0;i<a.size();i++){
        ret+=a[i];
    }
    return ret;
}

double avrg(vector<double> a){
    return sum(a)/a.size();
}

void print(vector<double> v){
    cout<<"[ ";
    for(int i=0;i<v.size();i++){
        cout<<v[i]<<" ";
    }
    cout<<"]";
}